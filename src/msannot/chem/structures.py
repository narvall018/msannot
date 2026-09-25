"""Structures chimiques avec RDKit : validation, masses, empreintes et similarité.

- L'**InChIKey** identifie une structure ; son premier bloc (14 caractères) encode la
  connectivité sans la stéréochimie. Deux spectres de stéréo-isomères, indiscernables en
  MS/MS classique, sont donc considérés comme le même composé.
- Les **empreintes de Morgan** (rayon 2, équivalentes aux ECFP4) codent les environnements
  atomiques ; la **similarité de Tanimoto** entre empreintes est la mesure standard de
  similarité structurale.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

# Avertissements RDKit coupés : nos propres messages signalent les structures invalides.
RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]  # absent des annotations RDKit

PROTON_MASS = 1.007276466812


@dataclass(frozen=True, slots=True)
class StructureInfo:
    """Propriétés d'une structure calculées par RDKit."""

    canonical_smiles: str
    inchikey: str
    formula: str
    exact_mass: float
    formal_charge: int

    @property
    def inchikey14(self) -> str:
        """Premier bloc de l'InChIKey (connectivité)."""
        return self.inchikey[:14]


@lru_cache(maxsize=65536)
def parse_structure(smiles: str) -> StructureInfo | None:
    """Analyse un SMILES ; ``None`` s'il est invalide ou sans InChIKey calculable."""
    if not smiles:
        return None
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        return None
    inchikey = Chem.MolToInchiKey(molecule)
    if not inchikey:
        return None
    return StructureInfo(
        canonical_smiles=Chem.MolToSmiles(molecule),
        inchikey=inchikey,
        formula=rdMolDescriptors.CalcMolFormula(molecule),
        exact_mass=float(rdMolDescriptors.CalcExactMolWt(molecule)),
        formal_charge=int(Chem.GetFormalCharge(molecule)),
    )


def protonated_mz(info: StructureInfo) -> float:
    """m/z théorique de l'ion [M+H]+ d'une molécule neutre."""
    return info.exact_mass + PROTON_MASS


def ppm_error(observed: float, expected: float) -> float:
    """Écart relatif en parties par million."""
    return (observed - expected) / expected * 1e6


class FingerprintIndex:
    """Empreintes de Morgan indexées par premier bloc d'InChIKey.

    Parameters
    ----------
    smiles_by_key
        Structure SMILES de chaque composé (clé : premier bloc d'InChIKey).
    radius, n_bits
        Paramètres des empreintes de Morgan.
    """

    def __init__(self, smiles_by_key: dict[str, str], radius: int = 2, n_bits: int = 2048) -> None:
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
        self.keys: list[str] = []
        self.fingerprints: list[DataStructs.ExplicitBitVect] = []
        for key, smiles in sorted(smiles_by_key.items()):
            molecule = Chem.MolFromSmiles(smiles)
            if molecule is None:
                continue
            self.keys.append(key)
            self.fingerprints.append(generator.GetFingerprint(molecule))
        self.position = {key: index for index, key in enumerate(self.keys)}

    def __contains__(self, key: object) -> bool:
        return key in self.position

    def __len__(self) -> int:
        return len(self.keys)

    def similarity(self, key_a: str, key_b: str) -> float:
        """Tanimoto entre deux composés."""
        return float(
            DataStructs.TanimotoSimilarity(
                self.fingerprints[self.position[key_a]], self.fingerprints[self.position[key_b]]
            )
        )

    def similarity_to_all(self, key: str) -> np.ndarray:
        """Tanimoto d'un composé avec tous les composés de l'index (ordre de ``keys``)."""
        values = DataStructs.BulkTanimotoSimilarity(
            self.fingerprints[self.position[key]], self.fingerprints
        )
        return np.asarray(values, dtype=np.float64)

    def indices(self, keys: Iterable[str]) -> np.ndarray:
        """Positions de composés dans l'index."""
        return np.array([self.position[key] for key in keys], dtype=np.int64)


def molecule_png(smiles: str, width: int = 300, height: int = 220) -> bytes | None:
    """Dessin 2D d'une molécule (PNG) ; ``None`` si le SMILES est invalide."""
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()
    return bytes(drawer.GetDrawingText())


def molecule_svg(smiles: str, width: int = 300, height: int = 220) -> str | None:
    """Dessin 2D d'une molécule (SVG) ; ``None`` si le SMILES est invalide."""
    molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return None
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, molecule)
    drawer.FinishDrawing()
    return str(drawer.GetDrawingText())


def unique_structures(pairs: Sequence[tuple[str, str]]) -> dict[str, str]:
    """Premier SMILES rencontré pour chaque premier bloc d'InChIKey."""
    structures: dict[str, str] = {}
    for key, smiles in pairs:
        if key and smiles:
            structures.setdefault(key, smiles)
    return structures
