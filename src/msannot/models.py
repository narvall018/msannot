"""Structure de données centrale : un spectre MS/MS et ses métadonnées.

Un spectre est une liste de pics (m/z, intensité) triés par m/z croissant, accompagnée de
la masse du précurseur et de métadonnées textuelles normalisées (clés en minuscules :
``name``, ``smiles``, ``inchikey``, ``contributor``, ``license``…).
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType

import numpy as np

UNKNOWN = "inconnu"


def _readonly(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


@dataclass(frozen=True, slots=True, eq=False)
class Spectrum:
    """Spectre MS/MS centroïdé.

    Attributes
    ----------
    identifier
        Identifiant unique (accession MassBank, titre MGF...).
    mz, intensities
        Pics, triés par m/z croissant (tableaux en lecture seule).
    precursor_mz
        m/z de l'ion précurseur.
    metadata
        Métadonnées normalisées (lecture seule).
    """

    identifier: str
    mz: np.ndarray
    intensities: np.ndarray
    precursor_mz: float
    metadata: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        mz = np.asarray(self.mz, dtype=np.float64).copy()
        intensities = np.asarray(self.intensities, dtype=np.float64).copy()
        if mz.ndim != 1 or mz.shape != intensities.shape:
            raise ValueError(f"{self.identifier} : m/z et intensités de tailles différentes")
        if not (np.all(np.isfinite(mz)) and np.all(np.isfinite(intensities))):
            raise ValueError(f"{self.identifier} : valeurs non finies dans les pics")
        if np.any(intensities < 0) or np.any(mz <= 0):
            raise ValueError(f"{self.identifier} : m/z ou intensité négatif")
        if not (np.isfinite(self.precursor_mz) and self.precursor_mz > 0):
            raise ValueError(f"{self.identifier} : m/z du précurseur invalide")
        if mz.size > 1 and np.any(np.diff(mz) < 0):
            order = np.argsort(mz, kind="stable")
            mz, intensities = mz[order], intensities[order]
        object.__setattr__(self, "mz", _readonly(mz))
        object.__setattr__(self, "intensities", _readonly(intensities))
        object.__setattr__(self, "precursor_mz", float(self.precursor_mz))
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def n_peaks(self) -> int:
        """Nombre de pics."""
        return int(self.mz.size)

    def get(self, key: str, default: str = "") -> str:
        """Métadonnée ``key`` (chaîne vide si absente)."""
        return self.metadata.get(key, default)

    @property
    def name(self) -> str:
        """Nom du composé."""
        return self.get("name", self.identifier)

    @property
    def smiles(self) -> str:
        """Structure SMILES (vide si inconnue)."""
        return self.get("smiles")

    @property
    def inchikey14(self) -> str:
        """Premier bloc de l'InChIKey : squelette moléculaire, sans stéréochimie."""
        key = self.get("inchikey14") or self.get("inchikey")
        return key[:14]

    @property
    def contributor(self) -> str:
        """Laboratoire ou contributeur à l'origine du spectre."""
        return self.get("contributor", UNKNOWN) or UNKNOWN

    def with_peaks(self, mz: np.ndarray, intensities: np.ndarray) -> Spectrum:
        """Copie du spectre avec d'autres pics (mêmes métadonnées)."""
        return Spectrum(self.identifier, mz, intensities, self.precursor_mz, self.metadata)

    def with_metadata(self, **updates: str) -> Spectrum:
        """Copie du spectre avec des métadonnées ajoutées ou remplacées."""
        return Spectrum(
            self.identifier,
            self.mz,
            self.intensities,
            self.precursor_mz,
            {**self.metadata, **updates},
        )
