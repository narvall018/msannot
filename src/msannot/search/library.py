"""Bibliothèque spectrale et moteur de recherche.

Deux modes de recherche :

- **identification** : seuls les spectres dont le précurseur est à moins de ``ppm`` de celui
  de la requête sont candidats. On cherche le même composé ;
- **analogues** : tous les spectres (jusqu'à un écart de masse maximal) sont candidats, et on
  cherche des structures voisines. Le modified cosine est la mesure adaptée à ce mode.

Classement : score décroissant, puis nombre de pics appariés décroissant, puis ordre de la
bibliothèque (résultat déterministe).
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Literal

import numpy as np
import pandas as pd

from msannot.config import Metric, SimilarityConfig
from msannot.models import Spectrum
from msannot.similarity.index import LibraryIndex

SearchMode = Literal["identity", "analog"]
HIT_COLUMNS = [
    "rank",
    "library_index",
    "identifier",
    "name",
    "inchikey14",
    "smiles",
    "contributor",
    "precursor_mz",
    "delta_mz",
    "score",
    "n_matches",
]


def duplicate_groups(spectra: Sequence[Spectrum], decimals: int = 3) -> np.ndarray:
    """Identifiant de groupe des spectres strictement identiques (pics arrondis à 1e-3).

    Un doublon exact (même fichier déposé deux fois, parfois sous deux annotations) n'est
    pas une mesure indépendante : l'évaluation l'exclut des candidats de la requête.
    """
    groups: dict[tuple[float, bytes, bytes], int] = {}
    labels = np.empty(len(spectra), dtype=np.int64)
    for position, spectrum in enumerate(spectra):
        key = (
            round(spectrum.precursor_mz, decimals),
            np.round(spectrum.mz, decimals).tobytes(),
            np.round(spectrum.intensities, decimals).tobytes(),
        )
        labels[position] = groups.setdefault(key, len(groups))
    return labels


def rank_order(scores: np.ndarray, matches: np.ndarray) -> np.ndarray:
    """Indices triés par score décroissant, puis appariements décroissants, puis position."""
    positions = np.arange(scores.size)
    return np.lexsort((positions, -matches, -scores))


class SpectralLibrary:
    """Bibliothèque de spectres nettoyés, indexée pour la recherche.

    Parameters
    ----------
    spectra
        Spectres (déjà nettoyés, avec métadonnées ``inchikey14``, ``smiles``…).
    config
        Paramètres de similarité.
    """

    def __init__(self, spectra: Sequence[Spectrum], config: SimilarityConfig) -> None:
        self.spectra = list(spectra)
        self.config = config
        self.index = LibraryIndex(self.spectra, config)
        self.precursors = self.index.precursors
        self.table = pd.DataFrame(
            {
                "identifier": [s.identifier for s in self.spectra],
                "name": [s.name for s in self.spectra],
                "inchikey14": [s.inchikey14 for s in self.spectra],
                "smiles": [s.smiles for s in self.spectra],
                "contributor": [s.contributor for s in self.spectra],
                "instrument_type": [s.get("instrument_type") for s in self.spectra],
                "collision_energy": [s.get("collision_energy") for s in self.spectra],
                "license": [s.get("license") for s in self.spectra],
                "precursor_mz": self.precursors,
                "n_peaks": [s.n_peaks for s in self.spectra],
            }
        )
        self.inchikeys = self.table["inchikey14"].to_numpy()
        self.contributors = self.table["contributor"].to_numpy()
        self.duplicates = duplicate_groups(self.spectra)

    def duplicate_summary(self) -> dict[str, int]:
        """Bilan des spectres strictement identiques dans la bibliothèque."""
        frame = pd.DataFrame(
            {"group": self.duplicates, "compound": self.inchikeys, "lab": self.contributors}
        )
        sizes = frame.groupby("group").agg(
            n=("compound", "size"), compounds=("compound", "nunique"), labs=("lab", "nunique")
        )
        duplicated = sizes[sizes["n"] > 1]
        return {
            "groups": len(duplicated),
            "spectra": int(duplicated["n"].sum()),
            "cross_contributor_groups": int((duplicated["labs"] > 1).sum()),
            "conflicting_annotation_groups": int((duplicated["compounds"] > 1).sum()),
        }

    def __len__(self) -> int:
        return len(self.spectra)

    # ------------------------------------------------------------------ candidats
    def precursor_window(self, precursor_mz: float, ppm: float) -> np.ndarray:
        """Masque des spectres dont le précurseur est à moins de ``ppm`` de ``precursor_mz``."""
        return np.abs(self.precursors - precursor_mz) <= precursor_mz * ppm * 1e-6

    def mass_shift_window(self, precursor_mz: float, max_shift: float | None) -> np.ndarray:
        """Masque des spectres dont l'écart de précurseur est au plus ``max_shift`` (Da)."""
        if max_shift is None:
            return np.ones(len(self), dtype=bool)
        return np.abs(self.precursors - precursor_mz) <= max_shift

    # ------------------------------------------------------------------ scores
    def scores(
        self, query: Spectrum, metric: Metric, candidates: np.ndarray | None = None
    ) -> tuple[np.ndarray, np.ndarray]:
        """Scores (et pics appariés) de la requête contre toute la bibliothèque.

        Les spectres hors ``candidates`` reçoivent un score de 0.
        """
        return self.index.query(query, metric, candidates).dense(len(self))

    def search(
        self,
        query: Spectrum,
        *,
        metric: Metric,
        mode: SearchMode = "identity",
        top_k: int = 10,
        ppm: float = 10.0,
        max_mass_shift: float | None = 200.0,
        exclude: np.ndarray | None = None,
    ) -> pd.DataFrame:
        """Meilleurs résultats pour une requête.

        Parameters
        ----------
        query
            Spectre requête (nettoyé comme la bibliothèque).
        metric
            Mesure de similarité.
        mode
            ``"identity"`` (fenêtre de précurseur en ppm) ou ``"analog"``.
        top_k
            Nombre de résultats renvoyés (seuls les scores > 0 sont conservés).
        exclude
            Masque de spectres à exclure (par exemple le composé recherché lui-même).
        """
        if mode == "identity":
            candidates = self.precursor_window(query.precursor_mz, ppm)
        else:
            candidates = self.mass_shift_window(query.precursor_mz, max_mass_shift)
        if exclude is not None:
            candidates &= ~exclude
        scores, matches = self.scores(query, metric, candidates)
        order = [index for index in rank_order(scores, matches)[:top_k] if scores[index] > 0]
        hits = self.table.iloc[order][
            ["identifier", "name", "inchikey14", "smiles", "contributor", "precursor_mz"]
        ].copy()
        hits.insert(0, "library_index", order)
        hits.insert(0, "rank", np.arange(1, len(order) + 1))
        hits["delta_mz"] = hits["precursor_mz"] - query.precursor_mz
        hits["score"] = scores[order]
        hits["n_matches"] = matches[order]
        return hits.reset_index(drop=True)[HIT_COLUMNS]
