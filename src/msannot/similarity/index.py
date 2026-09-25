"""Index vectorisé : compare une requête à toute une bibliothèque en une passe numpy.

Principe :

1. Tous les pics de la bibliothèque sont concaténés et triés par m/z. Pour chaque pic de
   la requête, ``searchsorted`` donne d'un coup tous les pics de la bibliothèque à moins
   de la tolérance.
2. Pour le **modified cosine**, la condition « fragment décalé de la différence de
   précurseurs » équivaut à l'égalité des **pertes neutres** (précurseur − fragment) :
   ``|(P_ref − m_ref) − (P_req − m_req)| ≤ tolérance``. Un second index, trié par perte
   neutre, les trouve de la même façon, alors que le décalage varie d'un spectre à l'autre.
3. Chaque paire reçoit sa contribution : produit des poids pour le cosine, terme
   d'entropie pour l'entropie spectrale (voir :mod:`msannot.similarity.weights`).
4. Pour un spectre de la bibliothèque, si aucun pic n'est impliqué dans deux paires, le
   score est la simple somme des contributions (``bincount``). Sinon, l'appariement
   glouton exact de :mod:`msannot.similarity.pairwise` est appliqué à ce spectre seulement.

Le résultat est identique aux fonctions de référence paire par paire, à l'arrondi près.
Ces fonctions sont elles-mêmes validées contre matchms et ms_entropy par les tests.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from msannot.config import Metric, SimilarityConfig
from msannot.models import Spectrum
from msannot.similarity.pairwise import greedy_select
from msannot.similarity.weights import (
    LN4,
    cosine_weights,
    entropy_pair_contribution,
    entropy_weights,
)

_WINDOW_MARGIN = 1e-6  # Da ; très supérieur aux erreurs d'arrondi flottant


@dataclass(frozen=True, slots=True)
class QueryScores:
    """Scores d'une requête contre la bibliothèque (spectres ayant au moins une paire)."""

    indices: np.ndarray
    scores: np.ndarray
    n_matches: np.ndarray

    def dense(self, size: int) -> tuple[np.ndarray, np.ndarray]:
        """Scores et nombres de pics appariés pour tous les spectres (0 si aucune paire)."""
        scores = np.zeros(size)
        matches = np.zeros(size, dtype=np.int64)
        scores[self.indices] = self.scores
        matches[self.indices] = self.n_matches
        return scores, matches


def _window_pairs(
    sorted_values: np.ndarray, targets: np.ndarray, tolerance: float
) -> tuple[np.ndarray, np.ndarray]:
    """Paires (cible, position triée) telles que ``|valeur − cible| ≤ tolérance``."""
    low = np.searchsorted(sorted_values, targets - tolerance, side="left")
    high = np.searchsorted(sorted_values, targets + tolerance, side="right")
    counts = np.maximum(high - low, 0)
    target_idx = np.repeat(np.arange(targets.size), counts)
    starts = np.repeat(low - np.cumsum(counts) + counts, counts)
    positions = np.arange(int(counts.sum())) + starts
    return target_idx, positions


class LibraryIndex:
    """Index de recherche sur une liste de spectres (déjà nettoyés).

    Parameters
    ----------
    spectra
        Spectres de la bibliothèque.
    config
        Tolérance et transformations d'intensité.
    """

    def __init__(self, spectra: Sequence[Spectrum], config: SimilarityConfig) -> None:
        if not spectra:
            raise ValueError("bibliothèque vide")
        self.config = config
        self.size = len(spectra)
        self.precursors = np.array([s.precursor_mz for s in spectra], dtype=np.float64)
        counts = np.array([s.n_peaks for s in spectra], dtype=np.int64)
        self.peak_spectrum = np.repeat(np.arange(self.size), counts)
        self.peak_local = np.concatenate([np.arange(n) for n in counts])
        mz = np.concatenate([s.mz for s in spectra])
        self.cosine_weights = np.concatenate(
            [cosine_weights(s.intensities, config.intensity_power) for s in spectra]
        )
        self.entropy_weights = np.concatenate(
            [entropy_weights(s.intensities, config.entropy_weighted) for s in spectra]
        )
        self._all_mz = mz
        self._fragment_order = np.argsort(mz, kind="stable")
        self._fragment_values = mz[self._fragment_order]
        losses = self.precursors[self.peak_spectrum] - mz
        self._loss_order = np.argsort(losses, kind="stable")
        self._loss_values = losses[self._loss_order]

    # ------------------------------------------------------------------ paires
    def _pairs(
        self, query: Spectrum, metric: Metric, candidates: np.ndarray | None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Paires (pic requête, pic bibliothèque global, décalé ?) après filtrage.

        Les m/z publiés sont arrondis (souvent à 4 décimales) : des écarts exactement égaux à
        la tolérance sont fréquents. La recherche se fait donc dans une fenêtre à peine
        élargie, puis chaque paire est revérifiée avec **exactement** les mêmes opérations
        flottantes que la référence (``m/z_req + décalage`` comparé à
        ``m/z_ref ± tolérance``), pour des résultats identiques aux cas limites près.
        """
        tolerance = self.config.tolerance
        window = tolerance + _WINDOW_MARGIN
        library_mz = self._all_mz
        q_idx, positions = _window_pairs(self._fragment_values, query.mz, window)
        peaks = self._fragment_order[positions]
        if candidates is not None:  # filtrage précoce : peu de candidats en identification
            keep = candidates[self.peak_spectrum[peaks]]
            q_idx, peaks = q_idx[keep], peaks[keep]
        values = query.mz[q_idx] + 0.0
        exact = (values >= library_mz[peaks] - tolerance) & (
            values <= library_mz[peaks] + tolerance
        )
        q_idx, peaks = q_idx[exact], peaks[exact]
        shifted = np.zeros(q_idx.size, dtype=bool)
        if metric == "modified_cosine":
            losses = query.precursor_mz - query.mz
            l_idx, l_positions = _window_pairs(self._loss_values, losses, window)
            l_peaks = self._loss_order[l_positions]
            if candidates is not None:
                keep = candidates[self.peak_spectrum[l_peaks]]
                l_idx, l_peaks = l_idx[keep], l_peaks[keep]
            shift = self.precursors[self.peak_spectrum[l_peaks]] - query.precursor_mz
            values = query.mz[l_idx] + shift
            keep = (
                (values >= library_mz[l_peaks] - tolerance)
                & (values <= library_mz[l_peaks] + tolerance)
                # Comme matchms : pas d'appariement décalé si les précurseurs sont confondus.
                & (np.abs(shift) > tolerance)
            )
            q_idx = np.concatenate((q_idx, l_idx[keep]))
            peaks = np.concatenate((peaks, l_peaks[keep]))
            shifted = np.concatenate((shifted, np.ones(int(keep.sum()), dtype=bool)))
            # Une même paire trouvée directement et par décalage ne compte qu'une fois.
            key = q_idx.astype(np.int64) * (self.peak_spectrum.size + 1) + peaks
            _, first = np.unique(key, return_index=True)
            first.sort()
            q_idx, peaks, shifted = q_idx[first], peaks[first], shifted[first]
        return q_idx, peaks, shifted

    # ------------------------------------------------------------------ requête
    def query(
        self, query: Spectrum, metric: Metric, candidates: np.ndarray | None = None
    ) -> QueryScores:
        """Scores de ``query`` contre les spectres de la bibliothèque.

        Parameters
        ----------
        query
            Spectre requête, nettoyé comme la bibliothèque.
        metric
            ``"cosine"``, ``"modified_cosine"`` ou ``"entropy"``.
        candidates
            Masque booléen des spectres autorisés (tous par défaut).
        """
        q_idx, peaks, shifted = self._pairs(query, metric, candidates)
        if q_idx.size == 0:
            empty = np.zeros(0, dtype=np.int64)
            return QueryScores(empty, np.zeros(0), empty)
        if metric == "entropy":
            q_weights = entropy_weights(query.intensities, self.config.entropy_weighted)
            contributions = entropy_pair_contribution(q_weights[q_idx], self.entropy_weights[peaks])
        else:
            q_weights = cosine_weights(query.intensities, self.config.intensity_power)
            contributions = q_weights[q_idx] * self.cosine_weights[peaks]
        spectrum = self.peak_spectrum[peaks]

        # Spectres « en conflit » : un pic de la requête ou de la bibliothèque est apparié
        # plusieurs fois ; l'appariement glouton exact est nécessaire pour eux seuls.
        query_key = spectrum.astype(np.int64) * (query.n_peaks + 1) + q_idx
        conflict = np.zeros(self.size, dtype=bool)
        for key in (query_key, peaks):
            values, counts = np.unique(key, return_counts=True)
            duplicated = np.isin(key, values[counts > 1])
            conflict[spectrum[duplicated]] = True

        simple = ~conflict[spectrum]
        # astype : bincount d'une liste vide renvoie des entiers (scores tronqués sinon)
        scores = np.bincount(
            spectrum[simple], weights=contributions[simple], minlength=self.size
        ).astype(np.float64)
        matches = np.bincount(spectrum[simple], minlength=self.size)
        for index in np.flatnonzero(conflict):
            rows = np.flatnonzero(spectrum == index)
            # Ordre canonique de matchms : paires directes puis décalées, par pic de
            # référence puis pic de requête, avant le tri glouton.
            order = np.lexsort((q_idx[rows], self.peak_local[peaks[rows]], shifted[rows]))
            rows = rows[order]
            selected = greedy_select(q_idx[rows], peaks[rows], contributions[rows])
            scores[index] = contributions[rows][selected].sum()
            matches[index] = int(selected.sum())

        present = np.flatnonzero(matches > 0)
        values = scores[present]
        if metric == "entropy":
            values = np.minimum(values / LN4, 1.0)
        return QueryScores(present, values, matches[present])
