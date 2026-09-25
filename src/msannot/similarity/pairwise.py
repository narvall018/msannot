"""Similarités entre deux spectres : implémentations de référence, lisibles et testées.

Conventions identiques à matchms (``CosineGreedy``, ``ModifiedCosineGreedy``) :

- deux pics correspondent si ``|m/z_a − (m/z_b + décalage)| ≤ tolérance`` ;
- chaque pic n'est apparié qu'une fois ; les paires sont choisies de façon gloutonne par
  contribution décroissante (en cas d'égalité, même ordre que matchms) ;
- le cosine est normalisé par la norme de **tous** les pics, appariés ou non.

Le **modified cosine** apparie aussi les pics décalés de la différence de masse des
précurseurs : un fragment qui porte la modification chimique d'un analogue est retrouvé.
Ces fonctions servent de référence au moteur vectorisé (:mod:`msannot.similarity.index`),
qui doit donner les mêmes scores.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from msannot.config import Metric, SimilarityConfig
from msannot.models import Spectrum
from msannot.similarity.weights import (
    LN4,
    cosine_weights,
    entropy_pair_contribution,
    entropy_weights,
)


@dataclass(frozen=True, slots=True)
class SimilarityResult:
    """Score de similarité et nombre de pics appariés."""

    score: float
    n_matches: int


def find_matches(
    mz_a: np.ndarray, mz_b: np.ndarray, tolerance: float, shift: float = 0.0
) -> tuple[np.ndarray, np.ndarray]:
    """Toutes les paires ``(i, j)`` avec ``|mz_a[i] − (mz_b[j] + shift)| ≤ tolérance``.

    Les paires sont ordonnées par ``i`` puis ``j`` croissants (ordre de matchms).
    """
    shifted = mz_b + shift
    low = np.searchsorted(shifted, mz_a - tolerance, side="left")
    high = np.searchsorted(shifted, mz_a + tolerance, side="right")
    counts = np.maximum(high - low, 0)
    idx_a = np.repeat(np.arange(mz_a.size), counts)
    starts = np.repeat(low - np.cumsum(counts) + counts, counts)
    idx_b = np.arange(int(counts.sum())) + starts
    return idx_a, idx_b


def greedy_select(idx_a: np.ndarray, idx_b: np.ndarray, scores: np.ndarray) -> np.ndarray:
    """Sélection gloutonne un-pour-un ; renvoie le masque des paires retenues."""
    order = np.argsort(scores, kind="mergesort")[::-1]
    used_a: set[int] = set()
    used_b: set[int] = set()
    selected = np.zeros(scores.size, dtype=bool)
    for position in order:
        a, b = int(idx_a[position]), int(idx_b[position])
        if a in used_a or b in used_b:
            continue
        used_a.add(a)
        used_b.add(b)
        selected[position] = True
    return selected


def _cosine_from_pairs(
    idx_a: np.ndarray, idx_b: np.ndarray, weights_a: np.ndarray, weights_b: np.ndarray
) -> SimilarityResult:
    if idx_a.size == 0:
        return SimilarityResult(0.0, 0)
    products = weights_a[idx_a] * weights_b[idx_b]
    selected = greedy_select(idx_a, idx_b, products)
    return SimilarityResult(float(products[selected].sum()), int(selected.sum()))


def cosine(a: Spectrum, b: Spectrum, tolerance: float, intensity_power: float) -> SimilarityResult:
    """Cosine glouton entre deux spectres."""
    idx_a, idx_b = find_matches(a.mz, b.mz, tolerance)
    return _cosine_from_pairs(
        idx_a,
        idx_b,
        cosine_weights(a.intensities, intensity_power),
        cosine_weights(b.intensities, intensity_power),
    )


def modified_cosine(
    a: Spectrum, b: Spectrum, tolerance: float, intensity_power: float
) -> SimilarityResult:
    """Modified cosine : appariement direct ou décalé de ``précurseur_a − précurseur_b``.

    Si la différence de précurseur est inférieure à la tolérance, c'est un cosine simple.
    """
    shift = a.precursor_mz - b.precursor_mz
    if abs(shift) <= tolerance:
        return cosine(a, b, tolerance, intensity_power)
    zero_a, zero_b = find_matches(a.mz, b.mz, tolerance)
    shift_a, shift_b = find_matches(a.mz, b.mz, tolerance, shift=shift)
    return _cosine_from_pairs(
        np.concatenate((zero_a, shift_a)),
        np.concatenate((zero_b, shift_b)),
        cosine_weights(a.intensities, intensity_power),
        cosine_weights(b.intensities, intensity_power),
    )


def entropy_similarity(
    a: Spectrum, b: Spectrum, tolerance: float, weighted: bool = True
) -> SimilarityResult:
    """Similarité d'entropie spectrale (pondérée par défaut), entre 0 et 1."""
    idx_a, idx_b = find_matches(a.mz, b.mz, tolerance)
    if idx_a.size == 0:
        return SimilarityResult(0.0, 0)
    weights_a = entropy_weights(a.intensities, weighted)
    weights_b = entropy_weights(b.intensities, weighted)
    contributions = entropy_pair_contribution(weights_a[idx_a], weights_b[idx_b])
    selected = greedy_select(idx_a, idx_b, contributions)
    score = float(contributions[selected].sum()) / LN4
    return SimilarityResult(min(score, 1.0), int(selected.sum()))


@dataclass(frozen=True, slots=True)
class PeakMatches:
    """Pics appariés entre un spectre de référence ``a`` et une requête ``b``."""

    idx_a: np.ndarray
    idx_b: np.ndarray
    shifted: np.ndarray


def matched_peaks(
    metric: Metric, a: Spectrum, b: Spectrum, config: SimilarityConfig
) -> PeakMatches:
    """Paires de pics retenues par la mesure (pour les spectres miroirs)."""
    tolerance = config.tolerance
    idx_a, idx_b = find_matches(a.mz, b.mz, tolerance)
    shifted = np.zeros(idx_a.size, dtype=bool)
    if metric == "modified_cosine" and abs(a.precursor_mz - b.precursor_mz) > tolerance:
        shift_a, shift_b = find_matches(a.mz, b.mz, tolerance, a.precursor_mz - b.precursor_mz)
        idx_a = np.concatenate((idx_a, shift_a))
        idx_b = np.concatenate((idx_b, shift_b))
        shifted = np.concatenate((shifted, np.ones(shift_a.size, dtype=bool)))
    if idx_a.size == 0:
        return PeakMatches(idx_a, idx_b, shifted)
    if metric == "entropy":
        weights_a = entropy_weights(a.intensities, config.entropy_weighted)
        weights_b = entropy_weights(b.intensities, config.entropy_weighted)
        scores = entropy_pair_contribution(weights_a[idx_a], weights_b[idx_b])
    else:
        scores = (
            cosine_weights(a.intensities, config.intensity_power)[idx_a]
            * cosine_weights(b.intensities, config.intensity_power)[idx_b]
        )
    selected = greedy_select(idx_a, idx_b, scores)
    return PeakMatches(idx_a[selected], idx_b[selected], shifted[selected])


def similarity(
    metric: Metric, a: Spectrum, b: Spectrum, config: SimilarityConfig
) -> SimilarityResult:
    """Calcule la similarité ``metric`` entre deux spectres."""
    if metric == "cosine":
        return cosine(a, b, config.tolerance, config.intensity_power)
    if metric == "modified_cosine":
        return modified_cosine(a, b, config.tolerance, config.intensity_power)
    return entropy_similarity(a, b, config.tolerance, config.entropy_weighted)
