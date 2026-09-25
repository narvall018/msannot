"""Transformations d'intensité propres à chaque mesure de similarité.

- **Cosine** : intensités élevées à la puissance ``p`` (0,5 par défaut, ce qui atténue la
  domination des pics majoritaires), puis normalisées en norme euclidienne.
- **Entropie** (Li et al., 2021) : intensités normalisées en somme 1. Si l'entropie
  spectrale ``S`` est inférieure à 3, elles sont repondérées par ``I^(0,25 + 0,25·S)``, ce qui
  renforce les petits pics des spectres pauvres, puis renormalisées.
"""

from __future__ import annotations

import numpy as np

ENTROPY_WEIGHT_CUTOFF = 3.0


def cosine_weights(intensities: np.ndarray, power: float) -> np.ndarray:
    """Intensités transformées pour le cosine (norme euclidienne = 1)."""
    powered = np.power(np.asarray(intensities, dtype=np.float64), power)
    norm = float(np.sqrt(np.sum(powered * powered)))
    return powered / norm if norm > 0 else powered


def spectral_entropy(intensities: np.ndarray) -> float:
    """Entropie de Shannon (en nats) de la distribution des intensités."""
    total = float(np.sum(intensities))
    if total <= 0:
        return 0.0
    probabilities = np.asarray(intensities, dtype=np.float64) / total
    probabilities = probabilities[probabilities > 0]
    return float(-np.sum(probabilities * np.log(probabilities)))


def entropy_weights(intensities: np.ndarray, weighted: bool = True) -> np.ndarray:
    """Intensités transformées pour la similarité d'entropie (somme = 1)."""
    values = np.asarray(intensities, dtype=np.float64)
    total = float(values.sum())
    if total <= 0:
        return values
    probabilities = values / total
    if weighted:
        entropy = spectral_entropy(probabilities)
        if entropy < ENTROPY_WEIGHT_CUTOFF:
            probabilities = np.power(probabilities, 0.25 + 0.25 * entropy)
            probabilities = probabilities / probabilities.sum()
    return probabilities


def entropy_pair_contribution(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Contribution d'une paire de pics appariés à la similarité d'entropie (numérateur).

    La similarité d'entropie ``1 − (2·S_AB − S_A − S_B) / ln 4`` se réécrit comme une somme
    sur les paires appariées de ``(a+b)·ln(a+b) − a·ln a − b·ln b``, divisée par ``ln 4``
    (décomposition utilisée par Flash entropy, Li & Fiehn, 2023) : les pics non appariés ne
    contribuent pas. C'est ce qui permet une recherche vectorisée.
    """
    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)
    total = a + b
    with np.errstate(divide="ignore", invalid="ignore"):
        value = total * np.log(total) - a * np.log(a) - b * np.log(b)
    return np.where((a > 0) & (b > 0), value, 0.0)


LN4 = float(np.log(4.0))
