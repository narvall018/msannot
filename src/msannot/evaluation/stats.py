"""Outils statistiques utilisés par l'évaluation.

- Intervalle de **Wilson** pour les proportions (précis pour les petits effectifs et les
  proportions extrêmes).
- **Test de McNemar exact** pour comparer deux méthodes évaluées sur les **mêmes** requêtes :
  seules les requêtes où les deux méthodes divergent sont informatives (test binomial).
- **Test des rangs signés de Wilcoxon** pour des mesures continues appariées.
- **Intervalle bootstrap** (percentile) pour une médiane.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import stats

Z_95 = 1.959963984540054


def wilson_interval(successes: int, total: int, z: float = Z_95) -> tuple[float, float]:
    """Intervalle de confiance de Wilson (95 % par défaut) d'une proportion."""
    if total <= 0:
        return math.nan, math.nan
    if not 0 <= successes <= total:
        raise ValueError(f"effectifs incohérents : {successes}/{total}")
    p = successes / total
    denominator = 1.0 + z * z / total
    centre = (p + z * z / (2 * total)) / denominator
    half = z * math.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return max(0.0, min(p, centre - half)), min(1.0, max(p, centre + half))


@dataclass(frozen=True, slots=True)
class Proportion:
    """Proportion avec son intervalle de Wilson à 95 %."""

    successes: int
    total: int

    @property
    def value(self) -> float:
        """Valeur (NaN si effectif nul)."""
        return self.successes / self.total if self.total else math.nan

    def to_dict(self) -> dict[str, Any]:
        """Représentation sérialisable."""
        low, high = wilson_interval(self.successes, self.total)
        return {
            "value": self.value,
            "successes": self.successes,
            "total": self.total,
            "ci95_low": low,
            "ci95_high": high,
        }

    def __str__(self) -> str:
        if not self.total:
            return "n/a (0/0)"
        low, high = wilson_interval(self.successes, self.total)
        return f"{self.value:.3f} [{low:.3f}–{high:.3f}] ({self.successes}/{self.total})"


def mcnemar_exact(first_correct: np.ndarray, second_correct: np.ndarray) -> dict[str, float]:
    """Test de McNemar exact entre deux méthodes sur les mêmes requêtes.

    Returns
    -------
    dict
        ``only_first``, ``only_second`` (requêtes discordantes) et ``p_value`` bilatérale.
    """
    first = np.asarray(first_correct, dtype=bool)
    second = np.asarray(second_correct, dtype=bool)
    if first.shape != second.shape:
        raise ValueError("les deux méthodes doivent être évaluées sur les mêmes requêtes")
    only_first = int(np.sum(first & ~second))
    only_second = int(np.sum(~first & second))
    discordant = only_first + only_second
    p_value = (
        1.0
        if discordant == 0
        else float(stats.binomtest(min(only_first, only_second), discordant, 0.5).pvalue)
    )
    return {"only_first": only_first, "only_second": only_second, "p_value": p_value}


def wilcoxon_paired(first: np.ndarray, second: np.ndarray) -> dict[str, float]:
    """Test des rangs signés de Wilcoxon (bilatéral) et médiane des différences."""
    x = np.asarray(first, dtype=float)
    y = np.asarray(second, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    differences = x[keep] - y[keep]
    if differences.size == 0 or np.all(differences == 0):
        return {"n": int(differences.size), "median_difference": 0.0, "p_value": 1.0}
    result = stats.wilcoxon(x[keep], y[keep], zero_method="wilcox")
    return {
        "n": int(differences.size),
        "median_difference": float(np.median(differences)),
        "p_value": float(result.pvalue),
    }


def bootstrap_median_ci(
    values: np.ndarray, n_resamples: int = 2000, seed: int = 0
) -> tuple[float, float]:
    """Intervalle de confiance bootstrap (percentile, 95 %) de la médiane."""
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return math.nan, math.nan
    rng = np.random.default_rng(seed)
    samples = rng.choice(data, size=(n_resamples, data.size), replace=True)
    medians = np.median(samples, axis=1)
    return float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def spearman(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    """Corrélation de rang de Spearman."""
    if len(x) < 3:
        return {"rho": math.nan, "p_value": math.nan, "n": len(x)}
    result = stats.spearmanr(x, y)
    return {"rho": float(result.statistic), "p_value": float(result.pvalue), "n": len(x)}
