"""Recherche d'analogues et relation entre similarité spectrale et similarité structurale.

Protocole :

- une requête par composé (le spectre d'accession la plus petite, choix arbitraire mais
  indépendant du contenu) ;
- **tous** les spectres du composé sont retirés de la bibliothèque : on simule une molécule
  inconnue ;
- candidats : spectres dont le précurseur diffère d'au plus ``max_mass_shift`` Da ;
- mesure : Tanimoto (Morgan) entre la requête et le **premier résultat**, comparé :
  - au **hasard** (Tanimoto moyen de tous les candidats) ;
  - à l'**oracle** (meilleure structure disponible parmi les candidats) ;
- une requête sans aucun pic commun avec la bibliothèque compte comme un échec pour le
  taux d'« analogues proches » (choix conservateur).

Les paires (requête, candidat) de score non nul servent aussi à mesurer si un score
spectral élevé annonce une structure proche (analyse par classes de score, Spearman).
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd

from msannot.chem.structures import FingerprintIndex
from msannot.config import BenchmarkConfig, Metric
from msannot.evaluation.stats import Proportion, bootstrap_median_ci, spearman, wilcoxon_paired
from msannot.logging_utils import get_logger
from msannot.search.library import SpectralLibrary, rank_order

logger = get_logger(__name__)

SCORE_BINS = np.round(np.arange(0.0, 1.0001, 0.1), 1)


@dataclass(frozen=True, slots=True)
class AnalogResult:
    """Résultats de la recherche d'analogues et de l'analyse spectre/structure."""

    queries: pd.DataFrame
    summary: pd.DataFrame
    tests: pd.DataFrame
    relation_bins: pd.DataFrame
    relation_stats: pd.DataFrame
    baseline_tanimoto: float


def select_analog_queries(
    library: SpectralLibrary, limit: int | None = None, seed: int = 0
) -> np.ndarray:
    """Un spectre par composé (accession la plus petite), échantillonné si ``limit``."""
    table = library.table.reset_index().rename(columns={"index": "position"})
    first = table.sort_values("identifier").groupby("inchikey14", sort=True)["position"].first()
    queries = np.sort(first.to_numpy())
    if limit is not None and queries.size > limit:
        rng = np.random.default_rng(seed)
        queries = np.sort(rng.choice(queries, size=limit, replace=False))
    return queries


def evaluate_analogs(
    library: SpectralLibrary, fingerprints: FingerprintIndex, config: BenchmarkConfig
) -> AnalogResult:
    """Évalue la recherche d'analogues pour toutes les mesures configurées."""
    metrics = config.similarity.metrics
    threshold = config.analogs.close_analog_threshold
    compound_of = fingerprints.indices(library.inchikeys)
    queries = select_analog_queries(library, config.analogs.n_queries, config.analogs.seed)
    logger.info("Analogues : %d requêtes (une par composé)", queries.size)

    rows: list[dict[str, object]] = []
    pair_scores: dict[str, list[np.ndarray]] = {metric: [] for metric in metrics}
    pair_tanimoto: dict[str, list[np.ndarray]] = {metric: [] for metric in metrics}
    baselines: list[float] = []
    for query_index in queries:
        query = library.spectra[int(query_index)]
        allowed = (library.inchikeys != query.inchikey14) & library.mass_shift_window(
            query.precursor_mz, config.analogs.max_mass_shift
        )
        if not allowed.any():
            continue
        tanimoto = fingerprints.similarity_to_all(query.inchikey14)[compound_of]
        oracle = float(tanimoto[allowed].max())
        random_mean = float(tanimoto[allowed].mean())
        random_close = float(np.mean(tanimoto[allowed] >= threshold))
        baselines.append(random_mean)
        for metric in metrics:
            scores, matches = library.scores(query, metric, allowed)
            masked = np.where(allowed, scores, -1.0)
            best = int(rank_order(masked, matches)[0])
            has_hit = scores[best] > 0
            nonzero = np.flatnonzero(scores > 0)
            pair_scores[metric].append(scores[nonzero])
            pair_tanimoto[metric].append(tanimoto[nonzero])
            rows.append(
                {
                    "query_index": int(query_index),
                    "query_id": query.identifier,
                    "inchikey14": query.inchikey14,
                    "metric": metric,
                    "n_candidates": int(allowed.sum()),
                    "has_hit": bool(has_hit),
                    "top1_id": library.spectra[best].identifier if has_hit else "",
                    "top1_score": float(scores[best]) if has_hit else 0.0,
                    "top1_matches": int(matches[best]) if has_hit else 0,
                    "top1_delta_mz": float(library.precursors[best] - query.precursor_mz)
                    if has_hit
                    else np.nan,
                    "top1_tanimoto": float(tanimoto[best]) if has_hit else np.nan,
                    "oracle_tanimoto": oracle,
                    "random_tanimoto": random_mean,
                    "random_close_fraction": random_close,
                }
            )
    table = pd.DataFrame(rows)
    relation_bins, relation_stats = _relation(pair_scores, pair_tanimoto, threshold)
    return AnalogResult(
        queries=table,
        summary=_summary(table, threshold, config.analogs.seed),
        tests=_tests(table, metrics),
        relation_bins=relation_bins,
        relation_stats=relation_stats,
        baseline_tanimoto=float(np.mean(baselines)) if baselines else np.nan,
    )


def _summary(table: pd.DataFrame, threshold: float, seed: int) -> pd.DataFrame:
    rows = []
    for metric, group in table.groupby("metric", sort=False):
        close = Proportion(int((group["top1_tanimoto"] >= threshold).sum()), len(group)).to_dict()
        low, high = bootstrap_median_ci(group["top1_tanimoto"].to_numpy(), seed=seed)
        rows.append(
            {
                "method": metric,
                "n_queries": len(group),
                "with_hit": float(group["has_hit"].mean()),
                "median_tanimoto": float(np.nanmedian(group["top1_tanimoto"])),
                "median_ci_low": low,
                "median_ci_high": high,
                "close_fraction": close["value"],
                "close_ci_low": close["ci95_low"],
                "close_ci_high": close["ci95_high"],
            }
        )
    reference = table.drop_duplicates("query_index")
    for method, column, close_value in (
        ("random", "random_tanimoto", float(reference["random_close_fraction"].mean())),
        (
            "oracle",
            "oracle_tanimoto",
            float((reference["oracle_tanimoto"] >= threshold).mean()),
        ),
    ):
        low, high = bootstrap_median_ci(reference[column].to_numpy(), seed=seed)
        rows.append(
            {
                "method": method,
                "n_queries": len(reference),
                "with_hit": 1.0,
                "median_tanimoto": float(reference[column].median()),
                "median_ci_low": low,
                "median_ci_high": high,
                "close_fraction": close_value,
                "close_ci_low": np.nan,
                "close_ci_high": np.nan,
            }
        )
    return pd.DataFrame(rows)


def _tests(table: pd.DataFrame, metrics: tuple[Metric, ...]) -> pd.DataFrame:
    wide = table.pivot(index="query_index", columns="metric", values="top1_tanimoto")
    random = table.drop_duplicates("query_index").set_index("query_index")["random_tanimoto"]
    rows: list[dict[str, object]] = []
    for first, second in combinations(metrics, 2):
        test = wilcoxon_paired(wide[first].to_numpy(), wide[second].to_numpy())
        rows.append({"a": first, "b": second, **test})
    for metric in metrics:
        test = wilcoxon_paired(wide[metric].to_numpy(), random.loc[wide.index].to_numpy())
        rows.append({"a": metric, "b": "random", **test})
    return pd.DataFrame(rows)


def _relation(
    pair_scores: dict[str, list[np.ndarray]],
    pair_tanimoto: dict[str, list[np.ndarray]],
    threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    bins_rows = []
    stats_rows = []
    for metric, chunks in pair_scores.items():
        scores = np.concatenate(chunks) if chunks else np.zeros(0)
        tanimoto = np.concatenate(pair_tanimoto[metric]) if chunks else np.zeros(0)
        classes = np.clip(np.digitize(scores, SCORE_BINS[1:-1]), 0, len(SCORE_BINS) - 2)
        for index in range(len(SCORE_BINS) - 1):
            values = tanimoto[classes == index]
            n = values.size
            mean = float(values.mean()) if n else np.nan
            half = 1.96 * float(values.std(ddof=1)) / np.sqrt(n) if n > 1 else np.nan
            bins_rows.append(
                {
                    "metric": metric,
                    "score_low": float(SCORE_BINS[index]),
                    "score_high": float(SCORE_BINS[index + 1]),
                    "n_pairs": int(n),
                    "mean_tanimoto": mean,
                    "ci_low": mean - half if n > 1 else np.nan,
                    "ci_high": mean + half if n > 1 else np.nan,
                    "fraction_close": float(np.mean(values >= threshold)) if n else np.nan,
                }
            )
        stats_rows.append({"metric": metric, **spearman(scores, tanimoto)})
    return pd.DataFrame(bins_rows), pd.DataFrame(stats_rows)
