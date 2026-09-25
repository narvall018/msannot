"""Évaluation de l'identification : le composé recherché est présent dans la bibliothèque.

Protocole (voir ``docs/methodology.md``) :

- **requêtes** : chaque spectre dont le composé (premier bloc d'InChIKey) possède au moins
  un autre spectre ; le spectre requête et ses **doublons exacts** (même liste de pics,
  déposée plusieurs fois) sont retirés des candidats ;
- **candidats** : spectres dont le précurseur est à moins de ``precursor_ppm`` ;
- deux scénarios :
  - ``all_labs`` : tous les candidats ;
  - ``other_labs`` : seuls les candidats d'**autres laboratoires**. C'est le cas réaliste : un
    spectre expérimental n'a pas de jumeau acquis sur le même instrument ;
- **rang pessimiste** : en cas d'égalité de score, le bon composé est classé après les
  candidats incorrects ex æquo (évite de flatter la méthode) ;
- **référence aléatoire** : probabilité de tomber sur le bon composé en choisissant un
  candidat au hasard dans la fenêtre de précurseur (ce que donne la masse seule) ;
- **requêtes avec concurrents** : requêtes dont la fenêtre de précurseur contient au moins
  un autre composé (isomère ou isobare). Pour les autres, la masse seule suffit : seules
  les premières mesurent le pouvoir discriminant des spectres. La calibration des scores
  est calculée sur ce sous-ensemble.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd

from msannot.config import BenchmarkConfig, Metric
from msannot.evaluation.stats import Proportion, mcnemar_exact
from msannot.logging_utils import get_logger
from msannot.search.library import SpectralLibrary

logger = get_logger(__name__)

SETTINGS = ("all_labs", "other_labs")
SUBSETS = ("toutes", "avec concurrents")
SETTING_LABELS = {
    "all_labs": "Tous laboratoires",
    "other_labs": "Autres laboratoires uniquement",
}
THRESHOLDS = np.round(np.arange(0.0, 1.0001, 0.05), 2)


@dataclass(frozen=True, slots=True)
class IdentificationResult:
    """Résultats détaillés et agrégés de l'identification."""

    queries: pd.DataFrame
    summary: pd.DataFrame
    calibration: pd.DataFrame
    mcnemar: pd.DataFrame
    n_eligible: int
    n_evaluated: int


def select_identification_queries(library: SpectralLibrary, config: BenchmarkConfig) -> np.ndarray:
    """Spectres dont le composé a au moins un autre spectre (échantillonnés si demandé)."""
    counts = pd.Series(library.inchikeys).value_counts()
    eligible = np.flatnonzero(pd.Series(library.inchikeys).map(counts).to_numpy() >= 2)
    limit = config.identification.max_queries
    if limit is not None and eligible.size > limit:
        rng = np.random.default_rng(config.identification.seed)
        eligible = np.sort(rng.choice(eligible, size=limit, replace=False))
    return eligible


def _query_rows(
    library: SpectralLibrary, query_index: int, metrics: tuple[Metric, ...], ppm: float
) -> list[dict[str, object]]:
    query = library.spectra[query_index]
    window = library.precursor_window(query.precursor_mz, ppm)
    # La requête et ses doublons exacts ne sont pas des candidats indépendants.
    window &= library.duplicates != library.duplicates[query_index]
    candidates = np.flatnonzero(window)
    correct = library.inchikeys[candidates] == query.inchikey14
    other_lab = library.contributors[candidates] != query.contributor
    rows: list[dict[str, object]] = []
    base = {
        "query_index": query_index,
        "query_id": query.identifier,
        "inchikey14": query.inchikey14,
        "contributor": query.contributor,
    }
    for metric in metrics:
        scores = library.scores(query, metric, window)[0][candidates]
        for setting, allowed in (("all_labs", np.ones_like(correct)), ("other_labs", other_lab)):
            allowed_correct = correct & allowed
            row: dict[str, object] = {
                **base,
                "metric": metric,
                "setting": setting,
                "n_candidates": int(allowed.sum()),
                "n_correct": int(allowed_correct.sum()),
                "reachable": bool(allowed_correct.any()),
            }
            competitors = allowed & ~correct
            row["n_competitor_compounds"] = int(
                np.unique(library.inchikeys[candidates][competitors]).size
            )
            if allowed_correct.any():
                best_correct = float(scores[allowed_correct].max())
                rank = 1 + int(np.sum(allowed & ~correct & (scores >= best_correct)))
                row.update(
                    rank=rank,
                    top1_correct=rank == 1,
                    top1_score=float(scores[allowed].max()),
                    correct_score=best_correct,
                    random_top1=float(allowed_correct.sum() / allowed.sum()),
                )
            rows.append(row)
    return rows


def with_competitors(queries: pd.DataFrame) -> pd.DataFrame:
    """Requêtes dont la fenêtre de précurseur contient au moins un autre composé."""
    return queries[queries["n_competitor_compounds"] > 0]


def _summarise(queries: pd.DataFrame, top_k: tuple[int, ...]) -> pd.DataFrame:
    frames = []
    for subset, data in zip(SUBSETS, (queries, with_competitors(queries)), strict=True):
        frame = _summarise_subset(data, top_k)
        frame.insert(1, "subset", subset)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True)


def _summarise_subset(queries: pd.DataFrame, top_k: tuple[int, ...]) -> pd.DataFrame:
    rows = []
    for (setting, metric), group in queries.groupby(["setting", "metric"], sort=False):
        row: dict[str, object] = {"setting": setting, "metric": metric, "n_queries": len(group)}
        for k in top_k:
            proportion = Proportion(int((group["rank"] <= k).sum()), len(group))
            details = proportion.to_dict()
            row[f"top{k}"] = details["value"]
            row[f"top{k}_ci_low"] = details["ci95_low"]
            row[f"top{k}_ci_high"] = details["ci95_high"]
        row["random_top1"] = float(group["random_top1"].mean())
        row["median_candidates"] = float(group["n_candidates"].median())
        rows.append(row)
    return pd.DataFrame(rows)


def _calibration(queries: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (setting, metric), group in queries.groupby(["setting", "metric"], sort=False):
        for threshold in THRESHOLDS:
            retained = group[(group["top1_score"] >= threshold) & (group["top1_score"] > 0)]
            proportion = Proportion(int(retained["top1_correct"].sum()), len(retained))
            details = proportion.to_dict()
            rows.append(
                {
                    "setting": setting,
                    "metric": metric,
                    "threshold": float(threshold),
                    "n_retained": len(retained),
                    "coverage": len(retained) / len(group) if len(group) else np.nan,
                    "precision": details["value"],
                    "precision_ci_low": details["ci95_low"],
                    "precision_ci_high": details["ci95_high"],
                }
            )
    return pd.DataFrame(rows)


def _mcnemar(queries: pd.DataFrame, metrics: tuple[Metric, ...]) -> pd.DataFrame:
    rows = []
    for setting, group in queries.groupby("setting", sort=False):
        wide = group.pivot(index="query_index", columns="metric", values="top1_correct")
        for first, second in combinations(metrics, 2):
            test = mcnemar_exact(wide[first].to_numpy(bool), wide[second].to_numpy(bool))
            rows.append(
                {
                    "setting": setting,
                    "metric_a": first,
                    "metric_b": second,
                    "top1_a": float(wide[first].mean()),
                    "top1_b": float(wide[second].mean()),
                    **test,
                }
            )
    return pd.DataFrame(rows)


def evaluate_identification(
    library: SpectralLibrary, config: BenchmarkConfig
) -> IdentificationResult:
    """Évalue l'identification pour toutes les mesures de similarité configurées."""
    metrics = config.similarity.metrics
    queries = select_identification_queries(library, config)
    logger.info("Identification : %d requêtes, mesures %s", queries.size, ", ".join(metrics))
    rows: list[dict[str, object]] = []
    for position, query_index in enumerate(queries, start=1):
        rows.extend(
            _query_rows(library, int(query_index), metrics, config.identification.precursor_ppm)
        )
        if position % 5000 == 0:
            logger.info("  %d / %d requêtes", position, queries.size)
    table = pd.DataFrame(rows)
    evaluated = table[table["reachable"]].copy()
    evaluated["rank"] = evaluated["rank"].astype(int)
    evaluated["top1_correct"] = evaluated["top1_correct"].astype(bool)
    return IdentificationResult(
        queries=table,
        summary=_summarise(evaluated, config.identification.top_k),
        calibration=_calibration(with_competitors(evaluated)),
        mcnemar=_mcnemar(with_competitors(evaluated), metrics),
        n_eligible=int(queries.size),
        n_evaluated=int(evaluated["query_index"].nunique()),
    )
