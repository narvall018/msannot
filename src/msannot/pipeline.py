"""Orchestration : préparation d'une bibliothèque et benchmark complet.

- :func:`prepare_dataset` : curation d'un export MassBank → bibliothèque MSP + bilan +
  fichier d'attribution (licences) ;
- :func:`load_library` / :func:`run_benchmark` : évaluation complète, écriture des tableaux,
  figures, résumé JSON et rapport HTML.

Le calcul (modules ``evaluation``) est séparé de l'écriture (ce module) ; le dashboard
réutilise directement :class:`~msannot.search.SpectralLibrary`.
"""

from __future__ import annotations

import platform
import sys
import time
from dataclasses import asdict, dataclass
from datetime import UTC, datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, cast

import numpy as np
import pandas as pd

from msannot import __version__
from msannot.chem.structures import FingerprintIndex, unique_structures
from msannot.config import BenchmarkConfig, Metric, PrepareConfig
from msannot.evaluation import (
    AnalogResult,
    IdentificationResult,
    build_network,
    evaluate_analogs,
    evaluate_identification,
    select_analog_queries,
)
from msannot.exceptions import AnalysisError
from msannot.io.common import read_json, sha256sum, write_json, write_tsv
from msannot.io.massbank import read_license_map
from msannot.io.msp import iter_msp_records, load_msp, write_msp
from msannot.logging_utils import get_logger
from msannot.plotting import (
    analog_figure,
    calibration_figure,
    curation_figure,
    dataset_figure,
    identification_figure,
    mirror_figure,
    network_figure,
    relation_figure,
    save_figure,
)
from msannot.processing.cleaning import clean_all
from msannot.processing.curation import CurationReport, curate
from msannot.search.library import SpectralLibrary, rank_order

logger = get_logger(__name__)

TRACKED_PACKAGES = ("rdkit", "numpy", "scipy", "pandas", "matplotlib", "networkx", "pydantic")


def environment_versions() -> dict[str, str]:
    """Versions de Python, du système et des dépendances principales."""
    versions = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "msannot": __version__,
    }
    for package in TRACKED_PACKAGES:
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = "absent"
    return versions


def sidecar(path: Path, suffix: str) -> Path:
    """Fichier associé à une bibliothèque : ``lib.msp.gz`` → ``lib.<suffix>``."""
    name = path.name
    for ending in (".msp.gz", ".msp", ".mgf"):
        if name.endswith(ending):
            name = name[: -len(ending)]
            break
    return path.with_name(f"{name}.{suffix}")


# ---------------------------------------------------------------------------- préparation
@dataclass(frozen=True, slots=True)
class PrepareOutcome:
    """Fichiers produits par la curation."""

    library: Path
    report_path: Path
    attribution_path: Path
    report: CurationReport


def prepare_dataset(config: PrepareConfig, exclude: frozenset[str] = frozenset()) -> PrepareOutcome:
    """Curation d'un export MSP (avec licences si fournies) → bibliothèque MSP compressée.

    Parameters
    ----------
    config
        Critères de curation.
    exclude
        Accessions à retirer de la bibliothèque (par exemple des requêtes de démonstration).
    """
    licenses = read_license_map(config.licenses_json) if config.licenses_json else None
    spectra, report = curate(iter_msp_records(config.input_msp), config, licenses)
    spectra = [spectrum for spectrum in spectra if spectrum.identifier not in exclude]
    write_msp(spectra, config.output)
    attribution = pd.DataFrame(
        {
            "identifier": [s.identifier for s in spectra],
            "name": [s.name for s in spectra],
            "contributor": [s.contributor for s in spectra],
            "license": [s.get("license") for s in spectra],
            "source": [s.get("source") for s in spectra],
        }
    )
    attribution_path = sidecar(config.output, "attribution.tsv.gz")
    attribution.to_csv(
        attribution_path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0}
    )
    report_path = sidecar(config.output, "curation.json")
    write_json(
        {
            **asdict(report),
            "n_written": len(spectra),
            "excluded_identifiers": sorted(exclude),
            "input": {"msp": config.input_msp.name, "sha256": sha256sum(config.input_msp)},
            "created_utc": datetime.now(UTC),
            "versions": environment_versions(),
        },
        report_path,
    )
    return PrepareOutcome(config.output, report_path, attribution_path, report)


# ---------------------------------------------------------------------------- benchmark
def load_library(config: BenchmarkConfig) -> tuple[SpectralLibrary, dict[str, int]]:
    """Charge, nettoie et indexe la bibliothèque (spectres avec structure uniquement)."""
    spectra = load_msp(config.dataset)
    with_structure = [s for s in spectra if s.inchikey14 and s.smiles]
    cleaned, rejected = clean_all(with_structure, config.cleaning)
    if len(cleaned) < 20:
        raise AnalysisError(
            f"trop peu de spectres utilisables ({len(cleaned)}) dans {config.dataset}"
        )
    info = {
        "n_loaded": len(spectra),
        "n_without_structure": len(spectra) - len(with_structure),
        "n_rejected_cleaning": rejected,
        "n_library": len(cleaned),
    }
    logger.info("Bibliothèque : %d spectres utilisables sur %d", len(cleaned), len(spectra))
    return SpectralLibrary(cleaned, config.similarity), info


@dataclass(frozen=True, slots=True)
class Example:
    """Exemple commenté (spectre miroir) : requête, résultat, mesure, contexte."""

    kind: str
    title: str
    query_index: int
    hit_index: int
    metric: str
    tanimoto: float


def _preferred(metrics: tuple[Metric, ...], choice: Metric) -> Metric:
    return choice if choice in metrics else metrics[0]


def _best_candidate(
    library: SpectralLibrary, query_index: int, metric: Metric, ppm: float, correct: bool
) -> int | None:
    """Meilleur candidat correct (ou incorrect) d'une requête, scénario autres laboratoires."""
    query = library.spectra[query_index]
    allowed = (
        library.precursor_window(query.precursor_mz, ppm)
        & (library.contributors != query.contributor)
        & (library.duplicates != library.duplicates[query_index])
        & ((library.inchikeys == query.inchikey14) == correct)
    )
    if not allowed.any():
        return None
    scores, matches = library.scores(query, metric, allowed)
    best = int(rank_order(np.where(allowed, scores, -1.0), matches)[0])
    return best if allowed[best] else None


def select_examples(
    library: SpectralLibrary,
    identification: IdentificationResult,
    analogs: AnalogResult,
    fingerprints: FingerprintIndex,
    config: BenchmarkConfig,
) -> list[Example]:
    """Choisit des exemples par des règles fixes (aucune sélection manuelle).

    - réussite « typique » : score du bon composé le plus proche de la médiane des réussites ;
    - erreur la plus nette : plus grand écart entre le meilleur mauvais composé et le bon ;
    - analogue : meilleur score parmi les analogues proches (0,5 ≤ Tanimoto < 0,95) portant
      une vraie modification (écart de masse ≥ 1 Da).
    Les exemples d'identification relèvent du scénario « autres laboratoires ».
    """
    examples: list[Example] = []
    metric = _preferred(config.similarity.metrics, "entropy")
    ppm = config.identification.precursor_ppm
    rows = identification.queries
    realistic = rows[
        (rows["setting"] == "other_labs")
        & (rows["metric"] == metric)
        & rows["reachable"]
        & (rows["n_competitor_compounds"] > 0)
    ]
    correct = realistic["top1_correct"].astype(bool)
    successes = realistic[correct]
    if not successes.empty:
        distance = (successes["correct_score"] - successes["correct_score"].median()).abs()
        successes = successes.assign(order=distance).sort_values(["order", "query_index"])
    failures = realistic[~correct]
    if not failures.empty:
        margin = failures["top1_score"] - failures["correct_score"]
        failures = failures.assign(order=-margin).sort_values(["order", "query_index"])
    for kind, title, subset, want_correct in (
        (
            "success",
            "Identification typique réussie malgré des composés de même masse",
            successes,
            True,
        ),
        (
            "failure",
            "Erreur la plus nette : un autre composé devance nettement le bon",
            failures,
            False,
        ),
    ):
        if subset.empty:
            continue
        query_index = int(subset.iloc[0]["query_index"])
        hit_index = _best_candidate(library, query_index, metric, ppm, want_correct)
        if hit_index is None:
            continue
        query = library.spectra[query_index]
        tanimoto = fingerprints.similarity(query.inchikey14, library.spectra[hit_index].inchikey14)
        examples.append(Example(kind, title, query_index, hit_index, metric, tanimoto))

    analog_metric = _preferred(config.similarity.metrics, "modified_cosine")
    candidates = analogs.queries[
        (analogs.queries["metric"] == analog_metric)
        & (analogs.queries["top1_tanimoto"] >= 0.5)
        & (analogs.queries["top1_tanimoto"] < 0.95)
        & (analogs.queries["top1_delta_mz"].abs() >= 1.0)  # vraie modification chimique
    ].sort_values(["top1_score", "query_index"], ascending=[False, True])
    if not candidates.empty:
        row = candidates.iloc[0]
        identifiers = library.table["identifier"].to_numpy()
        hit_index = int(np.flatnonzero(identifiers == row["top1_id"])[0])
        examples.append(
            Example(
                "analog",
                "Analogue structural trouvé (molécule recherchée absente de la bibliothèque)",
                int(row["query_index"]),
                hit_index,
                analog_metric,
                float(row["top1_tanimoto"]),
            )
        )
    return examples


@dataclass(frozen=True, slots=True)
class BenchmarkOutcome:
    """Bilan d'un benchmark."""

    output_dir: Path
    summary: dict[str, Any]
    report_path: Path | None


def _calibration_points(calibration: pd.DataFrame, target: float = 0.95) -> list[dict[str, Any]]:
    """Plus petit seuil donnant une précision ≥ ``target`` (scénario autres laboratoires)."""
    points = []
    data = calibration[calibration["setting"] == "other_labs"]
    for metric, group in data.groupby("metric", sort=False):
        ok = group[(group["precision"] >= target) & (group["n_retained"] >= 30)]
        if ok.empty:
            points.append(
                {"metric": metric, "threshold": None, "coverage": None, "precision": None}
            )
            continue
        best = ok.sort_values("threshold").iloc[0]
        points.append(
            {
                "metric": metric,
                "threshold": float(best["threshold"]),
                "coverage": float(best["coverage"]),
                "precision": float(best["precision"]),
                "n_retained": int(best["n_retained"]),
            }
        )
    return points


def run_benchmark(
    config: BenchmarkConfig,
    *,
    output_dir: Path | None = None,
    config_path: Path | None = None,
    make_report: bool = True,
) -> BenchmarkOutcome:
    """Exécute le benchmark complet et écrit tous les résultats."""
    output = (output_dir or config.output_dir).resolve()
    figures_dir = output / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    started = datetime.now(UTC)
    timings: dict[str, float] = {}
    clock = time.perf_counter()

    library, load_info = load_library(config)
    fingerprints = FingerprintIndex(
        unique_structures(list(zip(library.inchikeys, library.table["smiles"], strict=True))),
        radius=config.fingerprint.radius,
        n_bits=config.fingerprint.n_bits,
    )
    timings["chargement"] = time.perf_counter() - clock

    clock = time.perf_counter()
    identification = evaluate_identification(library, config)
    timings["identification"] = time.perf_counter() - clock
    clock = time.perf_counter()
    analogs = evaluate_analogs(library, fingerprints, config)
    timings["analogues"] = time.perf_counter() - clock
    network = None
    if config.network.enabled:
        clock = time.perf_counter()
        network = build_network(
            library,
            select_analog_queries(library),  # un nœud par composé
            fingerprints,
            config.similarity,
            config.network,
            config.analogs.close_analog_threshold,
            seed=config.analogs.seed,
        )
        timings["réseau"] = time.perf_counter() - clock

    clock = time.perf_counter()
    metrics = config.similarity.metrics
    curation = _curation(config)
    figures = {
        "dataset": dataset_figure(library.table),
        "identification": identification_figure(identification.summary, metrics),
        "calibration": calibration_figure(identification.calibration, metrics, "other_labs"),
        "analogs": analog_figure(analogs.queries, metrics, config.analogs.close_analog_threshold),
        "relation": relation_figure(analogs.relation_bins, metrics, analogs.baseline_tanimoto),
    }
    if curation is not None:
        figures["curation"] = curation_figure(
            curation["excluded"], curation["n_input"], curation["n_kept"]
        )
    if network is not None:
        figures["network"] = network_figure(
            network.graph, network.summary["random_tanimoto"], seed=config.analogs.seed
        )
    examples = select_examples(library, identification, analogs, fingerprints, config)
    for number, example in enumerate(examples, start=1):
        figures[f"example_{number}"] = mirror_figure(
            library.spectra[example.query_index],
            library.spectra[example.hit_index],
            cast(Metric, example.metric),
            config.similarity,
            title=example.title,
            tanimoto=example.tanimoto,
        )
    figure_files = {}
    for key, figure in figures.items():
        save_figure(figure, figures_dir / f"{key}.png", dpi=config.figure_dpi)
        figure_files[key] = f"figures/{key}.png"
    timings["figures"] = time.perf_counter() - clock

    write_tsv(identification.queries, output / "identification_queries.tsv")
    write_tsv(identification.summary, output / "identification_summary.tsv")
    write_tsv(identification.calibration, output / "identification_calibration.tsv")
    write_tsv(identification.mcnemar, output / "identification_mcnemar.tsv")
    write_tsv(analogs.queries, output / "analog_queries.tsv")
    write_tsv(analogs.summary, output / "analog_summary.tsv")
    write_tsv(analogs.tests, output / "analog_tests.tsv")
    write_tsv(analogs.relation_bins, output / "relation_bins.tsv")
    if network is not None:
        write_tsv(network.edges, output / "network_edges.tsv")

    summary: dict[str, Any] = {
        "created_utc": started,
        "finished_utc": datetime.now(UTC),
        "command": " ".join(
            Path(sys.argv[0]).name if i == 0 else arg for i, arg in enumerate(sys.argv)
        ),
        "config_file": config_path.name if config_path else None,
        "config": config.model_dump(mode="json"),
        "versions": environment_versions(),
        "dataset": {
            "file": config.dataset.name,
            "sha256": sha256sum(config.dataset),
            **load_info,
            "n_compounds": int(library.table["inchikey14"].nunique()),
            "n_contributors": int(library.table["contributor"].nunique()),
            "instrument_types": library.table["instrument_type"].value_counts().to_dict(),
            "licenses": library.table["license"].value_counts().to_dict(),
            "duplicates": library.duplicate_summary(),
        },
        "curation": curation,
        "identification": {
            "n_eligible": identification.n_eligible,
            "n_evaluated": identification.n_evaluated,
            "summary": identification.summary,
            "mcnemar": identification.mcnemar,
            "calibration_95": _calibration_points(identification.calibration),
        },
        "analogs": {
            "summary": analogs.summary,
            "tests": analogs.tests,
            "relation_stats": analogs.relation_stats,
            "relation_bins": analogs.relation_bins,
            "baseline_tanimoto": analogs.baseline_tanimoto,
        },
        "network": (
            {key: value for key, value in network.summary.items() if key != "random_tanimoto"}
            if network is not None
            else None
        ),
        "examples": [
            {
                **asdict(example),
                "query_id": library.spectra[example.query_index].identifier,
                "query_name": library.spectra[example.query_index].name,
                "hit_id": library.spectra[example.hit_index].identifier,
                "hit_name": library.spectra[example.hit_index].name,
            }
            for example in examples
        ],
        "figures": figure_files,
        "timings_seconds": timings,
    }
    write_json(summary, output / "benchmark_summary.json")
    report_path = None
    if make_report:
        from msannot.report.builder import build_report

        report_path = build_report(output)
    return BenchmarkOutcome(output, summary, report_path)


def _curation(config: BenchmarkConfig) -> dict[str, Any] | None:
    path = config.curation_report or sidecar(config.dataset, "curation.json")
    if not path.is_file():
        return None
    data = read_json(path)
    return data if isinstance(data, dict) else None
