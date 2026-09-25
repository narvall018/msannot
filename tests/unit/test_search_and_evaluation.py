"""Moteur de recherche, statistiques et protocoles d'évaluation sur une bibliothèque jouet."""

from pathlib import Path

import numpy as np
import pytest

from msannot.chem.structures import FingerprintIndex, unique_structures
from msannot.config import BenchmarkConfig, NetworkConfig
from msannot.evaluation import build_network, evaluate_analogs, evaluate_identification
from msannot.evaluation.stats import (
    Proportion,
    bootstrap_median_ci,
    mcnemar_exact,
    spearman,
    wilcoxon_paired,
    wilson_interval,
)
from msannot.search.library import duplicate_groups, rank_order

CONFIG = BenchmarkConfig(dataset=Path("jouet.msp"))


# --------------------------------------------------------------------------- recherche
def test_identity_search_finds_same_compound(toy_library):
    query = toy_library.spectra[0]
    exclude = np.zeros(len(toy_library), dtype=bool)
    exclude[0] = True
    hits = toy_library.search(query, metric="entropy", mode="identity", exclude=exclude)
    assert list(hits["rank"]) == list(range(1, len(hits) + 1))
    assert set(hits["inchikey14"]) == {"RYYVLZVUVIJVG"}  # seuls les spectres de même masse
    assert hits["score"].is_monotonic_decreasing


def test_analog_search_uses_mass_window(toy_library):
    query = toy_library.spectra[0]
    exclude = toy_library.inchikeys == query.inchikey14
    hits = toy_library.search(
        query, metric="modified_cosine", mode="analog", exclude=exclude, max_mass_shift=20
    )
    assert len(hits) > 0 and not (hits["inchikey14"] == query.inchikey14).any()
    assert (hits["delta_mz"].abs() <= 20).all()
    assert (
        len(
            toy_library.search(
                query, metric="cosine", mode="analog", exclude=exclude, max_mass_shift=1
            )
        )
        == 0
    )


def test_rank_order_breaks_ties_deterministically():
    order = rank_order(np.array([0.5, 0.9, 0.9, 0.1]), np.array([1, 2, 5, 1]))
    assert list(order) == [2, 1, 0, 3]


def test_duplicate_groups(make_spectrum):
    a = make_spectrum([100, 200], [1, 0.5], 300, "a")
    b = make_spectrum([100, 200], [1, 0.5], 300, "b")
    c = make_spectrum([100, 201], [1, 0.5], 300, "c")
    labels = duplicate_groups([a, b, c])
    assert labels[0] == labels[1] != labels[2]


# --------------------------------------------------------------------------- statistiques
def test_wilson_reference_values():
    low, high = wilson_interval(8, 10)
    assert (round(low, 4), round(high, 4)) == (0.4902, 0.9433)
    assert wilson_interval(0, 5)[0] == 0.0 and wilson_interval(5, 5)[1] == 1.0
    assert np.isnan(Proportion(0, 0).value)
    with pytest.raises(ValueError):
        wilson_interval(6, 5)


def test_mcnemar_exact():
    first = np.array([True] * 10 + [False] * 2 + [True] * 30)
    second = np.array([False] * 10 + [True] * 2 + [True] * 30)
    result = mcnemar_exact(first, second)
    assert (result["only_first"], result["only_second"]) == (10, 2)
    assert result["p_value"] == pytest.approx(0.0386, abs=1e-3)
    assert mcnemar_exact(first, first)["p_value"] == 1.0


def test_other_statistics():
    assert wilcoxon_paired(np.arange(10.0), np.arange(10.0))["p_value"] == 1.0
    shifted = wilcoxon_paired(np.arange(20.0) + 1, np.arange(20.0))
    assert shifted["median_difference"] == 1.0 and shifted["p_value"] < 0.001
    low, high = bootstrap_median_ci(np.arange(101.0), seed=1)
    assert low < 50 < high
    assert bootstrap_median_ci(np.arange(101.0), seed=1) == (low, high)
    assert spearman(np.arange(10), np.arange(10))["rho"] == pytest.approx(1.0)


# --------------------------------------------------------------------------- évaluation
def test_identification_protocol(toy_library):
    result = evaluate_identification(toy_library, CONFIG)
    rows = result.queries
    # la caféine de lab2 n'a de jumeau que dans lab1 : atteignable en « autres labos »
    caf2 = rows[(rows["query_id"] == "caf_lab2") & (rows["setting"] == "other_labs")]
    assert caf2["reachable"].all() and (caf2["n_candidates"] == 2).all()
    # théobromine : composé unique, non éligible
    assert "theo" not in set(rows["query_id"])
    # paraxanthine : un concurrent de même masse (théobromine)
    para = rows[(rows["query_id"] == "para") & (rows["setting"] == "all_labs")]
    assert (para["n_competitor_compounds"] == 1).all()
    assert (para["rank"] == 1).all()
    summary = result.summary
    assert set(summary["subset"]) == {"toutes", "avec concurrents"}
    assert summary["top1"].between(0, 1).all()


def test_analog_protocol_and_network(toy_library):
    fingerprints = FingerprintIndex(
        unique_structures(
            list(zip(toy_library.inchikeys, toy_library.table["smiles"], strict=True))
        )
    )
    result = evaluate_analogs(toy_library, fingerprints, CONFIG)
    queries = result.queries
    assert set(queries["metric"]) == set(CONFIG.similarity.metrics)
    assert (queries["top1_tanimoto"].dropna() <= queries["oracle_tanimoto"].max()).all()
    assert set(result.summary["method"]) >= {"random", "oracle"}
    network = build_network(
        toy_library,
        np.arange(len(toy_library)),
        fingerprints,
        CONFIG.similarity,
        NetworkConfig(min_score=0.5, min_matches=3),
        threshold=0.7,
    )
    assert network.summary["n_nodes"] == len(toy_library)
    assert network.summary["n_edges"] >= 1
    assert (network.edges["score"] >= 0.5).all()


def test_duplicate_summary(make_spectrum):
    from msannot.config import SimilarityConfig
    from msannot.search.library import SpectralLibrary

    peaks = ([100.0, 150.0, 200.0, 250.0, 280.0], [1.0, 0.5, 0.3, 0.2, 0.1])
    spectra = [
        make_spectrum(*peaks, 300.0, "a", inchikey14="AAAAAAAAAAAAAA", contributor="lab1"),
        make_spectrum(*peaks, 300.0, "b", inchikey14="AAAAAAAAAAAAAA", contributor="lab1"),
        # même spectre, autre laboratoire et autre annotation : conflit
        make_spectrum(*peaks, 300.0, "c", inchikey14="BBBBBBBBBBBBBB", contributor="lab2"),
        make_spectrum(
            [100.0, 160.0, 210.0, 260.0, 290.0],
            [1, 1, 1, 1, 1],
            300.0,
            "d",
            inchikey14="CCCCCCCCCCCCCC",
            contributor="lab2",
        ),
    ]
    summary = SpectralLibrary(spectra, SimilarityConfig()).duplicate_summary()
    assert summary == {
        "groups": 1,
        "spectra": 3,
        "cross_contributor_groups": 1,
        "conflicting_annotation_groups": 1,
    }
