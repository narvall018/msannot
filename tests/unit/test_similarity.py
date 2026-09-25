"""Mesures de similarité : propriétés mathématiques, cas limites, moteur vectorisé."""

import math

import numpy as np
import pytest
from tests.conftest import spectrum

from msannot.config import SimilarityConfig
from msannot.similarity import LibraryIndex, cosine, entropy_similarity, modified_cosine
from msannot.similarity.pairwise import find_matches, matched_peaks, similarity
from msannot.similarity.weights import (
    LN4,
    cosine_weights,
    entropy_pair_contribution,
    entropy_weights,
    spectral_entropy,
)

CONFIG = SimilarityConfig()


def test_identical_spectra_score_one(make_spectrum):
    s = make_spectrum([100, 150, 200, 250.5], [1, 0.5, 0.2, 0.9])
    for metric in CONFIG.metrics:
        assert similarity(metric, s, s, CONFIG).score == pytest.approx(1.0)


def test_disjoint_spectra_score_zero(make_spectrum):
    a = make_spectrum([100, 150], [1, 1])
    b = make_spectrum([101, 151], [1, 1], precursor=500.0)
    for metric in ("cosine", "entropy"):
        assert similarity(metric, a, b, CONFIG).score == 0.0


def test_cosine_hand_computed(make_spectrum):
    a = make_spectrum([100.0, 200.0], [1.0, 1.0])
    b = make_spectrum([100.0, 300.0], [1.0, 1.0])
    # un pic commun sur deux, intensités égales : cosine = 0,5
    assert cosine(a, b, 0.01, 0.5).score == pytest.approx(0.5)


def test_tolerance_boundary_is_inclusive(make_spectrum):
    a = make_spectrum([100.0], [1.0])
    b = make_spectrum([100.01], [1.0])
    assert cosine(a, b, 0.01, 1.0).n_matches == 1
    assert cosine(a, b, 0.009, 1.0).n_matches == 0


def test_each_peak_is_matched_once(make_spectrum):
    a = make_spectrum([100.000], [1.0])
    b = make_spectrum([99.995, 100.004], [1.0, 0.2])
    result = cosine(a, b, 0.01, 1.0)
    assert result.n_matches == 1
    assert result.score == pytest.approx(1.0 / math.sqrt(1.04))


def test_find_matches_order_and_shift():
    idx_a, idx_b = find_matches(np.array([100.0, 114.0]), np.array([86.0, 100.0]), 0.01, shift=14.0)
    assert list(zip(idx_a, idx_b, strict=True)) == [(0, 0), (1, 1)]


def test_modified_cosine_uses_precursor_shift(make_spectrum):
    # le spectre b porte une modification de +14,0157 (méthylation) sur deux fragments
    a = make_spectrum([80.0, 120.0, 150.0], [0.3, 1.0, 0.6], precursor=300.0)
    b = make_spectrum([80.0, 134.0157, 164.0157], [0.3, 1.0, 0.6], precursor=314.0157)
    plain = cosine(a, b, 0.01, 0.5)
    modified = modified_cosine(a, b, 0.01, 0.5)
    assert plain.n_matches == 1
    assert modified.n_matches == 3
    assert modified.score == pytest.approx(1.0)
    shifted = matched_peaks("modified_cosine", a, b, CONFIG).shifted
    assert shifted.sum() == 2


def test_modified_cosine_equals_cosine_for_same_precursor(make_spectrum):
    a = make_spectrum([80.0, 120.0, 150.0], [0.3, 1.0, 0.6], precursor=300.0)
    b = make_spectrum([80.0, 121.0, 150.0], [0.5, 1.0, 0.2], precursor=300.004)
    assert modified_cosine(a, b, 0.01, 0.5) == cosine(a, b, 0.01, 0.5)


def test_entropy_decomposition_matches_definition(make_spectrum):
    """1 − (2·S_AB − S_A − S_B)/ln 4 (définition) == somme des contributions / ln 4."""
    a_int = np.array([0.5, 0.3, 0.2])
    b_int = np.array([0.6, 0.4])
    a = make_spectrum([100, 150, 200], a_int)
    b = make_spectrum([100, 175], b_int)
    pa, pb = a_int / a_int.sum(), b_int / b_int.sum()
    merged = np.array([pa[0] + pb[0], pa[1], pa[2], pb[1]]) / 2
    definition = (
        1 - (2 * spectral_entropy(merged) - spectral_entropy(pa) - spectral_entropy(pb)) / LN4
    )
    result = entropy_similarity(a, b, 0.01, weighted=False)
    assert result.score == pytest.approx(definition)
    assert entropy_pair_contribution(np.array([0.5]), np.array([0.5]))[0] == pytest.approx(
        math.log(2)
    )


def test_entropy_weighting_only_below_three_nats():
    poor = np.array([1.0, 0.1, 0.05])
    rich = np.ones(40)
    assert spectral_entropy(poor) < 3 < spectral_entropy(rich)
    weighted = entropy_weights(poor, weighted=True)
    assert weighted.sum() == pytest.approx(1.0)
    assert weighted[1] > poor[1] / poor.sum()  # les petits pics sont renforcés
    assert np.allclose(entropy_weights(rich, weighted=True), rich / rich.sum())


def test_cosine_weights_unit_norm():
    weights = cosine_weights(np.array([4.0, 1.0]), 0.5)
    assert np.linalg.norm(weights) == pytest.approx(1.0)
    assert weights[0] / weights[1] == pytest.approx(2.0)


# --------------------------------------------------------------------------- index vectorisé
@pytest.mark.parametrize("metric", ["cosine", "modified_cosine", "entropy"])
def test_index_matches_reference_implementation(random_spectra, metric):
    index = LibraryIndex(random_spectra, CONFIG)
    for query in random_spectra[::7]:
        dense, matches = index.query(query, metric).dense(len(random_spectra))
        for position, reference in enumerate(random_spectra):
            expected = similarity(metric, reference, query, CONFIG)
            assert dense[position] == pytest.approx(expected.score, abs=1e-12)
            assert matches[position] == expected.n_matches


def test_index_respects_candidate_mask(random_spectra):
    index = LibraryIndex(random_spectra, CONFIG)
    mask = np.zeros(len(random_spectra), dtype=bool)
    mask[::3] = True
    result = index.query(random_spectra[0], "cosine", mask)
    assert set(result.indices) <= set(np.flatnonzero(mask))


def test_index_handles_conflicting_pairs():
    """Deux pics requête dans la tolérance d'un même pic : glouton exact requis."""
    library = [spectrum([100.0, 150.0], [1.0, 0.5], 300.0, "lib")]
    query = spectrum([99.996, 100.004, 150.0], [0.2, 1.0, 0.5], 300.0, "q")
    index = LibraryIndex(library, CONFIG)
    result = index.query(query, "cosine")
    reference = similarity("cosine", library[0], query, CONFIG)
    assert result.n_matches[0] == reference.n_matches == 2
    assert result.scores[0] == pytest.approx(reference.score)


def test_index_empty_result():
    index = LibraryIndex([spectrum([100.0, 150.0], [1, 1], 300.0)], CONFIG)
    result = index.query(spectrum([500.0, 600.0], [1, 1], 700.0), "entropy")
    assert result.indices.size == 0


def test_index_rejects_empty_library():
    with pytest.raises(ValueError):
        LibraryIndex([], CONFIG)
