"""Parité avec les implémentations de référence : matchms (cosine) et ms_entropy (entropie).

Ces tests utilisent de vrais spectres MassBank de la bibliothèque de démonstration. Ils sont
ignorés si les bibliothèques de référence ne sont pas installées (extra ``ref``).
"""

import numpy as np
import pytest

from msannot.config import SimilarityConfig
from msannot.similarity.index import LibraryIndex
from msannot.similarity.pairwise import similarity

CONFIG = SimilarityConfig()


@pytest.fixture(scope="module")
def sample(demo_spectra):
    rng = np.random.default_rng(3)
    library = [demo_spectra[i] for i in sorted(rng.choice(len(demo_spectra), 1500, replace=False))]
    return library, LibraryIndex(library, CONFIG)


def _related_pairs(sample, metric, n_queries=8, per_query=40):
    library, index = sample
    for query in library[:: len(library) // n_queries][:n_queries]:
        dense, _ = index.query(query, metric).dense(len(library))
        for position in np.flatnonzero(dense > 0)[:per_query]:
            yield library[position], query


def test_cosine_and_modified_cosine_match_matchms(sample):
    matchms = pytest.importorskip("matchms")
    from matchms import similarity as reference_similarity

    # ModifiedCosine a été renommé ModifiedCosineGreedy dans matchms 0.33.
    CosineGreedy = reference_similarity.CosineGreedy
    ModifiedCosineGreedy = (
        getattr(reference_similarity, "ModifiedCosineGreedy", None)
        or reference_similarity.ModifiedCosine
    )

    def convert(s):
        return matchms.Spectrum(
            mz=np.asarray(s.mz),
            intensities=np.asarray(s.intensities),
            metadata={"precursor_mz": s.precursor_mz},
        )

    references = {
        "cosine": CosineGreedy(tolerance=0.01, intensity_power=0.5),
        "modified_cosine": ModifiedCosineGreedy(tolerance=0.01, intensity_power=0.5),
    }
    checked = 0
    for metric, reference in references.items():
        for library_spectrum, query in _related_pairs(sample, "modified_cosine"):
            ours = similarity(metric, library_spectrum, query, CONFIG)
            theirs = reference.pair(convert(library_spectrum), convert(query))
            assert ours.score == pytest.approx(float(theirs["score"]), abs=1e-9)
            assert ours.n_matches == int(theirs["matches"])
            checked += 1
    assert checked > 200


def test_entropy_matches_ms_entropy_on_unambiguous_spectra(sample):
    """Accord à 1e-5 (float32 de ms_entropy) quand l'appariement n'est pas ambigu.

    ms_entropy suppose des pics séparés d'au moins deux tolérances ; et ses m/z en float32
    peuvent faire basculer les écarts exactement égaux à la tolérance. Ces deux cas, qui
    expliquent tous les écarts observés (voir docs/methodology.md), sont exclus ici.
    """
    ms_entropy = pytest.importorskip("ms_entropy")

    def unambiguous(a, b):
        spaced = all(s.n_peaks < 2 or np.min(np.diff(s.mz)) > 0.02 for s in (a, b))
        differences = np.abs(a.mz[:, None] - b.mz[None, :])
        return spaced and not np.any(np.abs(differences - 0.01) < 1e-3)

    def peaks(s):
        return np.column_stack([s.mz, s.intensities / s.intensities.sum()]).astype(np.float32)

    checked = 0
    for library_spectrum, query in _related_pairs(sample, "entropy", n_queries=10):
        if not unambiguous(library_spectrum, query):
            continue
        theirs = ms_entropy.calculate_entropy_similarity(
            peaks(library_spectrum), peaks(query), ms2_tolerance_in_da=0.01, clean_spectra=False
        )
        ours = similarity("entropy", library_spectrum, query, CONFIG).score
        assert ours == pytest.approx(theirs, abs=1e-5)
        checked += 1
    assert checked > 50
