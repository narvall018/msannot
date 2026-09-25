"""Rapport de validation des similarités contre les implémentations de référence.

- cosine et modified cosine comparés à **matchms** (CosineGreedy, ModifiedCosine[Greedy]) ;
- entropie spectrale comparée à **ms_entropy** (Li et al.) ;
- moteur vectorisé comparé aux fonctions de référence paire par paire.

Les paires testées sont celles que le moteur trouve avec un score non nul, pour des requêtes
tirées au hasard (graine fixe) dans la bibliothèque de démonstration. Résultat écrit dans
``results/validation/similarity_validation.json`` et cité dans ``docs/methodology.md``.

Usage : ``python scripts/validate_similarity.py`` (extra ``ref`` requis).
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from msannot.config import SimilarityConfig
from msannot.io.msp import load_msp
from msannot.similarity.index import LibraryIndex
from msannot.similarity.pairwise import similarity

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "data" / "demo" / "massbank_open.msp.gz"
N_LIBRARY, N_QUERIES, PER_QUERY, SEED = 6000, 60, 200, 1


def main() -> int:
    try:
        import matchms
        import ms_entropy
        from matchms import similarity as reference
    except ImportError:
        print("matchms et ms_entropy sont nécessaires : pip install -e '.[ref]'")
        return 1
    config = SimilarityConfig()
    rng = np.random.default_rng(SEED)
    spectra = load_msp(LIBRARY)
    library = [spectra[i] for i in sorted(rng.choice(len(spectra), N_LIBRARY, replace=False))]
    queries = [library[i] for i in rng.choice(len(library), N_QUERIES, replace=False)]
    index = LibraryIndex(library, config)
    modified_class = getattr(reference, "ModifiedCosineGreedy", None) or reference.ModifiedCosine
    cosine_ref = reference.CosineGreedy(
        tolerance=config.tolerance, intensity_power=config.intensity_power
    )
    modified_ref = modified_class(
        tolerance=config.tolerance, intensity_power=config.intensity_power
    )

    def to_matchms(s):
        return matchms.Spectrum(
            mz=np.asarray(s.mz),
            intensities=np.asarray(s.intensities),
            metadata={"precursor_mz": s.precursor_mz},
        )

    report: dict[str, dict[str, float | int]] = {}
    # 1) moteur vectorisé contre fonctions de référence
    for metric in config.metrics:
        differences, disagreements, pairs = [], 0, 0
        for query in queries:
            dense, matches = index.query(query, metric).dense(len(library))
            for position in np.flatnonzero(dense > 0)[:PER_QUERY]:
                expected = similarity(metric, library[position], query, config)
                differences.append(abs(expected.score - dense[position]))
                disagreements += int(expected.n_matches != matches[position])
                pairs += 1
        report[f"index_vs_reference_{metric}"] = {
            "pairs": pairs,
            "max_abs_difference": float(max(differences)),
            "match_count_disagreements": disagreements,
        }
    # 2) matchms
    for metric, ref in (("cosine", cosine_ref), ("modified_cosine", modified_ref)):
        differences, disagreements, pairs = [], 0, 0
        for query in queries[:30]:
            dense, _ = index.query(query, "modified_cosine").dense(len(library))
            for position in np.flatnonzero(dense > 0)[:PER_QUERY]:
                ours = similarity(metric, library[position], query, config)
                theirs = ref.pair(to_matchms(library[position]), to_matchms(query))
                differences.append(abs(ours.score - float(theirs["score"])))
                disagreements += int(ours.n_matches != int(theirs["matches"]))
                pairs += 1
        report[f"matchms_{metric}"] = {
            "pairs": pairs,
            "max_abs_difference": float(max(differences)),
            "match_count_disagreements": disagreements,
        }

    # 3) ms_entropy, avec diagnostic des écarts
    def peaks(s):
        return np.column_stack([s.mz, s.intensities / s.intensities.sum()]).astype(np.float32)

    rows = []
    for query in queries[:30]:
        dense, _ = index.query(query, "entropy").dense(len(library))
        for position in np.flatnonzero(dense > 0)[:PER_QUERY]:
            a = library[position]
            theirs = ms_entropy.calculate_entropy_similarity(
                peaks(a), peaks(query), ms2_tolerance_in_da=config.tolerance, clean_spectra=False
            )
            close = any(
                s.n_peaks > 1 and np.min(np.diff(s.mz)) <= 2 * config.tolerance for s in (a, query)
            )
            boundary = bool(
                np.any(np.abs(np.abs(a.mz[:, None] - query.mz[None, :]) - config.tolerance) < 1e-4)
            )
            rows.append((abs(theirs - dense[position]), close, boundary))
    differences = np.array([r[0] for r in rows])
    close = np.array([r[1] for r in rows])
    boundary = np.array([r[2] for r in rows])
    large = differences > 1e-3
    report["ms_entropy"] = {
        "pairs": len(rows),
        "median_abs_difference": float(np.median(differences)),
        "pairs_above_1e-3": int(large.sum()),
        "above_1e-3_with_close_peaks": int((large & close).sum()),
        "above_1e-3_at_tolerance_boundary_only": int((large & ~close & boundary).sum()),
        "above_1e-3_unexplained": int((large & ~close & ~boundary).sum()),
        "max_abs_difference_unambiguous": float(differences[~close & ~boundary].max()),
    }
    output = ROOT / "results" / "validation" / "similarity_validation.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    report["versions"] = {
        "matchms": matchms.__version__,
        "ms_entropy": getattr(ms_entropy, "__version__", "?"),
    }
    output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
