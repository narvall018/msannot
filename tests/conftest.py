"""Fixtures partagées : spectres synthétiques, petite bibliothèque, bibliothèque de démo."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from msannot.config import SimilarityConfig
from msannot.io.msp import load_msp
from msannot.models import Spectrum
from msannot.search.library import SpectralLibrary

ROOT = Path(__file__).resolve().parents[1]
DEMO_LIBRARY = ROOT / "data" / "demo" / "massbank_open.msp.gz"
EXAMPLE_QUERIES = ROOT / "data" / "demo" / "example_queries.mgf"

# Quelques molécules réelles (SMILES) pour les tests qui ont besoin de structures.
MOLECULES = {
    "caffeine": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "theobromine": "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
    "paraxanthine": "Cn1cnc2c1c(=O)n(C)c(=O)[nH]2",
    "carbamazepine": "NC(=O)N1c2ccccc2C=Cc2ccccc21",
    "phenol": "Oc1ccccc1",
}


def spectrum(
    mz: list[float] | np.ndarray,
    intensities: list[float] | np.ndarray,
    precursor: float = 500.0,
    identifier: str = "s",
    **metadata: str,
) -> Spectrum:
    """Construit un spectre de test."""
    return Spectrum(
        identifier, np.asarray(mz, float), np.asarray(intensities, float), precursor, metadata
    )


@pytest.fixture
def make_spectrum():
    return spectrum


def random_spectrum(rng: np.random.Generator, identifier: str, precursor: float) -> Spectrum:
    """Spectre aléatoire réaliste (m/z à 4 décimales, pics bien séparés ou non)."""
    n = int(rng.integers(5, 30))
    mz = np.round(np.sort(rng.uniform(50, precursor - 2, n)), 4)
    intensities = rng.gamma(0.8, 1.0, n)
    return Spectrum(identifier, mz, intensities / intensities.max(), precursor)


@pytest.fixture
def random_spectra() -> list[Spectrum]:
    rng = np.random.default_rng(42)
    spectra = []
    for index in range(80):
        precursor = float(np.round(rng.uniform(150, 600), 4))
        base = random_spectrum(rng, f"r{index}", precursor)
        spectra.append(base)
        # variantes proches (bruit sur m/z et intensités) pour créer des scores élevés
        noisy_mz = np.round(base.mz + rng.normal(0, 0.003, base.n_peaks), 4)
        noisy = base.intensities * rng.uniform(0.7, 1.3, base.n_peaks)
        spectra.append(
            Spectrum(f"r{index}b", np.sort(noisy_mz), noisy, precursor + rng.choice([0.0, 14.0157]))
        )
    return spectra


@pytest.fixture
def toy_library() -> SpectralLibrary:
    """Petite bibliothèque à structures connues : deux laboratoires, un isomère concurrent."""
    caffeine = MOLECULES["caffeine"]
    theobromine = MOLECULES["theobromine"]
    paraxanthine = MOLECULES["paraxanthine"]
    base = [110.07, 138.066, 123.043, 83.06, 69.045]
    spectra = [
        spectrum(
            base,
            [30, 100, 20, 10, 5],
            195.0877,
            "caf_lab1",
            name="caffeine",
            smiles=caffeine,
            inchikey14="RYYVLZVUVIJVG",
            contributor="lab1",
        ),
        spectrum(
            base,
            [28, 100, 22, 9, 6],
            195.0877,
            "caf_lab1b",
            name="caffeine",
            smiles=caffeine,
            inchikey14="RYYVLZVUVIJVG",
            contributor="lab1",
        ),
        spectrum(
            base,
            [35, 100, 15, 12, 4],
            195.0877,
            "caf_lab2",
            name="caffeine",
            smiles=caffeine,
            inchikey14="RYYVLZVUVIJVG",
            contributor="lab2",
        ),
        spectrum(
            [124.05, 138.066, 96.05, 67.04, 55.03],
            [100, 40, 30, 20, 10],
            181.0720,
            "theo",
            name="theobromine",
            smiles=theobromine,
            inchikey14="YAPQBXQYLJRXS",
            contributor="lab2",
        ),
        spectrum(
            [124.05, 138.066, 96.05, 69.045, 55.03],
            [100, 60, 20, 10, 10],
            181.0720,
            "para",
            name="paraxanthine",
            smiles=paraxanthine,
            inchikey14="QUNWUDVFRNGTC",
            contributor="lab1",
        ),
        spectrum(
            [124.05, 138.066, 96.05, 69.045, 55.03],
            [100, 50, 25, 10, 12],
            181.0720,
            "para2",
            name="paraxanthine",
            smiles=paraxanthine,
            inchikey14="QUNWUDVFRNGTC",
            contributor="lab2",
        ),
    ]
    return SpectralLibrary(spectra, SimilarityConfig())


@pytest.fixture(scope="session")
def demo_spectra() -> list[Spectrum]:
    return load_msp(DEMO_LIBRARY)
