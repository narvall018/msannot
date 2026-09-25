"""Chimie (RDKit), modèle de spectre, nettoyage et curation."""

import numpy as np
import pytest
from tests.conftest import MOLECULES

from msannot.chem.structures import (
    FingerprintIndex,
    molecule_png,
    molecule_svg,
    parse_structure,
    ppm_error,
    protonated_mz,
)
from msannot.config import CleaningConfig, PrepareConfig
from msannot.io.msp import RawRecord
from msannot.models import Spectrum
from msannot.processing.cleaning import clean_spectrum
from msannot.processing.curation import REASONS, curate


# --------------------------------------------------------------------------- chimie
def test_caffeine_properties():
    info = parse_structure(MOLECULES["caffeine"])
    assert info is not None
    assert info.inchikey == "RYYVLZVUVIJVGH-UHFFFAOYSA-N"
    assert info.formula == "C8H10N4O2"
    assert protonated_mz(info) == pytest.approx(195.08765, abs=1e-4)
    assert ppm_error(195.0877, protonated_mz(info)) == pytest.approx(0.0, abs=1.0)


def test_invalid_smiles():
    assert parse_structure("pas_un_smiles") is None
    assert parse_structure("") is None
    assert molecule_png("C1CC") is None


def test_fingerprint_similarity():
    index = FingerprintIndex(dict(MOLECULES))
    assert index.similarity("caffeine", "caffeine") == 1.0
    assert index.similarity("caffeine", "theobromine") == index.similarity(
        "theobromine", "caffeine"
    )
    assert index.similarity("caffeine", "theobromine") > index.similarity("caffeine", "phenol")
    values = index.similarity_to_all("caffeine")
    assert values.shape == (len(MOLECULES),)
    assert values[index.position["caffeine"]] == 1.0


def test_molecule_drawings():
    png = molecule_png(MOLECULES["phenol"])
    assert png is not None and png.startswith(b"\x89PNG")
    svg = molecule_svg(MOLECULES["phenol"])
    assert svg is not None and "<svg" in svg


# --------------------------------------------------------------------------- spectre
def test_spectrum_validation_and_sorting():
    s = Spectrum("x", np.array([200.0, 100.0]), np.array([0.5, 1.0]), 300.0)
    assert list(s.mz) == [100.0, 200.0] and list(s.intensities) == [1.0, 0.5]
    with pytest.raises(ValueError):
        s.mz[0] = 5.0  # tableaux en lecture seule
    for mz, intensities, precursor in (
        ([100.0], [1.0, 2.0], 300.0),
        ([100.0], [-1.0], 300.0),
        ([100.0], [np.nan], 300.0),
        ([100.0], [1.0], 0.0),
    ):
        with pytest.raises(ValueError):
            Spectrum("bad", np.array(mz), np.array(intensities), precursor)


# --------------------------------------------------------------------------- nettoyage
def test_cleaning_steps():
    config = CleaningConfig(max_peaks=3, min_peaks=2)
    s = Spectrum(
        "x",
        np.array([50.0, 60.0, 70.0, 80.0, 90.0, 299.0]),
        np.array([100.0, 0.5, 40.0, 30.0, 20.0, 500.0]),
        300.0,
    )
    cleaned = clean_spectrum(s, config)
    assert cleaned is not None
    assert list(cleaned.mz) == [50.0, 70.0, 80.0]  # précurseur, bruit et 4e pic retirés
    assert cleaned.intensities.max() == 1.0
    assert clean_spectrum(s, CleaningConfig(min_peaks=10)) is None


# --------------------------------------------------------------------------- curation
def _record(identifier, smiles, precursor, **extra):
    metadata = {
        "identifier": identifier,
        "smiles": smiles,
        "precursor_mz": str(precursor),
        "precursor_type": "[M+H]+",
        "ion_mode": "POSITIVE",
        "spectrum_type": "MS2",
        "instrument_type": "LC-ESI-QTOF",
        **extra,
    }
    peaks = [(50.0 + 10 * i, 100.0 - 5 * i) for i in range(8)]
    return RawRecord(metadata=metadata, peaks=peaks)


def test_curation_reasons_are_counted():
    caffeine = MOLECULES["caffeine"]
    records = [
        _record("MSBNK-Lab-1", caffeine, 195.0877),  # conservé
        _record("MSBNK-Lab-2", caffeine, 195.0877, ion_mode="NEGATIVE"),  # mode
        _record("MSBNK-Lab-3", caffeine, 195.0877, precursor_type="[M+Na]+"),  # adduit
        _record("MSBNK-Lab-4", caffeine, 195.0877, instrument_type="LC-ESI-QQ"),
        _record("MSBNK-Lab-5", "xyz", 195.0877),  # structure
        _record("MSBNK-Lab-6", caffeine, 195.0877, inchikey="AAAAAAAAAAAAAA-X-N"),
        _record("MSBNK-Lab-7", caffeine, 200.0),  # masse
        _record("MSBNK-Lab-8", "C[N+](C)(C)C", 74.0964),  # chargée
    ]
    config = PrepareConfig(input_msp="x.msp")
    kept, report = curate(records, config)
    assert [s.identifier for s in kept] == ["MSBNK-Lab-1"]
    assert kept[0].get("contributor") == "Lab"
    assert kept[0].inchikey14 == "RYYVLZVUVIJVGH"
    for reason in (
        REASONS[1],
        REASONS[2],
        REASONS[3],
        REASONS[6],
        REASONS[7],
        REASONS[8],
        REASONS[9],
    ):
        assert report.excluded[reason] == 1, reason
    assert report.n_input == 8 and report.n_kept == 1
    assert report.parameters["input_msp"] == "x.msp"


def test_license_filter_requires_license_map():
    config = PrepareConfig(input_msp="x.msp", allowed_licenses=("CC0 1.0",))
    with pytest.raises(ValueError, match="licences"):
        curate([_record("MSBNK-A-1", MOLECULES["caffeine"], 195.0877)], config)
    kept, report = curate(
        [_record("MSBNK-A-1", MOLECULES["caffeine"], 195.0877)],
        config,
        {"MSBNK-A-1": "https://creativecommons.org/licenses/by-nc/4.0"},
    )
    assert kept == [] and report.excluded[REASONS[4]] == 1
