"""Formats MSP/MGF, licences MassBank, téléchargement (simulé)."""

import gzip
import io

import numpy as np
import pytest

from msannot.exceptions import DataFetchError, InputFormatError
from msannot.io import massbank
from msannot.io.common import to_jsonable, write_json
from msannot.io.mgf import load_mgf, load_spectra, write_mgf
from msannot.io.msp import iter_msp_records, load_msp, write_msp
from msannot.models import Spectrum

MSP = """Name: Caffeine
DB#: MSBNK-TEST-0001
SMILES: Cn1c(=O)c2c(ncn2C)n(C)c1=O
PrecursorMZ: 195.0877
Precursor_type: [M+H]+
Ion_mode: POSITIVE
Num Peaks: 3
110.0712 30
138.0662 100
123.0430 20

Name: Sans précurseur
DB#: MSBNK-TEST-0002
Num Peaks: 1
100.0 1
"""


def test_msp_reading_normalises_keys(tmp_path):
    path = tmp_path / "lib.msp"
    path.write_text(MSP)
    spectra = load_msp(path, strict=False)
    assert len(spectra) == 1
    caffeine = spectra[0]
    assert caffeine.identifier == "MSBNK-TEST-0001"
    assert caffeine.get("precursor_type") == "[M+H]+"
    assert caffeine.precursor_mz == pytest.approx(195.0877)
    assert list(caffeine.mz) == sorted(caffeine.mz)
    with pytest.raises(InputFormatError, match="précurseur"):
        load_msp(path, strict=True)


def test_msp_round_trip_and_deterministic_gzip(tmp_path):
    original = Spectrum(
        "id1",
        np.array([100.0, 150.5]),
        np.array([1.0, 0.25]),
        300.1234,
        {"name": "x", "smiles": "CCO", "license": "CC0 1.0"},
    )
    first, second = tmp_path / "a.msp.gz", tmp_path / "b.msp.gz"
    write_msp([original], first)
    write_msp([original], second)
    assert first.read_bytes() == second.read_bytes()
    (loaded,) = load_msp(first)
    assert loaded.identifier == "id1"
    assert loaded.get("license") == "CC0 1.0"
    assert np.allclose(loaded.mz, original.mz) and np.allclose(
        loaded.intensities, original.intensities
    )


def test_msp_illegible_line(tmp_path):
    path = tmp_path / "bad.msp"
    path.write_text("Name: x\nPrecursorMZ: 100\nceci n'est pas une ligne MSP\n")
    with pytest.raises(InputFormatError, match="ligne 3"):
        list(iter_msp_records(path))


def test_mgf_round_trip(tmp_path):
    spectra = [
        Spectrum(
            "q1",
            np.array([50.0, 60.0, 70.0]),
            np.array([1.0, 0.5, 0.1]),
            200.0,
            {"name": "query", "smiles": "CCO"},
        )
    ]
    path = tmp_path / "q.mgf"
    write_mgf(spectra, path)
    (loaded,) = load_spectra(path)
    assert loaded.identifier == "q1" and loaded.name == "query" and loaded.smiles == "CCO"
    assert loaded.precursor_mz == 200.0


@pytest.mark.parametrize(
    ("content", "message"),
    [
        ("BEGIN IONS\nPEPMASS=100\n50 1\n", "non terminé"),
        ("BEGIN IONS\nTITLE=x\n50 1\nEND IONS\n", "PEPMASS"),
        ("BEGIN IONS\nPEPMASS=100\n50 abc\nEND IONS\n", "pic illisible"),
        ("", "aucun spectre"),
    ],
)
def test_mgf_errors(tmp_path, content, message):
    path = tmp_path / "q.mgf"
    path.write_text(content)
    with pytest.raises(InputFormatError, match=message):
        load_mgf(path)


def test_load_spectra_detects_format_by_content(tmp_path):
    path = tmp_path / "queries.txt"
    path.write_text("BEGIN IONS\nPEPMASS=100 1500\n50 1\nEND IONS\n")
    assert load_spectra(path)[0].precursor_mz == 100.0


def test_license_map_streaming_parser(tmp_path):
    sample = """[
  {"@type": "Dataset", "identifier": "MSBNK-A-1",
   "license": "https://creativecommons.org/licenses/by/4.0/"},
  {"@type": "ChemicalSubstance", "identifier": "MSBNK-A-1"},
  {"@type": "Dataset", "identifier": "MSBNK-B-2",
   "license": "https://creativecommons.org/publicdomain/zero/1.0/"}
]"""
    path = tmp_path / "MassBank.json.gz"
    path.write_bytes(gzip.compress(sample.replace(", ", ",\n").encode()))
    licenses = massbank.read_license_map(path)
    assert {k: massbank.license_label(v) for k, v in licenses.items()} == {
        "MSBNK-A-1": "CC BY 4.0",
        "MSBNK-B-2": "CC0 1.0",
    }
    empty = tmp_path / "empty.json"
    empty.write_text("[]")
    with pytest.raises(InputFormatError):
        massbank.read_license_map(empty)


def test_massbank_helpers():
    assert massbank.contributor_from_accession("MSBNK-Eawag-EQ000001") == "Eawag"
    assert massbank.contributor_from_accession("autre") == "inconnu"
    assert massbank.license_label(None) == "inconnue"
    assert "RecordDisplay?id=MSBNK-X-1" in massbank.record_url("MSBNK-X-1")


def test_fetch_release_with_simulated_network(tmp_path, monkeypatch):
    payloads = {"MassBank_NISTformat.msp": b"Name: x\n", "MassBank.json": b"[]"}

    def fake_urlopen(url, timeout):
        return io.BytesIO(payloads[url.rsplit("/", 1)[1]])

    monkeypatch.setattr(massbank.urllib.request, "urlopen", fake_urlopen)
    paths = massbank.fetch_release(tmp_path, "2026.03")
    assert paths["MassBank.json"].read_bytes() == b"[]"
    assert (tmp_path / "2026.03" / "manifest.json").is_file()
    with pytest.raises(DataFetchError, match="invalide"):
        massbank.fetch_release(tmp_path, "latest")


def test_download_failure_is_reported(tmp_path, monkeypatch):
    def failing(url, timeout):
        raise OSError("réseau indisponible")

    monkeypatch.setattr(massbank.urllib.request, "urlopen", failing)
    with pytest.raises(DataFetchError, match="échec"):
        massbank.download("https://exemple.org/x", tmp_path / "x")
    assert not (tmp_path / "x.part").exists()


def test_json_conversion(tmp_path):
    assert to_jsonable({"a": np.float64("nan"), "b": np.int64(2), "c": np.bool_(True)}) == {
        "a": None,
        "b": 2,
        "c": True,
    }
    write_json({"x": 1}, tmp_path / "out.json")
    assert (tmp_path / "out.json").read_text().strip().startswith("{")
