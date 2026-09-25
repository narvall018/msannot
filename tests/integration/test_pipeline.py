"""Benchmark réduit sur la bibliothèque de démonstration et préparation d'un export MassBank."""

import gzip
import json

import pandas as pd
import pytest
import yaml
from tests.conftest import DEMO_LIBRARY, MOLECULES

from msannot.config import BenchmarkConfig, PrepareConfig, load_config
from msannot.io.msp import load_msp
from msannot.pipeline import prepare_dataset, run_benchmark
from msannot.report import build_report

pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def outcome(tmp_path_factory):
    folder = tmp_path_factory.mktemp("bench")
    config_path = folder / "config.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "dataset": str(DEMO_LIBRARY),
                "output_dir": "out",
                "identification": {"max_queries": 400},
                "analogs": {"n_queries": 120},
            }
        )
    )
    return run_benchmark(load_config(config_path, BenchmarkConfig), config_path=config_path)


def test_expected_outputs(outcome):
    output = outcome.output_dir
    for name in (
        "identification_queries.tsv",
        "identification_summary.tsv",
        "identification_calibration.tsv",
        "identification_mcnemar.tsv",
        "analog_queries.tsv",
        "analog_summary.tsv",
        "relation_bins.tsv",
        "network_edges.tsv",
        "benchmark_summary.json",
        "report.html",
    ):
        assert (output / name).stat().st_size > 0, name
    figures = {path.stem for path in (output / "figures").glob("*.png")}
    assert figures >= {
        "dataset",
        "curation",
        "identification",
        "calibration",
        "analogs",
        "relation",
        "network",
        "example_1",
    }


def test_scientific_sanity(outcome):
    summary = json.loads((outcome.output_dir / "benchmark_summary.json").read_text())
    ident = pd.DataFrame(summary["identification"]["summary"])
    # cosine et modified cosine sont identiques quand les précurseurs sont confondus
    cos = ident[ident["metric"] == "cosine"].set_index(["setting", "subset"])["top1"]
    mod = ident[ident["metric"] == "modified_cosine"].set_index(["setting", "subset"])["top1"]
    assert (cos == mod).all()
    # le spectre fait mieux que la masse seule
    assert (ident["top1"] > ident["random_top1"]).all()
    analogs = pd.DataFrame(summary["analogs"]["summary"]).set_index("method")
    assert (
        analogs.loc["oracle", "median_tanimoto"]
        >= analogs.loc["modified_cosine", "median_tanimoto"]
    )
    assert (
        analogs.loc["modified_cosine", "median_tanimoto"] > analogs.loc["random", "median_tanimoto"]
    )
    network = summary["network"]
    assert network["edge_tanimoto_median"] > network["random_tanimoto_median"]
    assert "home" not in json.dumps(summary["dataset"])  # pas de chemin local


def test_report_rebuilds_from_files(outcome, tmp_path):
    html = build_report(outcome.output_dir, tmp_path / "r.html").read_text(encoding="utf-8")
    for section in (
        "Identification",
        "Recherche d'analogues",
        "Réseau moléculaire",
        "Limites",
        "Reproductibilité",
    ):
        assert section in html
    assert html.count("data:image/png;base64,") >= 8
    assert "{{" not in html


def test_prepare_dataset(tmp_path):
    caffeine = MOLECULES["caffeine"]
    block = (
        "Name: Caffeine\nDB#: MSBNK-LabA-{n}\nSMILES: {smiles}\nPrecursorMZ: 195.0877\n"
        "Precursor_type: [M+H]+\nIon_mode: POSITIVE\nSpectrum_type: MS2\n"
        "Instrument_type: LC-ESI-QTOF\nNum Peaks: 6\n"
        "110.0712 30\n138.0662 100\n123.043 20\n83.06 10\n69.045 5\n56.05 3\n\n"
    )
    msp = tmp_path / "MassBank_NISTformat.msp"
    msp.write_text("".join(block.format(n=n, smiles=caffeine) for n in range(1, 4)))
    licenses = tmp_path / "MassBank.json"
    licenses.write_text(
        '[\n{"identifier": "MSBNK-LabA-1",\n"license": "https://creativecommons.org/licenses/by/4.0/"},\n'
        '{"identifier": "MSBNK-LabA-2",\n"license": "https://creativecommons.org/licenses/by-nc/4.0"},\n'
        '{"identifier": "MSBNK-LabA-3",\n"license": "https://creativecommons.org/publicdomain/zero/1.0/"}\n]'
    )
    config = PrepareConfig(
        input_msp=msp,
        licenses_json=licenses,
        output=tmp_path / "lib.msp.gz",
        allowed_licenses=("CC BY 4.0", "CC0 1.0"),
    )
    result = prepare_dataset(config, exclude=frozenset({"MSBNK-LabA-3"}))
    spectra = load_msp(result.library)
    assert [s.identifier for s in spectra] == ["MSBNK-LabA-1"]
    assert spectra[0].get("license") == "CC BY 4.0"
    report = json.loads(result.report_path.read_text())
    assert report["excluded"]["licence non retenue"] == 1
    assert report["excluded_identifiers"] == ["MSBNK-LabA-3"]
    assert "/" not in report["parameters"]["input_msp"]
    with gzip.open(result.attribution_path, "rt") as handle:
        assert "MSBNK-LabA-1" in handle.read()
