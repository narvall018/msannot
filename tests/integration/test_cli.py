"""CLI exécutée comme un utilisateur le ferait (données de démonstration réelles)."""

import pandas as pd
import pytest
import yaml
from tests.conftest import DEMO_LIBRARY, EXAMPLE_QUERIES
from typer.testing import CliRunner

from msannot import __version__
from msannot.cli import app
from msannot.io.mgf import load_spectra

pytestmark = pytest.mark.integration
runner = CliRunner()


def test_version():
    result = runner.invoke(app, ["--version"])
    assert result.exit_code == 0 and __version__ in result.output


def test_validate_demo_library():
    result = runner.invoke(app, ["validate", str(DEMO_LIBRARY)])
    assert result.exit_code == 0, result.output
    assert "Fichier valide" in result.output
    assert "Laboratoires" in result.output and "24" in result.output


def test_validate_reports_errors_cleanly(tmp_path):
    bad = tmp_path / "bad.mgf"
    bad.write_text("BEGIN IONS\nTITLE=x\n")
    result = runner.invoke(app, ["validate", str(bad)])
    assert result.exit_code == 1 and "Traceback" not in result.output


@pytest.mark.parametrize(
    ("mode", "metric"), [("identity", "entropy"), ("analog", "modified_cosine")]
)
def test_search_example_queries(tmp_path, mode, metric):
    result = runner.invoke(
        app,
        [
            "-q",
            "search",
            str(EXAMPLE_QUERIES),
            "--library",
            str(DEMO_LIBRARY),
            "--mode",
            mode,
            "--metric",
            metric,
            "--top",
            "3",
            "-o",
            str(tmp_path),
        ],
    )
    assert result.exit_code == 0, result.output
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t")
    expected_ids = {spectrum.identifier for spectrum in load_spectra(EXAMPLE_QUERIES)}
    assert set(hits["query_id"]) == expected_ids
    assert not set(hits["identifier"]) & expected_ids  # une requête ne se retrouve pas elle-même
    assert len(list((tmp_path / "figures").glob("*.png"))) == hits["query_id"].nunique()
    if mode == "identity":
        # les requêtes d'exemple ont d'autres spectres du même composé dans la bibliothèque
        top = hits[hits["rank"] == 1].set_index("query_id")["name"]
        assert {"Caffeine", "Carbamazepine", "Venlafaxine"} <= set(top)


def test_search_rejects_unknown_metric():
    result = runner.invoke(
        app, ["search", str(EXAMPLE_QUERIES), "-l", str(DEMO_LIBRARY), "--metric", "tanimoto"]
    )
    assert result.exit_code == 1


def test_benchmark_and_report_commands(tmp_path):
    config = tmp_path / "c.yaml"
    config.write_text(
        yaml.safe_dump(
            {
                "dataset": str(DEMO_LIBRARY),
                "output_dir": str(tmp_path / "out"),
                "identification": {"max_queries": 150},
                "analogs": {"n_queries": 40},
                "network": {"enabled": False},
            }
        )
    )
    result = runner.invoke(app, ["-q", "benchmark", str(config)])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "out" / "report.html").is_file()
    rebuilt = runner.invoke(app, ["report", str(tmp_path / "out"), "-o", str(tmp_path / "r.html")])
    assert rebuilt.exit_code == 0 and (tmp_path / "r.html").is_file()


def test_benchmark_missing_dataset(tmp_path):
    config = tmp_path / "c.yaml"
    config.write_text(yaml.safe_dump({"dataset": "absent.msp.gz"}))
    result = runner.invoke(app, ["benchmark", str(config)])
    assert result.exit_code == 1


def test_prepare_requires_licenses_for_open_only(tmp_path):
    msp = tmp_path / "x.msp"
    msp.write_text("Name: x\nPrecursorMZ: 100\nNum Peaks: 1\n50 1\n")
    result = runner.invoke(
        app, ["prepare", "--msp", str(msp), "-o", str(tmp_path / "o.msp.gz"), "--open-only"]
    )
    assert result.exit_code == 1


def test_search_min_score(tmp_path):
    result = runner.invoke(
        app,
        [
            "-q",
            "search",
            str(EXAMPLE_QUERIES),
            "-l",
            str(DEMO_LIBRARY),
            "--top",
            "10",
            "--min-score",
            "0.9",
            "-o",
            str(tmp_path),
        ],
    )
    assert result.exit_code == 0, result.output
    hits = pd.read_csv(tmp_path / "hits.tsv", sep="\t")
    assert len(hits) > 0 and (hits["score"] >= 0.9).all()
    rejected = runner.invoke(
        app, ["search", str(EXAMPLE_QUERIES), "-l", str(DEMO_LIBRARY), "--min-score", "1.5"]
    )
    assert rejected.exit_code != 0


def test_search_without_figures(tmp_path):
    result = runner.invoke(
        app,
        [
            "-q",
            "search",
            str(EXAMPLE_QUERIES),
            "-l",
            str(DEMO_LIBRARY),
            "--no-figures",
            "-o",
            str(tmp_path),
        ],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "hits.tsv").is_file()
    assert not (tmp_path / "figures").exists()


def test_search_warns_about_unknown_adduct(tmp_path):
    query = tmp_path / "sodium.mgf"
    query.write_text(
        "BEGIN IONS\nTITLE=test\nPEPMASS=217.0692\nADDUCT=[M+Na]+\n"
        "110.0712 30\n138.0662 100\n123.043 20\n83.06 10\n69.045 5\nEND IONS\n"
    )
    result = runner.invoke(
        app,
        ["-q", "search", str(query), "-l", str(DEMO_LIBRARY), "--no-figures", "-o", str(tmp_path)],
    )
    assert result.exit_code == 0, result.output
    assert "adduit [M+Na]+ absent" in result.output
