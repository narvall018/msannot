from pathlib import Path

import pytest
import yaml

from msannot.config import BenchmarkConfig, SimilarityConfig, load_config
from msannot.exceptions import ConfigError

ROOT = Path(__file__).resolve().parents[2]


@pytest.mark.parametrize("name", ["demo.yaml", "full.yaml", "default.yaml"])
def test_shipped_configurations_are_valid(name):
    config = load_config(ROOT / "config" / name, BenchmarkConfig)
    assert config.dataset.is_absolute()
    assert config.dataset.is_file()


def test_default_yaml_documents_the_defaults():
    config = load_config(ROOT / "config" / "default.yaml", BenchmarkConfig)
    defaults = BenchmarkConfig(dataset=config.dataset, output_dir=config.output_dir)
    assert config.model_dump(exclude={"output_dir"}) == defaults.model_dump(exclude={"output_dir"})


@pytest.mark.parametrize(
    "data",
    [
        {},
        {"dataset": "x.msp", "similarity": {"metrics": ["cosine", "cosine"]}},
        {"dataset": "x.msp", "similarity": {"metrics": ["tanimoto"]}},
        {"dataset": "x.msp", "similarity": {"tolerance": -1}},
        {"dataset": "x.msp", "unknown": 1},
    ],
)
def test_invalid_configurations(tmp_path, data):
    path = tmp_path / "c.yaml"
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ConfigError):
        load_config(path, BenchmarkConfig)


def test_missing_and_malformed_files(tmp_path):
    with pytest.raises(ConfigError, match="introuvable"):
        load_config(tmp_path / "absent.yaml", BenchmarkConfig)
    bad = tmp_path / "bad.yaml"
    bad.write_text("- une\n- liste\n")
    with pytest.raises(ConfigError, match="dictionnaire"):
        load_config(bad, BenchmarkConfig)


def test_similarity_defaults():
    config = SimilarityConfig()
    assert config.metrics == ("cosine", "modified_cosine", "entropy")
    assert config.tolerance == 0.01
