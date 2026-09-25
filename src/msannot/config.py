"""Configuration : modèles pydantic et chargement des fichiers YAML.

Les chemins relatifs d'un fichier de configuration sont résolus par rapport au dossier qui
contient ce fichier. Les valeurs par défaut sont justifiées dans ``docs/methodology.md``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, TypeVar

import yaml
from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator

from msannot.exceptions import ConfigError

Metric = Literal["cosine", "modified_cosine", "entropy"]
METRIC_LABELS: dict[str, str] = {
    "cosine": "Cosine",
    "modified_cosine": "Modified cosine",
    "entropy": "Entropie spectrale",
}
HIGH_RESOLUTION_INSTRUMENTS = (
    "LC-ESI-QTOF",
    "LC-ESI-QFT",
    "LC-ESI-ITFT",
    "LC-ESI-TOF",
    "LC-ESI-FT",
    "ESI-QTOF",
    "ESI-TOF",
)
OPEN_LICENSES = ("CC BY 4.0", "CC0 1.0")


class _StrictModel(BaseModel):
    model_config = ConfigDict(extra="forbid", frozen=True)


class CleaningConfig(_StrictModel):
    """Nettoyage des spectres (appliqué à la bibliothèque et aux requêtes)."""

    precursor_margin: float = Field(1.6, ge=0, description="Retire les pics > précurseur − marge.")
    min_relative_intensity: float = Field(0.01, ge=0, lt=1)
    max_peaks: int | None = Field(100, ge=1)
    min_peaks: int = Field(5, ge=1)


class SimilarityConfig(_StrictModel):
    """Paramètres des mesures de similarité spectrale."""

    tolerance: float = Field(0.01, gt=0, le=1, description="Tolérance sur les fragments (Da).")
    intensity_power: float = Field(0.5, gt=0, le=1)
    entropy_weighted: bool = True
    metrics: tuple[Metric, ...] = ("cosine", "modified_cosine", "entropy")

    @field_validator("metrics")
    @classmethod
    def _unique(cls, value: tuple[str, ...]) -> tuple[str, ...]:
        if not value or len(set(value)) != len(value):
            raise ValueError("metrics : liste non vide et sans doublon attendue")
        return value


class IdentificationConfig(_StrictModel):
    """Scénario d'identification (le composé est dans la bibliothèque)."""

    precursor_ppm: float = Field(10.0, gt=0, le=100)
    top_k: tuple[int, ...] = (1, 3, 10)
    max_queries: int | None = Field(None, ge=10)
    seed: int = 0


class AnalogConfig(_StrictModel):
    """Scénario de recherche d'analogues (le composé est retiré de la bibliothèque)."""

    max_mass_shift: float | None = Field(200.0, gt=0)
    n_queries: int | None = Field(None, ge=10)
    seed: int = 0
    close_analog_threshold: float = Field(0.7, gt=0, le=1)


class FingerprintConfig(_StrictModel):
    """Empreintes moléculaires de Morgan (équivalent ECFP4 pour un rayon de 2)."""

    radius: int = Field(2, ge=1, le=4)
    n_bits: int = Field(2048, ge=256)


class NetworkConfig(_StrictModel):
    """Réseau moléculaire (inspiré des paramètres par défaut de GNPS)."""

    enabled: bool = True
    metric: Metric = "modified_cosine"
    min_score: float = Field(0.7, gt=0, le=1)
    min_matches: int = Field(6, ge=1)
    top_k: int = Field(10, ge=1)


class BenchmarkConfig(_StrictModel):
    """Configuration complète d'un benchmark."""

    dataset: Path
    output_dir: Path = Path("results/benchmark")
    curation_report: Path | None = None
    cleaning: CleaningConfig = Field(default_factory=CleaningConfig)
    similarity: SimilarityConfig = Field(default_factory=SimilarityConfig)
    identification: IdentificationConfig = Field(default_factory=IdentificationConfig)
    analogs: AnalogConfig = Field(default_factory=AnalogConfig)
    fingerprint: FingerprintConfig = Field(default_factory=FingerprintConfig)
    network: NetworkConfig = Field(default_factory=NetworkConfig)
    figure_dpi: int = Field(130, ge=50, le=600)


class PrepareConfig(_StrictModel):
    """Curation d'un export MassBank (MSP + JSON des licences)."""

    input_msp: Path
    licenses_json: Path | None = None
    output: Path = Path("data/processed/massbank_curated.msp.gz")
    allowed_licenses: tuple[str, ...] | None = None
    instrument_types: tuple[str, ...] | None = HIGH_RESOLUTION_INSTRUMENTS
    precursor_types: tuple[str, ...] = ("[M+H]+",)
    ion_mode: Literal["POSITIVE", "NEGATIVE"] = "POSITIVE"
    precursor_ppm: float = Field(10.0, gt=0, le=100)
    precursor_abs_tolerance: float = Field(0.005, ge=0, le=0.1)
    cleaning: CleaningConfig = Field(default_factory=CleaningConfig)


ModelT = TypeVar("ModelT", bound=BaseModel)


def _resolve(value: Any, base_dir: Path) -> Any:
    if isinstance(value, Path) and not value.is_absolute():
        return (base_dir / value).resolve()
    return value


def load_config(path: Path, model: type[ModelT]) -> ModelT:
    """Charge et valide un fichier YAML ; résout les chemins relatifs au fichier.

    Raises
    ------
    ConfigError
        Fichier absent, YAML invalide ou valeurs incorrectes.
    """
    if not path.is_file():
        raise ConfigError(f"fichier de configuration introuvable : {path}")
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise ConfigError(f"YAML invalide dans {path} : {exc}") from exc
    if not isinstance(raw, dict):
        raise ConfigError(f"{path} doit contenir un dictionnaire YAML")
    try:
        config = model.model_validate(raw)
    except ValidationError as exc:
        raise ConfigError(f"configuration invalide ({path}) :\n{exc}") from exc
    base_dir = path.resolve().parent
    updates = {
        name: _resolve(getattr(config, name), base_dir)
        for name in type(config).model_fields
        if isinstance(getattr(config, name), Path)
    }
    return config.model_copy(update=updates)
