"""Utilitaires d'entrée/sortie partagés : ouverture gzip, écriture JSON/TSV déterministe."""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import math
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import asdict, is_dataclass
from datetime import date, datetime
from pathlib import Path
from typing import IO, Any

import numpy as np
import pandas as pd


def is_gzip(path: Path) -> bool:
    """Détecte un fichier gzip par sa signature (et non par son extension)."""
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: Path) -> IO[str]:
    """Ouvre un fichier texte, compressé en gzip ou non."""
    if is_gzip(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


@contextmanager
def open_output_text(path: Path) -> Iterator[IO[str]]:
    """Ouvre un fichier en écriture ; gzip déterministe (sans horodatage) si ``.gz``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with (
            path.open("wb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed,
            io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as text,
        ):
            yield text
    else:
        with path.open("w", encoding="utf-8", newline="\n") as text:
            yield text


def sha256sum(path: Path) -> str:
    """Somme de contrôle SHA-256 d'un fichier."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def to_jsonable(value: Any) -> Any:
    """Convertit récursivement un objet en structure sérialisable en JSON."""
    if is_dataclass(value) and not isinstance(value, type):
        return to_jsonable(asdict(value))
    if isinstance(value, dict):
        return {str(key): to_jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [to_jsonable(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (datetime, date)):
        return value.isoformat()
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (float, np.floating)):
        number = float(value)
        return None if math.isnan(number) or math.isinf(number) else number
    if isinstance(value, np.ndarray):
        return to_jsonable(value.tolist())
    if isinstance(value, pd.DataFrame):
        return to_jsonable(value.to_dict(orient="records"))
    if value is pd.NA:
        return None
    return value


def write_json(data: Any, path: Path) -> None:
    """Écrit un objet en JSON lisible (NaN → null)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(to_jsonable(data), indent=2, ensure_ascii=False)
    path.write_text(text + "\n", encoding="utf-8")


def read_json(path: Path) -> Any:
    """Lit un fichier JSON."""
    return json.loads(path.read_text(encoding="utf-8"))


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Écrit un tableau TSV sans index."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
