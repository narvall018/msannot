"""Lecture et écriture de bibliothèques spectrales au format MSP (NIST).

Le format MSP est un format texte : un bloc par spectre, des lignes ``Clé: valeur`` puis
``Num Peaks: n`` suivi de ``n`` lignes ``m/z intensité``, et une ligne vide entre deux
spectres. Les clés sont normalisées (minuscules, noms canoniques) pour que les fichiers
de différentes sources soient traités de la même façon.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from msannot.exceptions import InputFormatError
from msannot.io.common import open_output_text, open_text
from msannot.logging_utils import get_logger
from msannot.models import Spectrum

logger = get_logger(__name__)

# Clé MSP (en minuscules) → clé de métadonnée normalisée.
KEY_ALIASES: dict[str, str] = {
    "name": "name",
    "db#": "identifier",
    "accession": "identifier",
    "spectrumid": "identifier",
    "inchikey": "inchikey",
    "smiles": "smiles",
    "inchi": "inchi",
    "precursormz": "precursor_mz",
    "precursor_mz": "precursor_mz",
    "precursor_type": "precursor_type",
    "precursortype": "precursor_type",
    "adduct": "precursor_type",
    "spectrum_type": "spectrum_type",
    "ion_mode": "ion_mode",
    "ionmode": "ion_mode",
    "instrument_type": "instrument_type",
    "instrument": "instrument",
    "collision_energy": "collision_energy",
    "formula": "formula",
    "contributor": "contributor",
    "license": "license",
    "source": "source",
    "inchikey14": "inchikey14",
    "comments": "comments",
}
WRITE_ORDER = (
    ("Name", "name"),
    ("DB#", None),
    ("InChIKey", "inchikey"),
    ("SMILES", "smiles"),
    ("Formula", "formula"),
    ("PrecursorMZ", None),
    ("Precursor_type", "precursor_type"),
    ("Spectrum_type", "spectrum_type"),
    ("Ion_mode", "ion_mode"),
    ("Instrument_type", "instrument_type"),
    ("Instrument", "instrument"),
    ("Collision_energy", "collision_energy"),
    ("Contributor", "contributor"),
    ("License", "license"),
    ("Source", "source"),
)
_PEAK_LINE = re.compile(r"^\s*([0-9.]+(?:[eE][-+]?\d+)?)[\s,;]+([0-9.]+(?:[eE][-+]?\d+)?)")


@dataclass(slots=True)
class RawRecord:
    """Bloc MSP brut : métadonnées normalisées et pics (avant validation)."""

    metadata: dict[str, str] = field(default_factory=dict)
    peaks: list[tuple[float, float]] = field(default_factory=list)
    line: int = 0


def _normalise_key(key: str) -> str:
    cleaned = key.strip().lower().replace(" ", "_")
    return KEY_ALIASES.get(cleaned, cleaned)


def iter_msp_records(path: Path) -> Iterator[RawRecord]:
    """Parcourt un fichier MSP bloc par bloc (lecture en flux, gzip accepté).

    Raises
    ------
    InputFormatError
        Fichier absent ou ligne de pic illisible.
    """
    if not path.is_file():
        raise InputFormatError(f"fichier introuvable : {path}")
    record = RawRecord()
    with open_text(path) as handle:
        for number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                if record.metadata or record.peaks:
                    yield record
                record = RawRecord(line=number + 1)
                continue
            if not record.metadata and not record.peaks:
                record.line = number
            match = _PEAK_LINE.match(line)
            if match and (record.peaks or "num_peaks" in record.metadata):
                record.peaks.append((float(match.group(1)), float(match.group(2))))
            elif ":" in line:
                key, value = line.split(":", 1)
                record.metadata.setdefault(_normalise_key(key), value.strip())
            else:
                raise InputFormatError(
                    f"{path}, ligne {number} : ligne MSP illisible ({line[:60]!r})"
                )
    if record.metadata or record.peaks:
        yield record


def record_to_spectrum(record: RawRecord, index: int) -> Spectrum:
    """Convertit un bloc brut en :class:`Spectrum`.

    Raises
    ------
    ValueError
        Précurseur absent ou invalide, pics incohérents.
    """
    metadata = dict(record.metadata)
    identifier = metadata.pop("identifier", "") or f"spectrum_{index}"
    precursor = metadata.pop("precursor_mz", "")
    try:
        precursor_mz = float(precursor.split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{identifier} : m/z du précurseur absent ou illisible") from exc
    metadata.pop("num_peaks", None)
    peaks = np.array(record.peaks, dtype=np.float64).reshape(-1, 2)
    return Spectrum(identifier, peaks[:, 0], peaks[:, 1], precursor_mz, metadata)


def load_msp(path: Path, *, strict: bool = True) -> list[Spectrum]:
    """Charge tous les spectres d'un fichier MSP.

    Parameters
    ----------
    path
        Fichier MSP (gzip accepté).
    strict
        Si vrai, un spectre invalide provoque une erreur ; sinon il est ignoré et compté.
    """
    spectra: list[Spectrum] = []
    skipped = 0
    for index, record in enumerate(iter_msp_records(path)):
        try:
            spectra.append(record_to_spectrum(record, index))
        except ValueError as exc:
            if strict:
                raise InputFormatError(f"{path}, ligne {record.line} : {exc}") from exc
            skipped += 1
    if skipped:
        logger.warning("%d spectre(s) invalide(s) ignoré(s) dans %s", skipped, path.name)
    if not spectra:
        raise InputFormatError(f"aucun spectre lisible dans {path}")
    return spectra


def _format_number(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def write_msp(spectra: Iterable[Spectrum], path: Path) -> int:
    """Écrit des spectres au format MSP (gzip déterministe si ``.gz``). Renvoie leur nombre."""
    count = 0
    with open_output_text(path) as handle:
        for spectrum in spectra:
            for label, key in WRITE_ORDER:
                if label == "DB#":
                    value = spectrum.identifier
                elif label == "PrecursorMZ":
                    value = _format_number(spectrum.precursor_mz)
                else:
                    value = spectrum.get(key or "")
                if value:
                    handle.write(f"{label}: {value}\n")
            handle.write(f"Num Peaks: {spectrum.n_peaks}\n")
            for mz, intensity in zip(spectrum.mz, spectrum.intensities, strict=True):
                handle.write(f"{mz:.5f} {intensity:.5g}\n")
            handle.write("\n")
            count += 1
    return count
