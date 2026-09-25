"""Lecture des spectres au format MGF (Mascot Generic Format), courant pour les requêtes.

Chaque spectre est délimité par ``BEGIN IONS`` / ``END IONS`` ; les métadonnées sont des
lignes ``CLÉ=valeur`` (``PEPMASS`` = m/z du précurseur) et les pics des lignes
``m/z intensité``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from msannot.exceptions import InputFormatError
from msannot.io.common import open_text
from msannot.models import Spectrum

MGF_KEYS: dict[str, str] = {
    "title": "name",
    "name": "name",
    "compound_name": "name",
    "smiles": "smiles",
    "inchikey": "inchikey",
    "ionmode": "ion_mode",
    "adduct": "precursor_type",
    "precursor_type": "precursor_type",
    "instrument_type": "instrument_type",
    "source_instrument": "instrument",
    "collision_energy": "collision_energy",
}


def load_mgf(path: Path) -> list[Spectrum]:
    """Charge les spectres d'un fichier MGF (gzip accepté).

    Raises
    ------
    InputFormatError
        Bloc non terminé, précurseur manquant, pic illisible ou fichier vide.
    """
    if not path.is_file():
        raise InputFormatError(f"fichier introuvable : {path}")
    spectra: list[Spectrum] = []
    inside = False
    metadata: dict[str, str] = {}
    peaks: list[tuple[float, float]] = []
    start = 0
    with open_text(path) as handle:
        for number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            upper = line.upper()
            if upper == "BEGIN IONS":
                inside, metadata, peaks, start = True, {}, [], number
            elif upper == "END IONS":
                if not inside:
                    raise InputFormatError(f"{path}, ligne {number} : END IONS sans BEGIN IONS")
                spectra.append(_to_spectrum(metadata, peaks, len(spectra), path, start))
                inside = False
            elif inside and "=" in line and not line[0].isdigit():
                key, value = line.split("=", 1)
                metadata[key.strip().lower()] = value.strip()
            elif inside:
                parts = line.replace(",", " ").split()
                try:
                    peaks.append((float(parts[0]), float(parts[1])))
                except (ValueError, IndexError) as exc:
                    raise InputFormatError(f"{path}, ligne {number} : pic illisible") from exc
    if inside:
        raise InputFormatError(f"{path} : bloc commencé ligne {start} non terminé (END IONS)")
    if not spectra:
        raise InputFormatError(f"aucun spectre dans {path}")
    return spectra


def _to_spectrum(
    metadata: dict[str, str], peaks: list[tuple[float, float]], index: int, path: Path, line: int
) -> Spectrum:
    try:
        precursor = float(metadata.get("pepmass", metadata.get("precursor_mz", "")).split()[0])
    except (ValueError, IndexError) as exc:
        raise InputFormatError(
            f"{path}, spectre ligne {line} : PEPMASS absent ou illisible"
        ) from exc
    identifier = metadata.get("spectrumid") or metadata.get("title") or f"requete_{index + 1}"
    normalised = {MGF_KEYS[key]: value for key, value in metadata.items() if key in MGF_KEYS}
    array = np.array(peaks, dtype=np.float64).reshape(-1, 2)
    try:
        return Spectrum(identifier, array[:, 0], array[:, 1], precursor, normalised)
    except ValueError as exc:
        raise InputFormatError(f"{path}, spectre ligne {line} : {exc}") from exc


def write_mgf(spectra: list[Spectrum], path: Path) -> None:
    """Écrit des spectres au format MGF (métadonnées principales conservées)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for spectrum in spectra:
        lines += ["BEGIN IONS", f"TITLE={spectrum.name}", f"SPECTRUMID={spectrum.identifier}"]
        lines.append(f"PEPMASS={spectrum.precursor_mz:.6f}".rstrip("0").rstrip("."))
        for key, label in (
            ("smiles", "SMILES"),
            ("inchikey", "INCHIKEY"),
            ("ion_mode", "IONMODE"),
            ("precursor_type", "ADDUCT"),
            ("instrument_type", "INSTRUMENT_TYPE"),
            ("collision_energy", "COLLISION_ENERGY"),
        ):
            if spectrum.get(key):
                lines.append(f"{label}={spectrum.get(key)}")
        lines += [
            f"{mz:.5f} {intensity:.5g}"
            for mz, intensity in zip(spectrum.mz, spectrum.intensities, strict=True)
        ]
        lines += ["END IONS", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def load_spectra(path: Path) -> list[Spectrum]:
    """Charge un fichier MSP ou MGF (détection par l'extension puis par le contenu)."""
    from msannot.io.msp import load_msp

    suffixes = [suffix.lower() for suffix in path.suffixes if suffix.lower() != ".gz"]
    if suffixes and suffixes[-1] == ".mgf":
        return load_mgf(path)
    if suffixes and suffixes[-1] == ".msp":
        return load_msp(path)
    if not path.is_file():
        raise InputFormatError(f"fichier introuvable : {path}")
    with open_text(path) as handle:
        head = handle.read(4096).upper()
    return load_mgf(path) if "BEGIN IONS" in head else load_msp(path)
