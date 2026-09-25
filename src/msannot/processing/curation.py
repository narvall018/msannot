"""Curation d'une bibliothèque spectrale publique (MassBank).

Une bibliothèque publique agrège des spectres de nombreux laboratoires : formats de
métadonnées hétérogènes, erreurs de structure, précurseurs incohérents. Chaque filtre
ci-dessous est explicite et compté, afin que la composition du jeu de données final soit
traçable (voir ``docs/data.md``).
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from msannot.chem.structures import parse_structure, ppm_error, protonated_mz
from msannot.config import PrepareConfig
from msannot.io.massbank import contributor_from_accession, license_label, record_url
from msannot.io.msp import RawRecord, record_to_spectrum
from msannot.logging_utils import get_logger
from msannot.models import Spectrum
from msannot.processing.cleaning import clean_spectrum

logger = get_logger(__name__)

REASONS = (
    "type de spectre différent de MS2",
    "mode d'ionisation non retenu",
    "adduit non retenu",
    "instrument non retenu",
    "licence non retenue",
    "précurseur absent ou illisible",
    "structure absente ou invalide",
    "InChIKey incohérent avec le SMILES",
    "molécule chargée",
    "masse du précurseur incohérente",
    "trop peu de pics après nettoyage",
)


@dataclass(slots=True)
class CurationReport:
    """Bilan de la curation : effectifs exclus par raison et composition finale."""

    n_input: int = 0
    excluded: dict[str, int] = field(default_factory=lambda: dict.fromkeys(REASONS, 0))
    n_kept: int = 0
    n_compounds: int = 0
    licenses: dict[str, int] = field(default_factory=dict)
    contributors: dict[str, int] = field(default_factory=dict)
    instrument_types: dict[str, int] = field(default_factory=dict)
    median_abs_ppm: float | None = None
    parameters: dict[str, Any] = field(default_factory=dict)

    def exclude(self, reason: str) -> None:
        """Compte un spectre exclu pour ``reason``."""
        self.excluded[reason] += 1


def curate(
    records: Iterable[RawRecord], config: PrepareConfig, licenses: dict[str, str] | None = None
) -> tuple[list[Spectrum], CurationReport]:
    """Filtre et normalise des spectres bruts.

    Parameters
    ----------
    records
        Blocs MSP bruts (:func:`msannot.io.msp.iter_msp_records`).
    config
        Critères de curation et de nettoyage.
    licenses
        Licence (URL) de chaque accession ; obligatoire si ``allowed_licenses`` est défini.

    Returns
    -------
    tuple
        Spectres conservés (métadonnées enrichies) et bilan de curation.
    """
    parameters = config.model_dump(mode="json", exclude={"input_msp", "licenses_json", "output"})
    parameters["input_msp"] = config.input_msp.name  # jamais de chemin local dans un bilan
    parameters["licenses_json"] = config.licenses_json.name if config.licenses_json else None
    report = CurationReport(parameters=parameters)
    kept: list[Spectrum] = []
    errors: list[float] = []
    instruments = set(config.instrument_types) if config.instrument_types else None
    allowed_licenses = set(config.allowed_licenses) if config.allowed_licenses else None
    if allowed_licenses and licenses is None:
        raise ValueError("un fichier de licences est nécessaire pour filtrer par licence")

    for index, record in enumerate(records):
        report.n_input += 1
        meta = record.metadata
        accession = meta.get("identifier", f"spectrum_{index}")
        if meta.get("spectrum_type", "MS2").upper() != "MS2":
            report.exclude(REASONS[0])
            continue
        if meta.get("ion_mode", "").upper() != config.ion_mode:
            report.exclude(REASONS[1])
            continue
        if meta.get("precursor_type", "") not in config.precursor_types:
            report.exclude(REASONS[2])
            continue
        if instruments is not None and meta.get("instrument_type", "") not in instruments:
            report.exclude(REASONS[3])
            continue
        license_name = license_label(licenses.get(accession)) if licenses is not None else ""
        if allowed_licenses is not None and license_name not in allowed_licenses:
            report.exclude(REASONS[4])
            continue
        try:
            spectrum = record_to_spectrum(record, index)
        except ValueError:
            report.exclude(REASONS[5])
            continue
        structure = parse_structure(meta.get("smiles", ""))
        if structure is None:
            report.exclude(REASONS[6])
            continue
        declared = meta.get("inchikey", "")
        if declared and declared[:14] != structure.inchikey14:
            report.exclude(REASONS[7])
            continue
        if structure.formal_charge != 0:
            report.exclude(REASONS[8])
            continue
        expected = protonated_mz(structure)
        error = ppm_error(spectrum.precursor_mz, expected)
        if (
            abs(error) > config.precursor_ppm
            and abs(spectrum.precursor_mz - expected) > config.precursor_abs_tolerance
        ):
            report.exclude(REASONS[9])
            continue
        cleaned = clean_spectrum(spectrum, config.cleaning)
        if cleaned is None:
            report.exclude(REASONS[10])
            continue
        errors.append(error)
        kept.append(
            cleaned.with_metadata(
                smiles=structure.canonical_smiles,
                inchikey=structure.inchikey,
                inchikey14=structure.inchikey14,
                formula=structure.formula,
                contributor=contributor_from_accession(accession),
                license=license_name,
                source=record_url(accession) if accession.startswith("MSBNK-") else "",
            )
        )

    report.n_kept = len(kept)
    report.n_compounds = len({spectrum.inchikey14 for spectrum in kept})
    report.licenses = dict(
        Counter(s.get("license") or "non renseignée" for s in kept).most_common()
    )
    report.contributors = dict(Counter(s.contributor for s in kept).most_common())
    report.instrument_types = dict(Counter(s.get("instrument_type") for s in kept).most_common())
    report.median_abs_ppm = float(np.median(np.abs(errors))) if errors else None
    logger.info(
        "Curation : %d spectres conservés sur %d (%d composés)",
        report.n_kept,
        report.n_input,
        report.n_compounds,
    )
    return kept, report
