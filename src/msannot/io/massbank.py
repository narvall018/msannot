"""Accès aux versions publiques de MassBank Europe et à la licence de chaque spectre.

Les versions sont publiées sur GitHub (dépôt ``MassBank/MassBank-data``). L'export MSP ne
contient pas les licences : elles sont lues dans l'export JSON-LD (un objet ``Dataset`` par
spectre, avec ``identifier`` puis ``license``). Ce fichier fait plusieurs centaines de Mo ;
il est lu en flux, ligne par ligne, sans être chargé en mémoire.
"""

from __future__ import annotations

import re
import shutil
import urllib.request
from datetime import UTC, datetime
from pathlib import Path
from urllib.error import URLError

from msannot.exceptions import DataFetchError, InputFormatError
from msannot.io.common import open_text, sha256sum, write_json
from msannot.logging_utils import get_logger

logger = get_logger(__name__)

RELEASE_URL = "https://github.com/MassBank/MassBank-data/releases/download/{release}/{asset}"
DEFAULT_RELEASE = "2026.03"
ASSETS = ("MassBank_NISTformat.msp", "MassBank.json")
RECORD_URL = "https://massbank.eu/MassBank/RecordDisplay?id={accession}"

LICENSE_LABELS: dict[str, str] = {
    "creativecommons.org/licenses/by/4.0": "CC BY 4.0",
    "creativecommons.org/licenses/by-sa/4.0": "CC BY-SA 4.0",
    "creativecommons.org/licenses/by-nc/4.0": "CC BY-NC 4.0",
    "creativecommons.org/licenses/by-nc-sa/4.0": "CC BY-NC-SA 4.0",
    "creativecommons.org/publicdomain/zero/1.0": "CC0 1.0",
    "www.govdata.de/dl-de/by-2-0": "DL-DE BY 2.0",
}
_IDENTIFIER = re.compile(r'"identifier"\s*:\s*"([^"]+)"')
_LICENSE = re.compile(r'"license"\s*:\s*"([^"]+)"')


def license_label(url: str | None) -> str:
    """Nom court d'une licence à partir de son URL (``inconnue`` si non reconnue)."""
    if not url:
        return "inconnue"
    key = re.sub(r"^https?://", "", url.strip()).rstrip("/")
    return LICENSE_LABELS.get(key, key)


def contributor_from_accession(accession: str) -> str:
    """Contributeur MassBank à partir de l'accession (``MSBNK-Eawag-EQ000001`` → ``Eawag``)."""
    parts = accession.split("-")
    return parts[1] if len(parts) >= 3 and parts[0] == "MSBNK" else "inconnu"


def record_url(accession: str) -> str:
    """Page publique d'un spectre MassBank."""
    return RECORD_URL.format(accession=accession)


def read_license_map(path: Path) -> dict[str, str]:
    """Associe chaque accession à l'URL de sa licence, à partir de l'export JSON-LD.

    Raises
    ------
    InputFormatError
        Si aucune licence n'est trouvée (fichier inattendu).
    """
    licenses: dict[str, str] = {}
    current: str | None = None
    with open_text(path) as handle:
        for line in handle:
            identifier = _IDENTIFIER.search(line)
            if identifier:
                current = identifier.group(1)
                continue
            found = _LICENSE.search(line)
            if found and current is not None:
                licenses.setdefault(current, found.group(1))
    if not licenses:
        raise InputFormatError(f"aucune licence trouvée dans {path} : export JSON MassBank attendu")
    logger.info("Licences lues pour %d spectres", len(licenses))
    return licenses


def download(url: str, destination: Path) -> Path:
    """Télécharge un fichier en flux (écriture atomique via un fichier temporaire)."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_suffix(destination.suffix + ".part")
    try:
        with urllib.request.urlopen(url, timeout=120) as response, partial.open("wb") as out:
            shutil.copyfileobj(response, out, length=1 << 20)
    except (URLError, OSError) as exc:
        partial.unlink(missing_ok=True)
        raise DataFetchError(f"échec du téléchargement de {url} : {exc}") from exc
    partial.replace(destination)
    return destination


def fetch_release(
    outdir: Path, release: str = DEFAULT_RELEASE, *, overwrite: bool = False
) -> dict[str, Path]:
    """Télécharge l'export MSP et l'export JSON d'une version de MassBank.

    Un manifeste (URL, taille, SHA-256, date) est écrit à côté des fichiers.
    """
    if not re.fullmatch(r"\d{4}\.\d{2}(\.\d+)?", release):
        raise DataFetchError(f"version MassBank invalide : {release!r} (ex. 2026.03)")
    paths: dict[str, Path] = {}
    manifest = []
    for asset in ASSETS:
        destination = outdir / release / asset
        url = RELEASE_URL.format(release=release, asset=asset)
        if destination.exists() and not overwrite:
            logger.info("%s déjà présent : téléchargement ignoré", destination)
        else:
            logger.info("Téléchargement de %s …", url)
            download(url, destination)
        paths[asset] = destination
        manifest.append(
            {
                "asset": asset,
                "url": url,
                "bytes": destination.stat().st_size,
                "sha256": sha256sum(destination),
            }
        )
    write_json(
        {"release": release, "downloaded_utc": datetime.now(UTC), "files": manifest},
        outdir / release / "manifest.json",
    )
    return paths
