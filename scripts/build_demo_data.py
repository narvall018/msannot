"""Construit les données de démonstration versionnées dans ``data/demo``.

Étapes (reproductibles, aucune sélection manuelle) :

1. curation de l'export MassBank (licences **CC BY 4.0 et CC0** uniquement) ;
2. choix de trois requêtes d'exemple : pour chaque molécule de ``EXAMPLES``, le spectre
   d'accession la plus petite. Ces spectres sont retirés de la bibliothèque et écrits dans
   ``example_queries.mgf`` : la bibliothèque contient d'autres spectres des mêmes
   molécules, qu'une recherche doit retrouver ;
3. sommes SHA-256 de tous les fichiers.

Usage :
    msannot fetch --release 2026.03 --outdir data/raw
    python scripts/build_demo_data.py data/raw/2026.03
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

from msannot.config import OPEN_LICENSES, PrepareConfig
from msannot.io.massbank import read_license_map
from msannot.io.mgf import write_mgf
from msannot.io.msp import iter_msp_records
from msannot.pipeline import prepare_dataset
from msannot.processing.curation import curate

ROOT = Path(__file__).resolve().parents[1]
DEMO = ROOT / "data" / "demo"
EXAMPLES = ("Caffeine", "Carbamazepine", "Venlafaxine")


def main(release_dir: Path) -> int:
    msp = release_dir / "MassBank_NISTformat.msp"
    licenses_json = release_dir / "MassBank.json"
    config = PrepareConfig(
        input_msp=msp,
        licenses_json=licenses_json,
        output=DEMO / "massbank_open.msp.gz",
        allowed_licenses=OPEN_LICENSES,
    )
    # Première passe : curation seule, pour choisir les requêtes d'exemple.
    spectra, _ = curate(iter_msp_records(msp), config, read_license_map(licenses_json))
    queries = []
    for name in EXAMPLES:
        matches = sorted((s for s in spectra if s.name == name), key=lambda s: s.identifier)
        if not matches:
            print(f"Molécule d'exemple absente : {name}")
            return 1
        queries.append(matches[0])
    write_mgf(queries, DEMO / "example_queries.mgf")
    outcome = prepare_dataset(config, exclude=frozenset(s.identifier for s in queries))
    print(f"{outcome.report.n_kept - len(queries)} spectres → {outcome.library}")

    files = sorted(path for path in DEMO.iterdir() if path.name != "SHA256SUMS")
    lines = [f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}" for path in files]
    (DEMO / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(Path(sys.argv[1])))
