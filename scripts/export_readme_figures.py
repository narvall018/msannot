"""Copie une sélection de figures d'une exécution réelle vers docs/images (README).

Usage : ``python scripts/export_readme_figures.py [dossier_de_résultats]``
(par défaut ``results/full``, produit par ``msannot benchmark config/full.yaml``).
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SELECTION = (
    "identification",
    "calibration",
    "analogs",
    "relation",
    "network",
    "dataset",
    "curation",
    "example_1",
    "example_2",
    "example_3",
)


def main(results: Path) -> int:
    destination = ROOT / "docs" / "images"
    destination.mkdir(parents=True, exist_ok=True)
    missing = [name for name in SELECTION if not (results / "figures" / f"{name}.png").is_file()]
    if missing:
        print(f"Figures absentes dans {results} : {missing}. Lancez d'abord le benchmark.")
        return 1
    for name in SELECTION:
        shutil.copyfile(results / "figures" / f"{name}.png", destination / f"{name}.png")
        print(f"figures/{name}.png -> docs/images/{name}.png")
    return 0


if __name__ == "__main__":
    sys.exit(main(Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "results" / "full"))
