"""Rapport HTML autonome, construit uniquement à partir des fichiers de résultats.

Le rapport lit ``benchmark_summary.json`` et les figures PNG (intégrées en base64) : il peut
être régénéré sans recalcul et n'affiche que des valeurs réellement calculées.
"""

from __future__ import annotations

import base64
import json
import math
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from jinja2 import Environment, PackageLoader, StrictUndefined, select_autoescape

from msannot.config import METRIC_LABELS
from msannot.evaluation.identification import SETTING_LABELS
from msannot.exceptions import MsannotError
from msannot.logging_utils import get_logger

logger = get_logger(__name__)

FIGURES: dict[str, tuple[str, str]] = {
    "curation": (
        "Curation du jeu de données",
        "Spectres écartés à chaque étape, dans l'ordre d'application des filtres. La majorité "
        "des exclusions tient au périmètre choisi (mode positif, [M+H]+, haute résolution, "
        "licence ouverte), pas à des erreurs.",
    ),
    "dataset": (
        "Composition du jeu de données",
        "Masse des précurseurs, nombre de pics par spectre après nettoyage et répartition par "
        "laboratoire. Un laboratoire majoritaire peut rendre l'identification « tous "
        "laboratoires » optimiste.",
    ),
    "identification": (
        "Identification",
        "Proportion de requêtes dont le bon composé est classé premier (IC de Wilson à 95 %). "
        "Trait gris : probabilité de réussir en choisissant un candidat de même masse au hasard. "
        "À droite, seules les requêtes où la fenêtre de masse contient d'autres composés : "
        "c'est là que le spectre doit faire la différence.",
    ),
    "calibration": (
        "Calibration des scores",
        "Si l'on n'accepte le premier résultat qu'au-dessus d'un seuil de score, quelle part "
        "est correcte (gauche) et quelle part des requêtes reste annotée (droite) ? Bandes : "
        "IC de Wilson à 95 %. Cosine et modified cosine sont confondus : à précurseur égal, "
        "le modified cosine se réduit au cosine. Seuils retenant moins de 30 requêtes non tracés.",
    ),
    "analogs": (
        "Recherche d'analogues",
        "Molécule retirée de la bibliothèque : part des requêtes dont le premier résultat "
        "atteint une similarité structurale donnée. L'oracle indique le mieux possible avec "
        "cette bibliothèque ; le gris, le niveau du hasard.",
    ),
    "relation": (
        "Similarité spectrale et similarité structurale",
        "Tanimoto moyen des paires (requête, candidat) selon leur score spectral. Un score "
        "élevé annonce une structure plus proche, mais la relation est faible pour les scores "
        "moyens, et les scores élevés sont rares (échelle log à droite).",
    ),
    "network": (
        "Réseau moléculaire",
        "Arêtes du réseau spectral (modified cosine, type GNPS) colorées par similarité "
        "structurale, et distribution de cette similarité pour les arêtes et pour des paires "
        "prises au hasard.",
    ),
}
EXAMPLE_CAPTION = (
    "Spectre miroir : requête en haut, résultat en bas ; pics appariés en bleu (en orange : "
    "appariés après décalage de la différence de masse). Exemple choisi automatiquement "
    "selon une règle fixe, et non sélectionné à la main."
)


def _image(path: Path) -> str | None:
    if not path.is_file():
        logger.warning("Figure absente : %s", path)
        return None
    return "data:image/png;base64," + base64.b64encode(path.read_bytes()).decode("ascii")


def _percent(value: float | None, digits: int = 1) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    return f"{100 * value:.{digits}f} %".replace(".", ",")


def _number(value: float | int | None, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    if isinstance(value, int):
        return f"{value:,}".replace(",", " ")
    return f"{value:.{digits}f}".replace(".", ",")


def _pvalue(value: float | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    if value < 1e-300:
        return "< 10⁻³⁰⁰"
    if value < 0.001:
        mantissa, exponent = f"{value:.1e}".split("e")
        return f"{mantissa.replace('.', ',')} × 10^{int(exponent)}"
    return f"{value:.3f}".replace(".", ",")


def _label(metric: str) -> str:
    return METRIC_LABELS.get(metric, {"random": "Hasard", "oracle": "Oracle"}.get(metric, metric))


def build_report(results_dir: Path, output: Path | None = None) -> Path:
    """Génère ``report.html`` à partir d'un dossier de résultats de benchmark."""
    results_dir = results_dir.resolve()
    summary_path = results_dir / "benchmark_summary.json"
    if not summary_path.is_file():
        raise MsannotError(
            f"{summary_path} introuvable : lancez d'abord 'msannot benchmark' ou 'msannot demo'."
        )
    try:
        summary: dict[str, Any] = json.loads(summary_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise MsannotError(f"{summary_path} illisible : {exc}") from exc

    figures = {
        key: {
            "title": FIGURES.get(key, ("", ""))[0],
            "caption": FIGURES.get(key, ("", EXAMPLE_CAPTION))[1],
            "src": _image(results_dir / relative),
        }
        for key, relative in summary.get("figures", {}).items()
    }
    environment = Environment(
        loader=PackageLoader("msannot.report", "templates"),
        autoescape=select_autoescape(["html", "j2"]),
        undefined=StrictUndefined,
        trim_blocks=True,
        lstrip_blocks=True,
    )
    environment.filters.update(percent=_percent, number=_number, pvalue=_pvalue, label=_label)
    html = environment.get_template("report.html.j2").render(
        s=summary,
        figures=figures,
        example_caption=EXAMPLE_CAPTION,
        settings=SETTING_LABELS,
        generated=datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC"),
    )
    output = output or results_dir / "report.html"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html, encoding="utf-8")
    logger.info("Rapport écrit : %s", output)
    return output
