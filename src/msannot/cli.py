"""Interface en ligne de commande (Typer).

Aucune logique scientifique ici : la CLI lit les arguments, appelle le package et affiche
les résultats. Les erreurs attendues (:class:`MsannotError`) donnent un message lisible et
un code de sortie 1.
"""

from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Annotated, cast

import numpy as np
import typer
from rich.console import Console
from rich.table import Table

from msannot import __version__
from msannot.config import (
    HIGH_RESOLUTION_INSTRUMENTS,
    OPEN_LICENSES,
    BenchmarkConfig,
    CleaningConfig,
    Metric,
    PrepareConfig,
    SimilarityConfig,
    load_config,
)
from msannot.exceptions import ConfigError, MsannotError
from msannot.logging_utils import setup_logging

app = typer.Typer(
    name="msannot",
    help="Annotation de petites molécules par spectres MS/MS et évaluation rigoureuse.",
    no_args_is_help=True,
    add_completion=False,
    rich_markup_mode="rich",
)
console = Console()
error_console = Console(stderr=True)
DEMO_CONFIG = Path("config") / "demo.yaml"


@contextmanager
def _handle_errors() -> Iterator[None]:
    try:
        yield
    except MsannotError as exc:
        error_console.print(f"[bold red]Erreur :[/] {exc}")
        raise typer.Exit(code=1) from exc


def _version_callback(value: bool) -> None:
    if value:
        console.print(f"msannot {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    verbose: Annotated[bool, typer.Option("--verbose", "-v", help="Messages de débogage.")] = False,
    quiet: Annotated[bool, typer.Option("--quiet", "-q", help="Avertissements seulement.")] = False,
    log_file: Annotated[Path | None, typer.Option(help="Copie le journal dans ce fichier.")] = None,
    version: Annotated[
        bool, typer.Option("--version", callback=_version_callback, is_eager=True, help="Version.")
    ] = False,
) -> None:
    """Recherche en bibliothèque spectrale, analogues structuraux et évaluation."""
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    setup_logging(level, log_file)


# ---------------------------------------------------------------------------- données
@app.command()
def fetch(
    release: Annotated[str, typer.Option(help="Version de MassBank (ex. 2026.03).")] = "2026.03",
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Dossier de sortie.")] = Path(
        "data/raw"
    ),
    force: Annotated[bool, typer.Option(help="Retélécharge les fichiers existants.")] = False,
) -> None:
    """Télécharge l'export MSP et l'export JSON (licences) d'une version de MassBank."""
    from msannot.io.massbank import fetch_release

    with _handle_errors():
        paths = fetch_release(outdir, release, overwrite=force)
    for path in paths.values():
        console.print(f"[green]✓[/] {path} ({path.stat().st_size / 1e6:.1f} Mo)")


@app.command()
def prepare(
    msp: Annotated[Path, typer.Option(help="Export MSP de MassBank.")],
    output: Annotated[
        Path, typer.Option("--output", "-o", help="Bibliothèque produite (.msp.gz).")
    ],
    licenses: Annotated[Path | None, typer.Option(help="Export JSON MassBank (licences).")] = None,
    open_only: Annotated[
        bool, typer.Option(help=f"Ne garde que les licences {', '.join(OPEN_LICENSES)}.")
    ] = False,
    all_instruments: Annotated[
        bool, typer.Option(help="Garde aussi les instruments basse résolution.")
    ] = False,
) -> None:
    """Cure un export MassBank : filtres, contrôles de structure, nettoyage des spectres."""
    from msannot.pipeline import prepare_dataset

    with _handle_errors():
        try:
            config = PrepareConfig(
                input_msp=msp,
                licenses_json=licenses,
                output=output,
                allowed_licenses=OPEN_LICENSES if open_only else None,
                instrument_types=None if all_instruments else HIGH_RESOLUTION_INSTRUMENTS,
            )
        except ValueError as exc:
            raise ConfigError(str(exc)) from exc
        if open_only and licenses is None:
            raise ConfigError("--open-only nécessite --licenses (export JSON de MassBank)")
        outcome = prepare_dataset(config)
    report = outcome.report
    table = Table(title="Curation", show_header=False)
    for reason, count in report.excluded.items():
        table.add_row(reason, str(count))
    console.print(table)
    console.print(
        f"[green]✓[/] {report.n_kept} spectres ({report.n_compounds} composés) → {outcome.library}"
    )


@app.command()
def validate(
    file: Annotated[Path, typer.Argument(help="Fichier MSP ou MGF (.gz accepté).")],
) -> None:
    """Vérifie un fichier de spectres et résume son contenu."""
    from msannot.io.mgf import load_spectra

    with _handle_errors():
        spectra = load_spectra(file)
    peaks = np.array([s.n_peaks for s in spectra])
    precursors = np.array([s.precursor_mz for s in spectra])
    table = Table(title=f"Validation de {file.name}", show_header=False)
    table.add_row("Spectres", str(len(spectra)))
    table.add_row("Avec structure (SMILES)", str(sum(bool(s.smiles) for s in spectra)))
    table.add_row("Composés (InChIKey, 1er bloc)", str(len({s.inchikey14 for s in spectra} - {""})))
    contributors = {s.get("contributor") for s in spectra} - {""}
    table.add_row("Laboratoires (contributeurs)", str(len(contributors)) if contributors else "—")
    table.add_row("m/z du précurseur", f"{precursors.min():.2f} – {precursors.max():.2f}")
    table.add_row("Pics par spectre (médiane)", f"{np.median(peaks):.0f}")
    for label, key in (
        ("Modes", "ion_mode"),
        ("Adduits", "precursor_type"),
        ("Licences", "license"),
    ):
        counts = Counter(s.get(key) or "—" for s in spectra).most_common(4)
        table.add_row(label, ", ".join(f"{k} ({v})" for k, v in counts))
    console.print(table)
    console.print("[green]✓ Fichier valide.[/]")


# ---------------------------------------------------------------------------- recherche
@app.command()
def search(
    queries: Annotated[Path, typer.Argument(help="Spectres requêtes (MGF ou MSP).")],
    library: Annotated[Path, typer.Option("--library", "-l", help="Bibliothèque MSP.")],
    mode: Annotated[str, typer.Option(help="identity ou analog.")] = "identity",
    metric: Annotated[str, typer.Option(help="cosine, modified_cosine ou entropy.")] = "entropy",
    top: Annotated[int, typer.Option(help="Résultats par requête.")] = 5,
    ppm: Annotated[float, typer.Option(help="Fenêtre de précurseur (identification).")] = 10.0,
    max_shift: Annotated[
        float, typer.Option(help="Écart de masse maximal (analogues, Da).")
    ] = 200.0,
    outdir: Annotated[Path, typer.Option("--outdir", "-o", help="Dossier de sortie.")] = Path(
        "results/search"
    ),
    min_score: Annotated[
        float,
        typer.Option(
            help="Score minimal des résultats conservés (voir la calibration).", min=0, max=1
        ),
    ] = 0.0,
) -> None:
    """Recherche des spectres dans une bibliothèque ; écrit résultats et spectres miroirs."""
    import pandas as pd

    from msannot.io.common import write_tsv
    from msannot.io.mgf import load_spectra
    from msannot.io.msp import load_msp
    from msannot.plotting import mirror_figure, save_figure
    from msannot.processing.cleaning import clean_all, clean_spectrum
    from msannot.search.library import SpectralLibrary

    if mode not in ("identity", "analog"):
        error_console.print("[bold red]Erreur :[/] --mode doit valoir identity ou analog")
        raise typer.Exit(1)
    if metric not in ("cosine", "modified_cosine", "entropy"):
        error_console.print("[bold red]Erreur :[/] --metric inconnue")
        raise typer.Exit(1)
    chosen_metric = cast(Metric, metric)
    cleaning = CleaningConfig()
    similarity = SimilarityConfig()
    with _handle_errors():
        library_spectra, _ = clean_all(load_msp(library), cleaning)
        spectral_library = SpectralLibrary(library_spectra, similarity)
        query_spectra = load_spectra(queries)
    all_hits = []
    for query in query_spectra:
        cleaned = clean_spectrum(query, cleaning)
        if cleaned is None:
            console.print(f"[yellow]⚠ {query.identifier} : trop peu de pics après nettoyage[/]")
            continue
        exclude = spectral_library.table["identifier"].to_numpy() == query.identifier
        hits = spectral_library.search(
            cleaned,
            metric=chosen_metric,
            mode="identity" if mode == "identity" else "analog",
            top_k=top,
            ppm=ppm,
            max_mass_shift=max_shift,
            exclude=exclude,
        )
        hits = hits[hits["score"] >= min_score].reset_index(drop=True)
        hits.insert(0, "query_id", query.identifier)
        all_hits.append(hits)
        table = Table(title=f"{query.identifier} — {query.name} (m/z {query.precursor_mz:.4f})")
        for column in ("rank", "name", "score", "n_matches", "delta_mz", "contributor"):
            table.add_column(column)
        for _, row in hits.iterrows():
            table.add_row(
                str(row["rank"]),
                str(row["name"])[:50],
                f"{row['score']:.3f}",
                str(row["n_matches"]),
                f"{row['delta_mz']:+.4f}",
                str(row["contributor"]),
            )
        console.print(table if len(hits) else f"{query.identifier} : aucun résultat")
        if len(hits):
            best = spectral_library.spectra[int(hits.iloc[0]["library_index"])]
            figure = mirror_figure(
                cleaned, best, chosen_metric, similarity, title=f"{query.name} — meilleur résultat"
            )
            safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in query.identifier)
            save_figure(figure, outdir / "figures" / f"{safe}.png")
    if all_hits:
        write_tsv(pd.concat(all_hits, ignore_index=True), outdir / "hits.tsv")
        console.print(f"[green]✓[/] Résultats : {outdir / 'hits.tsv'}")


# ---------------------------------------------------------------------------- benchmark
def _run(config_path: Path, output: Path | None, report: bool) -> None:
    from msannot.pipeline import run_benchmark

    with _handle_errors():
        config = load_config(config_path, BenchmarkConfig)
        if not config.dataset.is_file():
            raise ConfigError(f"jeu de données introuvable : {config.dataset}")
        outcome = run_benchmark(
            config, output_dir=output, config_path=config_path, make_report=report
        )
    rows = [
        r
        for r in outcome.summary["identification"]["summary"].to_dict(orient="records")
        if r["subset"] == "avec concurrents"
    ]
    table = Table(title="Identification — top-1 (requêtes avec concurrents)")
    for column in ("scénario", "mesure", "n", "top-1", "hasard"):
        table.add_column(column)
    for row in rows:
        table.add_row(
            row["setting"],
            row["metric"],
            str(row["n_queries"]),
            f"{row['top1']:.3f}",
            f"{row['random_top1']:.3f}",
        )
    console.print(table)
    analog = outcome.summary["analogs"]["summary"]
    table = Table(title="Analogues — Tanimoto du 1er résultat")
    for column in ("méthode", "médiane", "proches"):
        table.add_column(column)
    for _, row in analog.iterrows():
        table.add_row(
            str(row["method"]), f"{row['median_tanimoto']:.3f}", f"{row['close_fraction']:.3f}"
        )
    console.print(table)
    console.print(f"[green]✓[/] Résultats : {outcome.output_dir}")
    if outcome.report_path:
        console.print(f"[green]✓[/] Rapport : {outcome.report_path}")


@app.command()
def benchmark(
    config: Annotated[Path, typer.Argument(help="Fichier de configuration YAML.")],
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Dossier de sortie.")
    ] = None,
    report: Annotated[bool, typer.Option(help="Génère le rapport HTML.")] = True,
) -> None:
    """Évalue identification, analogues, relation spectre/structure et réseau moléculaire."""
    _run(config, output, report)


@app.command()
def report(
    results: Annotated[Path, typer.Argument(help="Dossier de résultats d'un benchmark.")],
    output: Annotated[Path | None, typer.Option("--output", "-o", help="Fichier HTML.")] = None,
) -> None:
    """Régénère le rapport HTML sans recalcul."""
    from msannot.report.builder import build_report

    with _handle_errors():
        path = build_report(results, output)
    console.print(f"[green]✓[/] Rapport : {path}")


def find_demo_config(start: Path | None = None) -> Path:
    """Cherche ``config/demo.yaml`` depuis le dossier courant puis ses parents."""
    here = (start or Path.cwd()).resolve()
    for folder in (here, *here.parents):
        if (folder / DEMO_CONFIG).is_file():
            return folder / DEMO_CONFIG
    source_tree = Path(__file__).resolve().parents[2] / DEMO_CONFIG
    if source_tree.is_file():
        return source_tree
    raise MsannotError("config/demo.yaml introuvable : lancez la commande depuis le dépôt msannot.")


@app.command()
def demo(
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Dossier de sortie.")
    ] = None,
) -> None:
    """Lance le benchmark de démonstration sur la bibliothèque fournie (data/demo)."""
    with _handle_errors():
        config_path = find_demo_config()
    console.print(f"Configuration de démonstration : {config_path}")
    _run(config_path, output, report=True)
