"""Figures décrivant le jeu de données et sa curation."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from matplotlib.ticker import MaxNLocator

from msannot.plotting.style import (
    BLUE,
    INK_SECONDARY,
    MUTED,
    ORANGE,
    new_figure,
    set_suptitle,
    set_title,
    style_axes,
)


def dataset_figure(table: pd.DataFrame) -> Figure:
    """Masse des précurseurs, pics par spectre et spectres par laboratoire."""
    figure = new_figure(12, 3.9)
    first, second, third = figure.subplots(1, 3)
    first.hist(table["precursor_mz"], bins=40, color=BLUE, edgecolor="white", linewidth=0.5)
    first.set_xlabel("m/z du précurseur [M+H]+")
    first.set_ylabel("Spectres")
    set_title(first, "Masse des molécules")

    peaks = table["n_peaks"].clip(upper=60)
    second.hist(peaks, bins=np.arange(4.5, 61.5, 2), color=BLUE, edgecolor="white", linewidth=0.5)
    second.set_xlabel("Pics par spectre après nettoyage (60 = 60 et plus)")
    set_title(second, "Richesse des spectres")

    counts = table["contributor"].value_counts()
    top = counts.head(10)[::-1]
    third.barh(range(len(top)), top.to_numpy(), color=BLUE, height=0.6)
    third.set_yticks(range(len(top)), top.index)
    third.set_xlabel("Spectres")
    set_title(third, f"10 premiers laboratoires sur {len(counts)}")
    style_axes(first)
    style_axes(second)
    style_axes(third, grid="x")
    n_spectra = f"{len(table):,}".replace(",", " ")
    n_compounds = f"{table['inchikey14'].nunique():,}".replace(",", " ")
    set_suptitle(figure, f"Jeu de données : {n_spectra} spectres, {n_compounds} composés")
    return figure


def curation_figure(excluded: Mapping[str, int], n_input: int, n_kept: int) -> Figure:
    """Spectres exclus par étape de curation."""
    figure = new_figure(9, 4.4)
    ax = figure.subplots()
    labels = list(excluded)[::-1]
    values = [excluded[label] for label in labels]
    ax.barh(range(len(labels)), values, color=ORANGE, height=0.6, zorder=2)
    for position, value in enumerate(values):
        ax.annotate(
            f"{value:,}".replace(",", " "),
            xy=(value, position),
            xytext=(4, 0),
            textcoords="offset points",
            va="center",
            fontsize=8.5,
            color=INK_SECONDARY,
        )
    ax.set_yticks(range(len(labels)), labels)
    ax.set_xlim(0, max([*values, 1]) * 1.2)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("Spectres exclus (filtres appliqués dans cet ordre)")
    style_axes(ax, grid="x")
    ax.spines["left"].set_color(MUTED)
    set_title(
        ax,
        f"Curation : {n_kept:,} spectres conservés sur {n_input:,}".replace(",", " "),
    )
    return figure
