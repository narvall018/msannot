"""Spectres miroirs : requête en haut, spectre de bibliothèque en bas, pics appariés en couleur.

C'est la représentation standard pour juger visuellement une annotation. Les structures
(dessinées par RDKit) sont affichées à côté, avec le score et la similarité de Tanimoto.
"""

from __future__ import annotations

import io

import matplotlib.image as mpimg
import numpy as np
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

from msannot.chem.structures import molecule_png
from msannot.config import METRIC_LABELS, Metric, SimilarityConfig
from msannot.models import Spectrum
from msannot.plotting.style import (
    COLOR_MATCHED,
    COLOR_SHIFTED,
    INK,
    INK_SECONDARY,
    MUTED,
    new_figure,
    set_suptitle,
    style_axes,
)
from msannot.similarity.pairwise import matched_peaks, similarity


def _stems(ax: Axes, mz: np.ndarray, heights: np.ndarray, colors: list[str]) -> None:
    ax.vlines(mz, 0, heights, colors=colors, linewidth=1.3)


def _structure(ax: Axes, smiles: str, caption: str) -> None:
    ax.set_axis_off()
    image = molecule_png(smiles, 360, 260)
    if image is not None:
        ax.imshow(mpimg.imread(io.BytesIO(image), format="png"))
    ax.set_title(caption, fontsize=8.5, color=INK_SECONDARY, wrap=True)


def mirror_figure(
    query: Spectrum,
    reference: Spectrum,
    metric: Metric,
    config: SimilarityConfig,
    *,
    title: str,
    tanimoto: float | None = None,
) -> Figure:
    """Spectre miroir requête/bibliothèque avec structures et score."""
    figure = new_figure(11, 4.8)
    grid = figure.add_gridspec(2, 2, width_ratios=[3.2, 1.0])
    ax = figure.add_subplot(grid[:, 0])
    matches = matched_peaks(metric, reference, query, config)
    result = similarity(metric, reference, query, config)

    query_colors = [MUTED] * query.n_peaks
    reference_colors = [MUTED] * reference.n_peaks
    for ref_idx, query_idx, shifted in zip(
        matches.idx_a, matches.idx_b, matches.shifted, strict=True
    ):
        color = COLOR_SHIFTED if shifted else COLOR_MATCHED
        query_colors[int(query_idx)] = color
        reference_colors[int(ref_idx)] = color
    _stems(ax, query.mz, 100 * query.intensities / query.intensities.max(), query_colors)
    _stems(
        ax,
        reference.mz,
        -100 * reference.intensities / reference.intensities.max(),
        reference_colors,
    )
    ax.axhline(0, color=INK_SECONDARY, linewidth=0.8)
    ax.set_ylim(-112, 112)
    ax.set_yticks([-100, -50, 0, 50, 100], ["100", "50", "0", "50", "100"])
    ax.set_ylabel("Intensité relative (%)")
    ax.set_xlabel("m/z")
    ax.text(
        0.01,
        0.97,
        f"Requête — {query.identifier}",
        transform=ax.transAxes,
        va="top",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    ax.text(
        0.01,
        0.03,
        f"Bibliothèque — {reference.identifier}",
        transform=ax.transAxes,
        va="bottom",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    handles = [Line2D([0], [0], color=COLOR_MATCHED, linewidth=2, label="Pic apparié")]
    if matches.shifted.any():
        handles.append(
            Line2D([0], [0], color=COLOR_SHIFTED, linewidth=2, label="Apparié par décalage")
        )
    handles.append(Line2D([0], [0], color=MUTED, linewidth=2, label="Non apparié"))
    ax.legend(
        handles=handles, frameon=False, fontsize=8, loc="upper right", labelcolor=INK_SECONDARY
    )
    style_axes(ax, grid=None)
    summary = (
        f"{METRIC_LABELS[metric]} = {result.score:.3f} · {result.n_matches} pics appariés"
        + (f" · Tanimoto = {tanimoto:.2f}" if tanimoto is not None else "")
        + f" · Δ précurseur = {reference.precursor_mz - query.precursor_mz:+.4f}"
    )
    ax.set_title(summary, loc="left", fontsize=9, color=INK)
    _structure(figure.add_subplot(grid[0, 1]), query.smiles, f"Requête : {query.name[:40]}")
    _structure(
        figure.add_subplot(grid[1, 1]), reference.smiles, f"Résultat : {reference.name[:40]}"
    )
    set_suptitle(figure, title)
    return figure
