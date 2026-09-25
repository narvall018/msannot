"""Figures d'évaluation : identification, calibration, analogues, spectre/structure."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from msannot.config import METRIC_LABELS
from msannot.evaluation.identification import SETTING_LABELS
from msannot.plotting.style import (
    COLOR_ORACLE,
    COLOR_RANDOM,
    INK_SECONDARY,
    LINE_WIDTH,
    METRIC_COLORS,
    SURFACE,
    new_figure,
    set_suptitle,
    set_title,
    style_axes,
    style_legend,
)


def identification_figure(summary: pd.DataFrame, metrics: Sequence[str]) -> Figure:
    """Exactitude top-1 par mesure et scénario, avec IC de Wilson et référence aléatoire."""
    subsets = list(dict.fromkeys(summary["subset"]))
    figure = new_figure(11, 4.4)
    axes = figure.subplots(1, len(subsets), sharey=True, squeeze=False)[0]
    settings = list(SETTING_LABELS)
    width = 0.8 / len(metrics)
    for ax, subset in zip(axes, subsets, strict=True):
        data = summary[summary["subset"] == subset]
        for position, setting in enumerate(settings):
            rows = data[data["setting"] == setting].set_index("metric")
            if rows.empty:
                continue
            for offset, metric in enumerate(metrics):
                row = rows.loc[metric]
                x = position - 0.4 + width * (offset + 0.5)
                ax.bar(
                    x,
                    row["top1"],
                    width=width * 0.85,
                    color=METRIC_COLORS[metric],
                    zorder=2,
                    label=METRIC_LABELS[metric] if position == 0 else None,
                )
                ax.errorbar(
                    x,
                    row["top1"],
                    yerr=[[row["top1"] - row["top1_ci_low"]], [row["top1_ci_high"] - row["top1"]]],
                    fmt="none",
                    ecolor=INK_SECONDARY,
                    elinewidth=1,
                    capsize=2.5,
                )
            random = float(rows["random_top1"].iloc[0])
            ax.hlines(
                random,
                position - 0.42,
                position + 0.42,
                color=COLOR_RANDOM,
                linewidth=2,
                label="Candidat au hasard (masse seule)" if position == 0 else None,
                zorder=3,
            )
            ax.annotate(
                f"n = {int(rows['n_queries'].iloc[0])}",
                xy=(position, 0),
                xytext=(0, -28),
                textcoords="offset points",
                ha="center",
                fontsize=8,
                color=INK_SECONDARY,
                annotation_clip=False,
            )
        ax.set_xticks(range(len(settings)), [SETTING_LABELS[s] for s in settings])
        ax.set_ylim(0, 1.05)
        style_axes(ax)
        set_title(ax, f"Requêtes : {subset}")
    axes[0].set_ylabel("Exactitude top-1 (bon composé en 1re position)")
    handles, labels = axes[0].get_legend_handles_labels()
    figure.legend(
        handles,
        labels,
        loc="outside lower center",
        ncols=len(labels),
        frameon=False,
        fontsize=8.5,
        labelcolor=INK_SECONDARY,
    )
    set_suptitle(figure, "Identification : le bon composé est-il classé premier ?")
    return figure


MIN_RETAINED = 30  # seuils trop sélectifs (moins de 30 requêtes) non tracés


def calibration_figure(calibration: pd.DataFrame, metrics: Sequence[str], setting: str) -> Figure:
    """Précision et couverture du premier résultat en fonction d'un seuil de score.

    Les seuils retenant moins de 30 requêtes ne sont pas tracés : la précision y est trop
    incertaine pour être lue.
    """
    data = calibration[calibration["setting"] == setting]
    figure = new_figure(11, 4.2)
    left, right = figure.subplots(1, 2, sharex=True)
    for metric in metrics:
        rows = data[(data["metric"] == metric) & (data["n_retained"] >= MIN_RETAINED)]
        color = METRIC_COLORS[metric]
        left.plot(
            rows["threshold"],
            rows["precision"],
            color=color,
            linewidth=LINE_WIDTH,
            marker="o",
            markersize=4,
            label=METRIC_LABELS[metric],
        )
        left.fill_between(
            rows["threshold"],
            rows["precision_ci_low"],
            rows["precision_ci_high"],
            color=color,
            alpha=0.12,
            linewidth=0,
        )
        right.plot(
            rows["threshold"],
            rows["coverage"],
            color=color,
            linewidth=LINE_WIDTH,
            marker="o",
            markersize=4,
            label=METRIC_LABELS[metric],
        )
    left.set_ylabel("Précision (1er résultat correct)")
    right.set_ylabel("Couverture (requêtes au-dessus du seuil)")
    for ax in (left, right):
        ax.set_xlabel("Seuil sur le score du 1er résultat")
        ax.set_ylim(0, 1.03)
        style_axes(ax, grid="both")
    set_title(left, "Fiabilité selon le seuil")
    set_title(right, "Requêtes conservées selon le seuil")
    style_legend(left, loc="lower right")
    set_suptitle(
        figure,
        f"Calibration des scores — {SETTING_LABELS[setting].lower()}, requêtes avec concurrents",
    )
    return figure


def analog_figure(queries: pd.DataFrame, metrics: Sequence[str], threshold: float) -> Figure:
    """Part des requêtes dont le premier analogue atteint au moins une similarité donnée."""
    figure = new_figure(8.5, 5.4)
    ax = figure.subplots()
    grid = np.linspace(0, 1, 201)
    reference = queries.drop_duplicates("query_index")

    def survival(values: np.ndarray) -> np.ndarray:
        values = np.nan_to_num(values, nan=-1.0)
        return np.array([(values >= x).mean() for x in grid])

    ax.plot(
        grid,
        survival(reference["oracle_tanimoto"].to_numpy()),
        color=COLOR_ORACLE,
        linewidth=LINE_WIDTH,
        label="Oracle (meilleure structure disponible)",
    )
    for metric in metrics:
        values = queries.loc[queries["metric"] == metric, "top1_tanimoto"].to_numpy()
        ax.plot(
            grid,
            survival(values),
            color=METRIC_COLORS[metric],
            linewidth=LINE_WIDTH,
            label=f"1er résultat — {METRIC_LABELS[metric]}",
        )
    ax.plot(
        grid,
        survival(reference["random_tanimoto"].to_numpy()),
        color=COLOR_RANDOM,
        linewidth=LINE_WIDTH,
        label="Tanimoto moyen des candidats (hasard)",
    )
    ax.axvline(threshold, color=INK_SECONDARY, linewidth=0.9)
    ax.annotate(
        f"analogue proche (≥ {threshold:g})",
        xy=(threshold, 1.0),
        xycoords=("data", "axes fraction"),
        xytext=(4, -12),
        textcoords="offset points",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Similarité structurale de Tanimoto avec la molécule recherchée")
    ax.set_ylabel("Part des requêtes ≥ cette similarité")
    style_axes(ax, grid="both")
    figure.legend(
        *ax.get_legend_handles_labels(),
        loc="outside lower center",
        ncols=2,
        frameon=False,
        fontsize=8.5,
        labelcolor=INK_SECONDARY,
    )
    set_title(ax, "Recherche d'analogues : molécule absente de la bibliothèque")
    return figure


def relation_figure(bins: pd.DataFrame, metrics: Sequence[str], baseline: float) -> Figure:
    """Similarité structurale moyenne selon la similarité spectrale (paires score > 0)."""
    figure = new_figure(11, 4.3)
    left, right = figure.subplots(1, 2, sharex=True)
    for metric in metrics:
        rows = bins[(bins["metric"] == metric) & (bins["n_pairs"] > 0)]
        centres = (rows["score_low"] + rows["score_high"]) / 2
        color = METRIC_COLORS[metric]
        left.plot(
            centres,
            rows["mean_tanimoto"],
            color=color,
            linewidth=LINE_WIDTH,
            marker="o",
            markersize=5,
            markeredgecolor=SURFACE,
            label=METRIC_LABELS[metric],
        )
        left.fill_between(
            centres, rows["ci_low"], rows["ci_high"], color=color, alpha=0.15, linewidth=0
        )
        right.plot(
            centres,
            rows["n_pairs"],
            color=color,
            linewidth=LINE_WIDTH,
            marker="o",
            markersize=5,
            markeredgecolor=SURFACE,
            label=METRIC_LABELS[metric],
        )
    left.axhline(baseline, color=COLOR_RANDOM, linewidth=LINE_WIDTH)
    left.annotate(
        "moyenne de tous les candidats",
        xy=(0.3, baseline),
        xytext=(0, -14),
        textcoords="offset points",
        fontsize=8.5,
        color=INK_SECONDARY,
    )
    left.set_ylabel("Tanimoto moyen (IC 95 %)")
    right.set_ylabel("Nombre de paires (échelle log)")
    right.set_yscale("log")
    for ax in (left, right):
        ax.set_xlabel("Score spectral (classes de 0,1)")
        ax.set_xlim(0, 1)
        style_axes(ax, grid="both")
    set_title(left, "Structure selon le score spectral")
    set_title(right, "Rareté des scores élevés")
    style_legend(left, loc="upper left")
    set_suptitle(figure, "Un spectre proche annonce-t-il une structure proche ?")
    return figure
