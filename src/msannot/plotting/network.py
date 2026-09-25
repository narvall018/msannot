"""Réseau moléculaire : dessin des plus grandes composantes et cohérence structurale."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.figure import Figure

from msannot.plotting.style import (
    BLUE,
    COLOR_RANDOM,
    INK_SECONDARY,
    SURFACE,
    new_figure,
    set_suptitle,
    set_title,
    style_axes,
    style_legend,
)

# Rampe séquentielle à une teinte (bleu clair → bleu foncé) pour la similarité des arêtes.
TANIMOTO_CMAP = LinearSegmentedColormap.from_list("tanimoto", ["#cde2fb", "#3987e5", "#0d366b"])


def network_figure(
    graph: Any, random_tanimoto: np.ndarray, max_nodes: int = 260, seed: int = 0
) -> Figure:
    """Plus grandes composantes du réseau (arêtes colorées par Tanimoto) et distributions."""
    figure = new_figure(12, 5.4)
    grid = figure.add_gridspec(1, 2, width_ratios=[1.5, 1.0])
    left = figure.add_subplot(grid[0, 0])
    right = figure.add_subplot(grid[0, 1])

    components = sorted(
        (c for c in nx.connected_components(graph) if len(c) > 1), key=len, reverse=True
    )
    selected: list[int] = []
    for component in components:
        if len(selected) + len(component) > max_nodes and selected:
            break
        selected.extend(component)
    subgraph = graph.subgraph(selected)
    if subgraph.number_of_nodes():
        layout = nx.spring_layout(subgraph, seed=seed, k=0.35, iterations=100)
        tanimoto = [data["tanimoto"] for _, _, data in subgraph.edges(data=True)]
        nx.draw_networkx_edges(
            subgraph,
            layout,
            ax=left,
            edge_color=tanimoto,
            edge_cmap=TANIMOTO_CMAP,
            edge_vmin=0,
            edge_vmax=1,
            width=1.6,
        )
        nx.draw_networkx_nodes(
            subgraph,
            layout,
            ax=left,
            node_size=14,
            node_color=INK_SECONDARY,
            edgecolors=SURFACE,
            linewidths=0.6,
        )
        colorbar = figure.colorbar(
            ScalarMappable(norm=Normalize(0, 1), cmap=TANIMOTO_CMAP), ax=left, shrink=0.6, pad=0.01
        )
        colorbar.set_label("Tanimoto de l'arête", color=INK_SECONDARY)
    left.set_axis_off()
    set_title(left, f"{len(selected)} nœuds des plus grandes composantes")

    edges = np.array([data["tanimoto"] for _, _, data in graph.edges(data=True)])
    bins = np.linspace(0, 1, 26).tolist()
    right.hist(
        random_tanimoto,
        bins=bins,
        density=True,
        histtype="step",
        color=COLOR_RANDOM,
        linewidth=2,
        label=f"Paires au hasard (n = {len(random_tanimoto)})",
    )
    if edges.size:
        right.hist(
            edges,
            bins=bins,
            density=True,
            histtype="step",
            color=BLUE,
            linewidth=2,
            label=f"Arêtes du réseau (n = {edges.size})",
        )
    right.set_xlabel("Similarité structurale de Tanimoto")
    right.set_ylabel("Densité")
    style_axes(right)
    style_legend(right, loc="upper right")
    set_title(right, "Les arêtes relient-elles des structures proches ?")
    set_suptitle(figure, "Réseau moléculaire (modified cosine, type GNPS)")
    return figure
