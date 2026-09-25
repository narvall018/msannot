"""Réseau moléculaire spectral et sa cohérence structurale.

Construction inspirée de GNPS (Wang et al., 2016) : un nœud par composé (un spectre
représentatif), une arête si le modified cosine dépasse ``min_score`` avec au moins
``min_matches`` pics appariés, et si chaque nœud figure parmi les ``top_k`` meilleurs
voisins de l'autre.

Évaluation : les arêtes relient-elles des structures plus proches que des paires prises
au hasard ? (Tanimoto des arêtes contre Tanimoto de paires aléatoires ; test de
Mann-Whitney.)
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from msannot.chem.structures import FingerprintIndex
from msannot.config import NetworkConfig, SimilarityConfig
from msannot.evaluation.stats import Proportion
from msannot.search.library import SpectralLibrary
from msannot.similarity.index import LibraryIndex


@dataclass(frozen=True, slots=True)
class NetworkResult:
    """Graphe, table des arêtes et bilan."""

    graph: Any
    edges: pd.DataFrame
    summary: dict[str, Any]


def _top_neighbours(
    index: LibraryIndex, spectra: Sequence[Any], config: NetworkConfig
) -> list[dict[int, tuple[float, int]]]:
    neighbours: list[dict[int, tuple[float, int]]] = []
    for position, spectrum in enumerate(spectra):
        result = index.query(spectrum, config.metric)
        keep = (
            (result.indices != position)
            & (result.scores >= config.min_score)
            & (result.n_matches >= config.min_matches)
        )
        indices, scores, matches = result.indices[keep], result.scores[keep], result.n_matches[keep]
        order = np.lexsort((indices, -matches, -scores))[: config.top_k]
        neighbours.append({int(indices[i]): (float(scores[i]), int(matches[i])) for i in order})
    return neighbours


def build_network(
    library: SpectralLibrary,
    nodes: np.ndarray,
    fingerprints: FingerprintIndex,
    similarity_config: SimilarityConfig,
    config: NetworkConfig,
    threshold: float,
    seed: int = 0,
) -> NetworkResult:
    """Construit le réseau sur les spectres ``nodes`` et mesure sa cohérence structurale."""
    spectra = [library.spectra[int(i)] for i in nodes]
    index = LibraryIndex(spectra, similarity_config)
    neighbours = _top_neighbours(index, spectra, config)

    graph = nx.Graph()
    for position, spectrum in enumerate(spectra):
        graph.add_node(
            position,
            identifier=spectrum.identifier,
            name=spectrum.name,
            inchikey14=spectrum.inchikey14,
        )
    rows = []
    for a, candidates in enumerate(neighbours):
        for b, (score, matches) in candidates.items():
            if a < b and a in neighbours[b]:
                tanimoto = fingerprints.similarity(spectra[a].inchikey14, spectra[b].inchikey14)
                graph.add_edge(a, b, score=score, matches=matches, tanimoto=tanimoto)
                rows.append(
                    {
                        "node_a": spectra[a].identifier,
                        "node_b": spectra[b].identifier,
                        "name_a": spectra[a].name,
                        "name_b": spectra[b].name,
                        "score": score,
                        "matches": matches,
                        "tanimoto": tanimoto,
                    }
                )
    edges = pd.DataFrame(
        rows, columns=["node_a", "node_b", "name_a", "name_b", "score", "matches", "tanimoto"]
    )

    rng = np.random.default_rng(seed)
    n_random = max(len(edges), 1000)
    first = rng.integers(0, len(spectra), n_random)
    second = rng.integers(0, len(spectra), n_random)
    distinct = first != second
    random_tanimoto = np.array(
        [
            fingerprints.similarity(spectra[a].inchikey14, spectra[b].inchikey14)
            for a, b in zip(first[distinct], second[distinct], strict=True)
        ]
    )
    components = [c for c in nx.connected_components(graph) if len(c) > 1]
    edge_tanimoto = edges["tanimoto"].to_numpy()
    test = (
        stats.mannwhitneyu(edge_tanimoto, random_tanimoto, alternative="greater")
        if edge_tanimoto.size
        else None
    )
    summary = {
        "n_nodes": graph.number_of_nodes(),
        "n_edges": graph.number_of_edges(),
        "n_connected_nodes": int(sum(len(c) for c in components)),
        "n_components": len(components),
        "largest_component": max((len(c) for c in components), default=0),
        "edge_tanimoto_median": float(np.median(edge_tanimoto)) if edge_tanimoto.size else None,
        "edge_close": Proportion(
            int((edge_tanimoto >= threshold).sum()), int(edge_tanimoto.size)
        ).to_dict(),
        "random_tanimoto_median": float(np.median(random_tanimoto)),
        "random_close": Proportion(
            int((random_tanimoto >= threshold).sum()), int(random_tanimoto.size)
        ).to_dict(),
        "mann_whitney_p": float(test.pvalue) if test is not None else None,
        "random_tanimoto": random_tanimoto,
    }
    return NetworkResult(graph=graph, edges=edges, summary=summary)
