"""Figures matplotlib. Ce module ne fait aucun calcul métier : il reçoit des résultats."""

from msannot.plotting.dataset import curation_figure, dataset_figure
from msannot.plotting.evaluation import (
    analog_figure,
    calibration_figure,
    identification_figure,
    relation_figure,
)
from msannot.plotting.network import network_figure
from msannot.plotting.spectra import mirror_figure
from msannot.plotting.style import save_figure

__all__ = [
    "analog_figure",
    "calibration_figure",
    "curation_figure",
    "dataset_figure",
    "identification_figure",
    "mirror_figure",
    "network_figure",
    "relation_figure",
    "save_figure",
]
