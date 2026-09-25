"""Nettoyage des spectres avant comparaison.

Étapes (paramètres dans :class:`~msannot.config.CleaningConfig`) :

1. retrait des pics de m/z supérieur à ``précurseur − marge`` (1,6 Da par défaut, comme
   Li et al., 2021) : le précurseur résiduel et ses isotopes ne sont pas des fragments, et
   ils gonfleraient artificiellement le modified cosine ;
2. retrait des pics d'intensité relative inférieure au seuil (1 % du pic de base) ;
3. conservation des ``max_peaks`` pics les plus intenses ;
4. rejet des spectres trop pauvres (moins de ``min_peaks`` pics) ;
5. normalisation des intensités (pic de base = 1).
"""

from __future__ import annotations

import numpy as np

from msannot.config import CleaningConfig
from msannot.models import Spectrum


def clean_spectrum(spectrum: Spectrum, config: CleaningConfig) -> Spectrum | None:
    """Renvoie une copie nettoyée du spectre, ou ``None`` s'il reste trop peu de pics."""
    mz, intensities = spectrum.mz, spectrum.intensities
    keep = (mz <= spectrum.precursor_mz - config.precursor_margin) & (intensities > 0)
    mz, intensities = mz[keep], intensities[keep]
    if mz.size == 0:
        return None
    keep = intensities >= config.min_relative_intensity * intensities.max()
    mz, intensities = mz[keep], intensities[keep]
    if config.max_peaks is not None and mz.size > config.max_peaks:
        # plus intenses d'abord ; à intensité égale, le plus petit m/z (déterministe)
        order = np.lexsort((mz, -intensities))[: config.max_peaks]
        order.sort()
        mz, intensities = mz[order], intensities[order]
    if mz.size < config.min_peaks:
        return None
    return spectrum.with_peaks(mz, intensities / intensities.max())


def clean_all(spectra: list[Spectrum], config: CleaningConfig) -> tuple[list[Spectrum], int]:
    """Nettoie une liste de spectres ; renvoie les spectres conservés et le nombre rejeté."""
    kept = [cleaned for spectrum in spectra if (cleaned := clean_spectrum(spectrum, config))]
    return kept, len(spectra) - len(kept)
