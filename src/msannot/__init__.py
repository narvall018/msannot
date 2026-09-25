"""msannot — annotation de petites molécules par spectres MS/MS.

Organisation en couches indépendantes :

- :mod:`msannot.io` : lecture/écriture MSP et MGF, données MassBank ;
- :mod:`msannot.chem` : structures chimiques (RDKit) ;
- :mod:`msannot.processing` : nettoyage des spectres, curation du jeu de données ;
- :mod:`msannot.similarity` : cosine, modified cosine, entropie spectrale, index vectorisé ;
- :mod:`msannot.search` : recherche en bibliothèque (identification, analogues) ;
- :mod:`msannot.evaluation` : évaluation statistique ;
- :mod:`msannot.plotting`, :mod:`msannot.report` : figures et rapport ;
- :mod:`msannot.pipeline` : orchestration utilisée par la CLI et le dashboard.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("msannot")
except PackageNotFoundError:  # pragma: no cover - package non installé
    __version__ = "0.0.0+unknown"

__all__ = ["__version__"]
