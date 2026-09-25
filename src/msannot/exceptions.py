"""Exceptions propres à msannot.

Toutes dérivent de :class:`MsannotError`, ce qui permet à la CLI de distinguer les erreurs
attendues (entrée invalide, configuration incorrecte) des bogues.
"""


class MsannotError(Exception):
    """Erreur de base du package."""


class InputFormatError(MsannotError):
    """Fichier d'entrée illisible, mal formé ou non pris en charge."""


class ConfigError(MsannotError):
    """Configuration invalide ou incohérente."""


class AnalysisError(MsannotError):
    """Analyse impossible avec les données et paramètres fournis."""


class DataFetchError(MsannotError):
    """Échec du téléchargement de données publiques."""
