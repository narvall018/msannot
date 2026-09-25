"""Le dashboard s'exécute sans erreur (streamlit.testing) sur les données de démonstration."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

pytestmark = pytest.mark.integration
APP = Path(__file__).resolve().parents[2] / "app" / "streamlit_app.py"


def test_dashboard_example_search():
    app = AppTest.from_file(str(APP), default_timeout=240)
    app.run()
    assert not app.exception, app.exception
    assert len(app.tabs) == 4
    assert len(app.dataframe) >= 1  # table des résultats


def test_dashboard_analog_mode_and_pasted_peaks():
    app = AppTest.from_file(str(APP), default_timeout=240)
    app.run()
    app.sidebar.radio[1].set_value("analog").run()
    assert not app.exception, app.exception
    app.sidebar.radio[0].set_value("Pics collés").run()
    assert not app.exception, app.exception
