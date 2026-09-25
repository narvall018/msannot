"""Dashboard Streamlit de msannot.

Couche d'interface uniquement : chargement, paramètres, affichage. La recherche, les
similarités et les figures viennent du package ``msannot`` (aucune logique dupliquée).

Lancement : ``streamlit run app/streamlit_app.py`` (ou ``make run``).
"""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import streamlit as st

from msannot import __version__
from msannot.chem.structures import FingerprintIndex, molecule_svg, parse_structure
from msannot.config import METRIC_LABELS, CleaningConfig, SimilarityConfig
from msannot.exceptions import MsannotError
from msannot.io.mgf import load_spectra
from msannot.io.msp import load_msp
from msannot.models import Spectrum
from msannot.plotting import dataset_figure, mirror_figure
from msannot.processing.cleaning import clean_all, clean_spectrum
from msannot.search.library import SpectralLibrary

ROOT = Path(__file__).resolve().parents[1]
DEMO_LIBRARY = ROOT / "data" / "demo" / "massbank_open.msp.gz"
EXAMPLE_QUERIES = ROOT / "data" / "demo" / "example_queries.mgf"
DEMO_RESULTS = ROOT / "results" / "demo"
METRICS = list(METRIC_LABELS)

st.set_page_config(page_title="msannot", page_icon="🧪", layout="wide")


# ---------------------------------------------------------------------------- cache
@st.cache_resource(show_spinner="Chargement et indexation de la bibliothèque…", max_entries=3)
def load_library(path: str, modified: float) -> SpectralLibrary:
    """Bibliothèque nettoyée et indexée (clé : chemin + date de modification)."""
    spectra, _ = clean_all(load_msp(Path(path)), CleaningConfig())
    return SpectralLibrary(spectra, SimilarityConfig())


def parse_pasted(text: str, precursor: float) -> Spectrum | None:
    """Pics collés sous la forme « m/z intensité », un pic par ligne."""
    peaks = []
    for line in text.strip().splitlines():
        parts = line.replace(",", " ").replace(";", " ").split()
        if len(parts) >= 2:
            try:
                peaks.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    if not peaks or precursor <= 0:
        return None
    array = np.array(peaks)
    return Spectrum(
        "requête collée", array[:, 0], array[:, 1], precursor, {"name": "requête collée"}
    )


def uploaded_path(uploaded: Any) -> Path:
    suffix = "".join(Path(uploaded.name).suffixes) or ".mgf"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as handle:
        handle.write(uploaded.getvalue())
    return Path(handle.name)


# ---------------------------------------------------------------------------- barre latérale
st.sidebar.title("🧪 msannot")
st.sidebar.caption(f"v{__version__} — annotation de spectres MS/MS")
library_file = st.sidebar.file_uploader(
    "Bibliothèque MSP (défaut : démonstration)", type=["msp", "gz"]
)
library_path = uploaded_path(library_file) if library_file else DEMO_LIBRARY
try:
    library = load_library(str(library_path), library_path.stat().st_mtime)
except MsannotError as exc:
    st.error(f"Bibliothèque illisible : {exc}")
    st.stop()

st.sidebar.subheader("Requête")
source = st.sidebar.radio("Source", ["Exemples", "Fichier MGF/MSP", "Pics collés"])
queries: list[Spectrum] = []
if source == "Exemples":
    queries = load_spectra(EXAMPLE_QUERIES)
elif source == "Fichier MGF/MSP":
    uploaded = st.sidebar.file_uploader("Spectres requêtes", type=["mgf", "msp", "gz"])
    if uploaded is not None:
        try:
            queries = load_spectra(uploaded_path(uploaded))
        except MsannotError as exc:
            st.sidebar.error(str(exc))
else:
    precursor = st.sidebar.number_input(
        "m/z du précurseur", min_value=0.0, value=195.0877, format="%.4f"
    )
    text = st.sidebar.text_area(
        "Pics (m/z intensité)", "110.0712 30\n138.0662 100\n123.0430 20\n83.0604 10\n69.0447 5"
    )
    pasted = parse_pasted(text, precursor)
    queries = [pasted] if pasted else []

query = None
if queries:
    labels = [f"{q.name} ({q.identifier})" for q in queries]
    query = queries[
        st.sidebar.selectbox("Spectre", range(len(queries)), format_func=lambda i: labels[i])
    ]

st.sidebar.subheader("Paramètres")
mode = st.sidebar.radio(
    "Mode",
    ["identity", "analog"],
    horizontal=True,
    format_func=lambda m: "Identification" if m == "identity" else "Analogues",
)
metric = st.sidebar.selectbox(
    "Mesure", METRICS, index=2 if mode == "identity" else 1, format_func=lambda m: METRIC_LABELS[m]
)
top_k = st.sidebar.slider("Résultats", 1, 25, 10)
ppm = st.sidebar.number_input("Fenêtre de précurseur (ppm)", 1.0, 100.0, 10.0)
max_shift = st.sidebar.number_input("Écart de masse maximal (Da, analogues)", 1.0, 1000.0, 200.0)
exclude_same = st.sidebar.checkbox(
    "Exclure le composé de la requête (simule une inconnue)", value=mode == "analog"
)

# ---------------------------------------------------------------------------- onglets
st.title("Annotation de petites molécules par spectres MS/MS")
tabs = st.tabs(["Recherche", "Bibliothèque", "Benchmark", "Méthode"])

with tabs[0]:
    if query is None:
        st.info(
            "Choisissez un exemple, importez un fichier ou collez des pics dans la barre latérale."
        )
    else:
        cleaned = clean_spectrum(query, CleaningConfig())
        if cleaned is None:
            st.warning("Moins de 5 pics après nettoyage : spectre trop pauvre pour une recherche.")
        else:
            exclude = library.table["identifier"].to_numpy() == query.identifier
            if exclude_same and query.inchikey14:
                exclude |= library.inchikeys == query.inchikey14
            hits = library.search(
                cleaned,
                metric=metric,
                mode=mode,
                top_k=top_k,
                ppm=ppm,
                max_mass_shift=max_shift,
                exclude=exclude,
            )
            if hits.empty:
                st.warning("Aucun spectre de la bibliothèque ne partage de pic avec la requête.")
            else:
                if query.smiles and parse_structure(query.smiles):
                    structures = {
                        "__query__": query.smiles,
                        **dict(zip(hits["identifier"], hits["smiles"], strict=True)),
                    }
                    fingerprints = FingerprintIndex(structures)
                    hits["tanimoto"] = [
                        fingerprints.similarity("__query__", identifier)
                        if identifier in fingerprints
                        else np.nan
                        for identifier in hits["identifier"]
                    ]
                st.dataframe(hits.drop(columns=["library_index"]), hide_index=True)
                st.download_button(
                    "Télécharger les résultats (TSV)",
                    hits.to_csv(sep="\t", index=False),
                    "hits.tsv",
                    "text/tab-separated-values",
                )
                rank = st.selectbox("Spectre miroir du résultat n°", list(hits["rank"]))
                row = hits[hits["rank"] == rank].iloc[0]
                hit = library.spectra[int(row["library_index"])]
                tanimoto = (
                    float(row["tanimoto"])
                    if "tanimoto" in hits and pd.notna(row.get("tanimoto"))
                    else None
                )
                st.pyplot(
                    mirror_figure(
                        cleaned,
                        hit,
                        metric,
                        SimilarityConfig(),
                        title=f"{query.name} ↔ {hit.name}",
                        tanimoto=tanimoto,
                    )
                )
                svg = molecule_svg(hit.smiles, 320, 220)
                if svg:
                    st.caption(f"Structure proposée : {hit.name} ({hit.get('formula')})")
                    st.image(svg, width=320)

with tabs[1]:
    table = library.table
    st.metric("Spectres", f"{len(table):,}".replace(",", " "))
    st.pyplot(dataset_figure(table))
    search_text = st.text_input("Filtrer par nom de composé")
    shown = (
        table[table["name"].str.contains(search_text, case=False, na=False)]
        if search_text
        else table
    )
    st.dataframe(shown.head(500), hide_index=True)
    st.caption(f"{len(shown)} spectres correspondants (500 premiers affichés).")

with tabs[2]:
    summary_path = DEMO_RESULTS / "benchmark_summary.json"
    if not summary_path.is_file():
        st.info(
            "Aucun résultat de benchmark : lancez `msannot demo` (ou `make demo`), "
            "puis rechargez la page."
        )
    else:
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        st.caption(
            f"Résultats de `msannot demo` du {summary['created_utc'][:16].replace('T', ' ')} UTC"
        )
        identification = pd.DataFrame(summary["identification"]["summary"])
        st.subheader("Identification (top-1)")
        st.dataframe(
            identification[
                [
                    "setting",
                    "subset",
                    "metric",
                    "n_queries",
                    "top1",
                    "top1_ci_low",
                    "top1_ci_high",
                    "random_top1",
                ]
            ].round(3),
            hide_index=True,
        )
        st.subheader("Analogues")
        st.dataframe(pd.DataFrame(summary["analogs"]["summary"]).round(3), hide_index=True)
        for key in ("identification", "calibration", "analogs", "relation", "network"):
            image = DEMO_RESULTS / summary["figures"].get(key, "")
            if key in summary["figures"] and image.is_file():
                st.image(str(image))

with tabs[3]:
    st.markdown(
        """
        - **Nettoyage** : pics au-delà de précurseur − 1,6 Da et sous 1 % du pic de base retirés.
        - **Cosine / modified cosine** : appariement glouton (tolérance 0,01 Da), conventions de
          matchms ; le modified cosine apparie aussi les pics décalés de la différence de masse.
        - **Entropie spectrale** pondérée (Li et al., 2021).
        - **Identification** : candidats de même masse (± ppm) ; **analogues** : composé exclu,
          candidats à ± 200 Da ; similarité structurale par Tanimoto (Morgan, RDKit).

        Une annotation par bibliothèque reste une hypothèse (niveau 2 de la Metabolomics
        Standards Initiative) : elle doit être confirmée par un standard. Détails :
        `docs/methodology.md`.
        """
    )
