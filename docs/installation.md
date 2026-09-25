# Installation

## Prérequis

| Outil | Version | Remarque |
|---|---|---|
| Python | **3.11, 3.12 ou 3.13** | testé sur les trois en CI |
| Git | récent | cloner le dépôt |
| make | optionnel | raccourcis `make install`, `make demo`… |
| Docker | optionnel | exécution sans installer Python |

Tout s'installe avec `pip`, **RDKit compris** : les roues PyPI de `rdkit` suffisent, conda
n'est pas nécessaire.

> **Attention aux Python trop anciens.** Un `python3` fourni par une distribution Anaconda
> ancienne peut être en 3.9 ou 3.10. Le Makefile cherche `python3.12`, puis `python3.13`,
> `python3.11` et enfin `python3`, et s'arrête avec un message clair si la version ne convient
> pas.

## Installation recommandée

```bash
git clone https://github.com/narvall018/msannot.git
cd msannot
make install              # .venv + msannot + dashboard + outils de dev + matchms/ms_entropy
source .venv/bin/activate
msannot --version
```

## Installation manuelle

```bash
python3.12 -m venv .venv && source .venv/bin/activate
pip install -e ".[app,dev,ref]"     # ou "pip install -e ." pour la CLI seule
```

| Extra | Contenu |
|---|---|
| *(aucun)* | CLI, recherche, évaluation, figures, rapport |
| `app` | dashboard Streamlit |
| `dev` | pytest, ruff, mypy, pre-commit, stubs |
| `ref` | matchms et ms_entropy : implémentations de référence, pour les tests de parité seulement |

## Versions exactes

```bash
make install-locked      # requirements.txt : versions testées (Python 3.12, Linux)
```

## Docker

```bash
docker build -t msannot .
docker run --rm -v "$PWD/results/docker:/app/results" msannot     # démonstration
```

## Vérification

```bash
msannot validate data/demo/massbank_open.msp.gz
make test          # 77 tests, moins d'une minute
```
