# Commandes courantes de msannot. `make help` liste les cibles.
# Python >= 3.11 requis : on privilégie python3.12, puis 3.13, 3.11, puis python3.

PYTHON ?= $(shell command -v python3.12 || command -v python3.13 || command -v python3.11 || command -v python3)
VENV   ?= .venv
BIN    := $(VENV)/bin
PORT   ?= 8501
IMAGE  ?= msannot:latest
RELEASE ?= 2026.03

.DEFAULT_GOAL := help
.PHONY: help venv install install-locked demo full test test-fast coverage lint format typecheck check \
        pre-commit run report figures data validate-similarity docker-build docker-demo docker-app clean

help: ## Affiche cette aide
	@grep -E '^[a-zA-Z_-]+:.*?## ' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-19s\033[0m %s\n", $$1, $$2}'

$(BIN)/python:
	@$(PYTHON) -c 'import sys; sys.exit(0 if sys.version_info >= (3, 11) else "Python >= 3.11 requis (trouvé : " + sys.version.split()[0] + ")")'
	$(PYTHON) -m venv $(VENV)
	$(BIN)/python -m pip install --quiet --upgrade pip

venv: $(BIN)/python ## Crée l'environnement virtuel (.venv)

install: venv ## Installe msannot + dashboard + outils de dev + bibliothèques de référence
	$(BIN)/pip install --quiet -e ".[app,dev,ref]"
	@echo "OK : activez l'environnement avec 'source $(VENV)/bin/activate'"

install-locked: venv ## Installe les versions exactes testées (requirements.txt)
	$(BIN)/pip install --quiet -r requirements.txt
	$(BIN)/pip install --quiet --no-deps -e .
	$(BIN)/pip install --quiet -e ".[dev,ref]"

demo: ## Démonstration (≈ 1 min) : résultats et rapport dans results/demo
	$(BIN)/msannot demo

full: ## Évaluation complète, toutes les requêtes (≈ 4 min) : results/full
	$(BIN)/msannot benchmark config/full.yaml

test: ## Lance tous les tests
	$(BIN)/pytest

test-fast: ## Tests unitaires uniquement
	$(BIN)/pytest -m "not integration"

coverage: ## Tests avec couverture
	$(BIN)/pytest --cov=msannot --cov-report=term-missing

lint: ## Style (ruff) et types (mypy)
	$(BIN)/ruff check src tests app scripts
	$(BIN)/ruff format --check src tests app scripts
	$(BIN)/mypy

check: lint test ## Lint, types et tests : à lancer avant un commit ou une PR

format: ## Reformate le code
	$(BIN)/ruff format src tests app scripts
	$(BIN)/ruff check --fix src tests app scripts

typecheck: ## Types uniquement
	$(BIN)/mypy

pre-commit: ## Exécute tous les hooks pre-commit
	PATH="$(abspath $(BIN)):$$PATH" $(BIN)/pre-commit run --all-files

run: ## Dashboard Streamlit (http://localhost:8501)
	$(BIN)/streamlit run app/streamlit_app.py --server.port $(PORT)

report: ## Régénère le rapport de la démo sans recalcul
	$(BIN)/msannot report results/demo

figures: ## Met à jour les figures du README depuis results/full
	$(BIN)/python scripts/export_readme_figures.py results/full

validate-similarity: ## Rapport de parité avec matchms et ms_entropy
	$(BIN)/python scripts/validate_similarity.py

data: ## Retélécharge MassBank et reconstruit data/demo (≈ 700 Mo téléchargés)
	$(BIN)/msannot fetch --release $(RELEASE) --outdir data/raw
	$(BIN)/python scripts/build_demo_data.py data/raw/$(RELEASE)

docker-build: ## Construit l'image Docker
	docker build -t $(IMAGE) .

docker-demo: ## Démonstration dans Docker (résultats dans results/docker)
	mkdir -p results/docker
	docker run --rm -v "$(CURDIR)/results/docker:/app/results" $(IMAGE)

docker-app: ## Dashboard dans Docker (http://localhost:8501)
	docker run --rm -p $(PORT):8501 --entrypoint streamlit $(IMAGE) \
		run app/streamlit_app.py --server.address 0.0.0.0 --server.headless true

clean: ## Supprime caches, builds et résultats générés
	rm -rf build dist .pytest_cache .mypy_cache .ruff_cache .coverage coverage.xml htmlcov
	find . -path ./$(VENV) -prune -o -name "__pycache__" -type d -print -exec rm -rf {} +
	find results -mindepth 1 ! -name .gitkeep -exec rm -rf {} +
