# Journal des modifications

Format inspiré de [Keep a Changelog](https://keepachangelog.com/fr/1.1.0/) ; versions selon
[Semantic Versioning](https://semver.org/lang/fr/).

## [Non publié]

## [0.1.1] — 2026-09-25

### Ajouté

- `msannot search --min-score` : ne garde que les résultats au-dessus d'un seuil de score.
- `msannot search --no-figures` : résultats TSV seulement.
- `msannot search` avertit quand l'adduit d'une requête est absent de la bibliothèque.
- `msannot validate` affiche le nombre de laboratoires contributeurs.
- `make check` : lint, types et tests en une commande.
- Marqueur `py.typed` (PEP 561), `.editorconfig`, `.gitattributes`.
- Dependabot pour les GitHub Actions.
- Tests : recherche de `config/demo.yaml`, bilan des doublons exacts.

### Corrigé

- Image Docker : ajout de `libexpat1`, requise par le rendu des molécules de RDKit. La démo
  échouait dans le conteneur ; le problème a été détecté par la CI.

### Modifié

- CI : durée maximale par job.
- Dépannage : bibliothèques système de RDKit et vérification avec `ldd`.

## [0.1.0] — 2026-09-25

Première version publique.

### Ajouté

- **Données** :
  - téléchargement de MassBank (`fetch`) ;
  - licences lues en flux dans l'export JSON-LD ;
  - curation avec RDKit (`prepare`) et bilan des exclusions ;
  - bibliothèque de démonstration de 23 491 spectres CC BY 4.0 / CC0, avec l'attribution de
    chaque spectre.
- **Similarités** :
  - cosine et modified cosine gloutons, identiques à matchms ;
  - entropie spectrale pondérée (Li et al., 2021), décomposée en somme sur les paires ;
  - moteur de recherche vectorisé sur les fragments et les pertes neutres.
- **Recherche** (`search`) : identification (fenêtre en ppm) et analogues (écart de masse),
  spectres miroirs avec structures RDKit.
- **Évaluation** (`benchmark`, `demo`) :
  - identification tous laboratoires et autres laboratoires, requêtes avec isomères
    concurrents ;
  - rang pessimiste, doublons exacts exclus ;
  - calibration des seuils de score, tests de McNemar ;
  - analogues face au hasard et à l'oracle, tests de Wilcoxon ;
  - relation spectre/structure, réseau moléculaire type GNPS.
- **Sorties** : TSV, JSON, figures, rapport HTML autonome (`report`).
- **Dashboard** Streamlit.
- **Qualité** :
  - 77 tests, dont la parité avec matchms (0.31, 0.33) et ms_entropy ;
  - ruff, mypy, pre-commit ;
  - CI GitHub Actions (Python 3.11–3.13, démo, Docker).
- **Documentation** : README, 9 documents dans `docs/`, script de validation des
  similarités.

[Non publié]: https://github.com/narvall018/msannot/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/narvall018/msannot/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/narvall018/msannot/releases/tag/v0.1.0
