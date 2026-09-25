# Contribuer à msannot

## Environnement

```bash
git clone https://github.com/narvall018/msannot.git && cd msannot
make install && source .venv/bin/activate   # inclut matchms et ms_entropy (tests de parité)
pre-commit install
make test && make lint
```

## Déroulement

1. Ouvrez une *issue* avant un changement important.
2. Travaillez sur une branche (`git switch -c feat/…`), avec de petits commits explicites.
3. Vérifiez : `make lint`, `make test`, et `msannot demo` si les résultats peuvent changer.
4. Ouvrez une *pull request* en remplissant le modèle. La CI doit être verte.

## Conventions

| Sujet | Règle |
|---|---|
| Langue | documentation, docstrings et messages en français ; identifiants en anglais |
| Style | `ruff format` (100 colonnes), `ruff check` sans erreur |
| Types | annotations complètes dans `src/`, `mypy` sans erreur |
| Erreurs | sous-classes de `MsannotError` pour les erreurs attendues |
| Journal | `get_logger(__name__)`, jamais de `print` dans le package (sauf la CLI) |
| Architecture | aucune logique scientifique dans `cli.py` ni dans `app/` ; les figures ne calculent rien |

## Tests

- Toute correction de bug s'accompagne d'un test qui échouait avant.
- Une nouvelle similarité doit être testée contre sa version de référence paire par paire
  (`test_similarity.py`), et contre une implémentation externe si elle existe.
- Pas de réseau dans les tests : simulez les téléchargements (`monkeypatch`).
- Si les résultats de l'évaluation changent, mettez à jour `docs/results.md`, les figures
  (`make full && make figures`) et le `CHANGELOG.md`.

## Ajouter une mesure de similarité

Voir [docs/architecture.md](docs/architecture.md#ajouter-une-mesure-de-similarité) : la
transformation d'intensité et la contribution d'une paire suffisent pour une mesure
décomposable. L'évaluation, les figures, le rapport et le dashboard l'intègrent sans autre
modification.

## Règles scientifiques

- Aucun chiffre écrit à la main : tout résultat documenté vient d'une exécution du code.
- Pas de réglage de paramètres sur les données d'évaluation sans le signaler.
- Distinguer ce qui est mesuré de ce qui est interprété.
