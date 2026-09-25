# Reproductibilité

## Ce qui est garanti

| Élément | Mécanisme |
|---|---|
| Données identiques | bibliothèque versionnée, `SHA256SUMS`, gzip déterministe, version MassBank fixée (2026.03) |
| Provenance | SHA-256 de l'export d'origine dans le bilan de curation ; attribution de chaque spectre |
| Calcul identique | graines fixes pour les échantillonnages, tris stables, rang pessimiste explicite |
| Implémentations validées | parité avec matchms et ms_entropy (`scripts/validate_similarity.py`, tests) |
| Paramètres traçables | configuration complète dans `benchmark_summary.json` et dans le rapport |
| Environnement traçable | versions de Python, RDKit, numpy, scipy, pandas, matplotlib et networkx dans chaque rapport |
| Environnement reconstructible | `requirements.txt` (versions exactes), `Dockerfile` |
| Figures du README | copiées d'une exécution réelle (`make figures`) |
| Vérification continue | GitHub Actions : lint, typage, tests sous Python 3.11, 3.12 et 3.13, démo, image Docker |

## Reproduire les résultats documentés

```bash
make install-locked && source .venv/bin/activate
cd data/demo && sha256sum -c SHA256SUMS && cd ../..
msannot benchmark config/full.yaml     # docs/results.md (environ 3 min)
python scripts/validate_similarity.py  # tableau de validation de docs/methodology.md
```

Sans Python :

```bash
docker build -t msannot . && docker run --rm -v "$PWD/results/docker:/app/results" msannot
```

## Ce qui peut varier

- **Versions des dépendances.** Sous Python 3.13, les contraintes de matchms font installer
  RDKit 2024.09 et numpy 2.2 au lieu de 2026.03 et 2.5. Les tests passent avec les deux
  jeux de versions.
- **Temps de calcul.** Ils dépendent de la machine et sont donnés à titre indicatif.
- **MassBank en ligne.** Une nouvelle version de MassBank modifie la bibliothèque. Les
  résultats documentés portent sur la version 2026.03.

## Environnement de référence

```
Python 3.12.3 — Linux 6.8 (glibc 2.39), x86_64
msannot 0.1.0 · rdkit 2026.03.6 · numpy 2.5.3 · scipy 1.16.3 · pandas 3.0.6
matplotlib 3.11.2 · networkx 3.4.2 · pydantic 2.13.5
matchms 0.33.1 · ms_entropy 1.5.2 (tests de parité)
```

Les tests passent aussi sous Python 3.11.16 (RDKit 2026.03, matchms 0.33) et 3.13.15
(RDKit 2024.09, matchms 0.31). La démo donne des résultats identiques avec les versions de
`requirements.txt` (scipy 1.18, networkx 3.7).
