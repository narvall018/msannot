# Démarrage rapide

```bash
make install && source .venv/bin/activate
msannot demo                        # environ 1 minute
xdg-open results/demo/report.html   # macOS : open
```

## 1. Rechercher un spectre

Trois spectres d'exemple (caféine, carbamazépine, venlafaxine) ont été retirés de la
bibliothèque. D'autres spectres des mêmes molécules, acquis par d'autres laboratoires, y
restent : la recherche doit les retrouver.

```bash
msannot search data/demo/example_queries.mgf --library data/demo/massbank_open.msp.gz
```

```text
        MSBNK-ACES_SU-AS000913 — Caffeine (m/z 195.0873)
┃ rank ┃ name     ┃ score ┃ n_matches ┃ delta_mz ┃ contributor ┃
│ 1    │ Caffeine │ 0.951 │ 7         │ +0.0004  │ LCSB        │
│ 2    │ Caffeine │ 0.949 │ 7         │ +0.0004  │ Eawag       │
...
```

Les résultats sont écrits dans `results/search/hits.tsv`, avec un spectre miroir par requête
dans `results/search/figures/`.

Recherche d'analogues, c'est-à-dire de molécules voisines quand le composé n'est pas dans la
bibliothèque :

```bash
msannot search mes_spectres.mgf -l data/demo/massbank_open.msp.gz --mode analog --metric modified_cosine
```

## 2. Évaluer les méthodes

`msannot demo` évalue sur un échantillon de requêtes :

- l'**identification**, avec la molécule présente dans la bibliothèque ;
- la **recherche d'analogues**, avec la molécule retirée ;
- la **relation spectre/structure** ;
- le **réseau moléculaire**.

```
results/demo/
├── report.html                    # rapport autonome : données, paramètres, résultats, limites
├── benchmark_summary.json         # toutes les métriques, paramètres, versions
├── identification_*.tsv           # par requête, résumé, calibration, tests de McNemar
├── analog_*.tsv                   # par requête, résumé, tests de Wilcoxon
├── relation_bins.tsv, network_edges.tsv
└── figures/*.png                  # 10 figures, dont 3 spectres miroirs
```

L'évaluation complète, sur toutes les requêtes, prend environ 3 minutes :
`msannot benchmark config/full.yaml` (ou `make full`).

## 3. Explorer

```bash
make run    # dashboard : recherche interactive, bibliothèque, résultats du benchmark
```
