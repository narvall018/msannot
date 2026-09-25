# Utilisation

## Ligne de commande

```text
msannot [--verbose | --quiet] [--log-file FICHIER] COMMANDE [OPTIONS]
```

| Commande | Rôle |
|---|---|
| `fetch` | télécharge une version de MassBank (export MSP + export JSON des licences) |
| `prepare` | cure un export MassBank : filtres, contrôles de structure, nettoyage |
| `validate FICHIER` | vérifie un fichier MSP/MGF et résume son contenu |
| `search REQUÊTES -l BIBLIOTHÈQUE` | recherche des spectres ; écrit les résultats et des spectres miroirs |
| `benchmark CONFIG.yaml` | évaluation complète et rapport HTML |
| `report DOSSIER` | régénère le rapport sans recalcul |
| `demo` | benchmark de démonstration (`config/demo.yaml`) |

Les erreurs attendues (fichier illisible, paramètre invalide…) produisent un message en
français et un code de sortie 1, sans trace Python.

### Construire sa propre bibliothèque MassBank

```bash
msannot fetch --release 2026.03 --outdir data/raw             # environ 670 Mo
msannot prepare --msp data/raw/2026.03/MassBank_NISTformat.msp \
                --licenses data/raw/2026.03/MassBank.json \
                -o data/processed/massbank_all.msp.gz          # toutes licences
msannot prepare ... --open-only                                # CC BY 4.0 / CC0 seulement
msannot prepare ... --all-instruments                          # garde aussi la basse résolution
```

Seules les bibliothèques sous licence ouverte peuvent être redistribuées ; les autres restent
dans `data/processed/` (ignoré par Git).

### Rechercher

```bash
msannot search requetes.mgf -l bibliotheque.msp.gz \
    --mode identity --metric entropy --ppm 10 --top 5 -o results/search
msannot search requetes.mgf -l bibliotheque.msp.gz \
    --mode analog --metric modified_cosine --max-shift 200
```

- `identity` : candidats de même masse (± ppm).
- `analog` : candidats dont la masse diffère d'au plus `--max-shift` Da.
- Mesures disponibles : `cosine`, `modified_cosine`, `entropy`.
- `--min-score 0.85` : ne garde que les résultats au-dessus d'un seuil. Dans l'évaluation
  complète (autres laboratoires, requêtes avec isomères), 0,85 en entropie ou 0,90 en cosine
  donnent environ 95 % de premiers résultats corrects (`docs/results.md`).
- `--no-figures` : résultats TSV seulement, sans spectre miroir (utile pour de nombreuses
  requêtes).
- Un avertissement s'affiche si l'adduit d'une requête (par exemple [M+Na]+) est absent de
  la bibliothèque : la fenêtre de précurseur ne peut alors pas retrouver le bon composé.

Formats de requête acceptés : MGF (`BEGIN IONS` / `PEPMASS=`) ou MSP, compressés en gzip ou
non. Les requêtes sont nettoyées comme la bibliothèque.

## Configuration du benchmark

[`config/default.yaml`](../config/default.yaml) documente **toutes** les options et leur
valeur par défaut :

```yaml
dataset: ../data/demo/massbank_open.msp.gz
output_dir: ../results/mon_benchmark
similarity:
  tolerance: 0.01
  metrics: [cosine, modified_cosine, entropy]
identification:
  precursor_ppm: 10
  max_queries: null        # null = toutes
analogs:
  max_mass_shift: 200
  close_analog_threshold: 0.7
network:
  enabled: true
```

Les chemins relatifs sont résolus par rapport au fichier de configuration. Une clé inconnue
ou une valeur invalide est refusée avant tout calcul.

**Bibliothèque requise.** Un fichier MSP dont chaque spectre porte une structure (`SMILES`)
et un identifiant de composé (`InChIKey`). Le laboratoire (`Contributor`) est utilisé pour
le scénario inter-laboratoires. Les bibliothèques produites par `msannot prepare` ont toutes
ces métadonnées.

## Sorties du benchmark

| Fichier | Contenu |
|---|---|
| `identification_queries.tsv` | par requête, mesure et scénario : candidats, concurrents, rang (pessimiste), scores |
| `identification_summary.tsv` | exactitude top-1/3/10 avec IC de Wilson, référence aléatoire |
| `identification_calibration.tsv` | précision et couverture en fonction du seuil de score |
| `identification_mcnemar.tsv` | comparaisons appariées des mesures |
| `analog_queries.tsv` | par requête et mesure : 1er résultat, Tanimoto, oracle, hasard |
| `analog_summary.tsv`, `analog_tests.tsv` | médianes, IC bootstrap, tests de Wilcoxon |
| `relation_bins.tsv` | Tanimoto moyen par classe de score |
| `network_edges.tsv` | arêtes du réseau (score, pics appariés, Tanimoto) |
| `benchmark_summary.json` | tout ce qui précède, avec paramètres, versions et SHA-256 des données |
| `report.html` | rapport autonome |

## Dashboard

```bash
make run     # ou : streamlit run app/streamlit_app.py
```

- **Recherche** : exemples, fichier importé ou pics collés ; tableau des résultats ;
  similarité de Tanimoto quand la structure de la requête est connue ; spectre miroir ;
  structure proposée.
- **Bibliothèque** : composition et filtrage par nom.
- **Benchmark** : résultats de `msannot demo`.
- **Méthode.**

## Bibliothèque Python

```python
from pathlib import Path
from msannot.config import CleaningConfig, SimilarityConfig
from msannot.io import load_msp, load_spectra
from msannot.processing import clean_all, clean_spectrum
from msannot.search import SpectralLibrary

spectra, _ = clean_all(load_msp(Path("data/demo/massbank_open.msp.gz")), CleaningConfig())
library = SpectralLibrary(spectra, SimilarityConfig())
query = clean_spectrum(load_spectra(Path("data/demo/example_queries.mgf"))[0], CleaningConfig())
print(library.search(query, metric="entropy", mode="identity", top_k=3))
```
