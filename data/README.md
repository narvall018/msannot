# Données

| Dossier | Contenu | Versionné |
|---|---|---|
| `demo/` | bibliothèque MassBank curée (CC BY 4.0 / CC0), attribution, bilan de curation, requêtes d'exemple, `SHA256SUMS` | oui |
| `raw/` | exports MassBank téléchargés par `msannot fetch` (environ 670 Mo) | non |
| `processed/` | bibliothèques produites localement, par exemple avec toutes les licences | non |

## Bibliothèque de démonstration

- **Source** : MassBank Europe, version 2026.03 (`MassBank/MassBank-data` sur GitHub).
- **Contenu** : 23 491 spectres MS/MS, [M+H]+, mode positif, haute résolution ; 3 958
  composés ; 24 contributeurs.
- **Licences** : uniquement **CC BY 4.0** (19 577 spectres) et **CC0 1.0** (3 914 spectres).
  L'attribution de chaque spectre (accession, contributeur, licence, lien vers la fiche
  MassBank) est dans `demo/massbank_open.attribution.tsv.gz`.
- **Vérifier** : `cd data/demo && sha256sum -c SHA256SUMS`.
- **Reconstruire** : `make data` (téléchargement de MassBank + `scripts/build_demo_data.py`).

Provenance détaillée, étapes de curation et problèmes de qualité connus :
[../docs/data.md](../docs/data.md).
