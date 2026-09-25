# Données

## Bibliothèque de démonstration (`data/demo/`)

| Fichier | Contenu |
|---|---|
| `massbank_open.msp.gz` | 23 491 spectres MS/MS curés (4,9 Mo) |
| `massbank_open.attribution.tsv.gz` | pour chaque spectre : accession, nom, contributeur, licence, URL de la fiche MassBank |
| `massbank_open.curation.json` | bilan de curation : exclusions par raison, paramètres, SHA-256 de l'export d'origine, versions |
| `example_queries.mgf` | 3 spectres requêtes (caféine, carbamazépine, venlafaxine) retirés de la bibliothèque |
| `SHA256SUMS` | sommes de contrôle (`cd data/demo && sha256sum -c SHA256SUMS`) |

**Source.** MassBank Europe, version **2026.03** (publiée le 15/04/2026, dépôt GitHub
`MassBank/MassBank-data`), téléchargée le 25/09/2026. SHA-256 de l'export MSP :
`cf49e1b2…` (valeur complète dans `massbank_open.curation.json`).

**Licences.** Chaque spectre de MassBank a sa propre licence. La bibliothèque de
démonstration ne contient que des spectres sous **CC BY 4.0** (19 577) ou **CC0 1.0** (3 914),
qui autorisent la redistribution. CC BY impose de citer la source : l'attribution de chaque
spectre (contributeur, licence, lien vers la fiche d'origine) figure dans
`massbank_open.attribution.tsv.gz`. Les spectres sous licence non commerciale (CC BY-NC,
CC BY-NC-SA), à partage à l'identique (CC BY-SA) ou sous licence allemande DL-DE BY ne sont
pas redistribués. Ils restent utilisables localement via `msannot fetch` puis
`msannot prepare`.

**Contenu.** 3 958 composés, 24 contributeurs. Les plus représentés sont Eawag, mFam, LCSB,
NaToxAq, Athens_Univ, HBM4EU et UFZ. On y trouve surtout des médicaments, des pesticides, des
contaminants environnementaux et des produits naturels. Tous les spectres sont en mode
positif, adduit [M+H]+, sur des instruments haute résolution (Orbitrap/FT 70 %, QTOF/TOF 30 %).

**Reconstruction.**

```bash
make data     # = msannot fetch --release 2026.03 + scripts/build_demo_data.py
```

La bibliothèque et le fichier d'attribution sont reproduits octet pour octet (gzip sans
horodatage). Le bilan JSON, lui, contient la date d'exécution.

## Qualité des données : points connus

- **Doublons exacts.** 360 groupes de spectres identiques ; 55 sont annotés comme des
  composés différents. Ils sont exclus des candidats de leur propre requête lors de
  l'évaluation.
- **Contributeurs liés.** Certains contributeurs partagent des données : NaToxAq est un projet
  hébergé à l'UFZ, et 5 groupes de doublons exacts relient les deux.
- **Énergies de collision.** Un même composé est souvent mesuré à plusieurs énergies par le
  même laboratoire, d'où le scénario d'évaluation « autres laboratoires ».

## Utiliser ses propres données

| Entrée | Format | Exigences |
|---|---|---|
| Requêtes | MGF ou MSP, gzip accepté | `PEPMASS` / `PrecursorMZ` ; au moins 5 pics après nettoyage |
| Bibliothèque pour la recherche | MSP | précurseur ; structure (`SMILES`) et `InChIKey` recommandés |
| Bibliothèque pour le benchmark | MSP | `SMILES` + `InChIKey` obligatoires ; `Contributor` pour le scénario inter-laboratoires |

Les dossiers `data/raw/` (téléchargements) et `data/processed/` (bibliothèques non
redistribuables) sont ignorés par Git.
