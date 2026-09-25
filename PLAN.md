# PLAN — msannot

Feuille de route du projet : objectifs, choix de conception et étapes. Document vivant.

## 1. Question scientifique

> Avec quelle fiabilité peut-on identifier une petite molécule à partir de son spectre de
> fragmentation (MS/MS) en le comparant à une bibliothèque ? Et quand la molécule est
> absente de la bibliothèque, les spectres les plus proches désignent-ils des structures
> chimiquement proches ?

Trois niveaux d'évaluation :

1. **Identification** : la molécule est présente dans la bibliothèque (autres spectres du
   même composé). On mesure le rang du bon composé, dans deux scénarios :
   - *tous laboratoires* : le spectre correct peut venir du même laboratoire ;
   - *autres laboratoires* : seuls les spectres d'autres laboratoires sont autorisés, ce
     qui est plus réaliste et plus difficile.
2. **Recherche d'analogues** : la molécule est retirée de la bibliothèque. On mesure la
   similarité structurale (Tanimoto, empreintes de Morgan) entre la molécule et le premier
   résultat, comparée au hasard et à la meilleure structure disponible (oracle).
3. **Relation spectre/structure** : la similarité spectrale prédit-elle la similarité
   structurale ? Est-ce vrai aussi pour le réseau moléculaire (type GNPS) ?

## 2. Données

- **MassBank Europe**, version 2026.03 (export NIST MSP + métadonnées JSON-LD avec la licence
  de chaque spectre).
- **Curation** :
  - MS2, mode positif, adduit [M+H]+, instruments haute résolution ;
  - SMILES valide (RDKit), InChIKey cohérent, masse du précurseur cohérente avec la
    structure (10 ppm) ;
  - au moins 5 pics après nettoyage.
- **Jeu de démonstration versionné** : uniquement les spectres sous **CC BY 4.0** ou
  **CC0**, avec l'attribution de chaque spectre (accession, laboratoire, licence, URL).
- **Jeu complet** (toutes licences) : reproductible par `msannot fetch` + `msannot prepare`,
  mais non redistribué.

## 3. Méthodes

| Élément | Choix |
|---|---|
| Nettoyage | pics > précurseur − 1,6 Da retirés, intensité relative < 1 %, 100 pics max., ≥ 5 pics |
| Cosine | appariement glouton (tolérance 0,01 Da), intensités^0,5, conventions de matchms |
| Modified cosine | appariement des fragments **ou** des pertes neutres (précurseur − fragment) |
| Entropie | similarité d'entropie spectrale pondérée (Li et al., 2021), décomposée en somme sur les paires appariées (Li & Fiehn, 2023) |
| Moteur de recherche | index vectorisé (fragments + pertes neutres) sur toute la bibliothèque, glouton exact seulement pour les paires en conflit |
| Structures | RDKit : InChIKey, masse exacte, empreintes de Morgan (r = 2, 2048 bits), Tanimoto |
| Statistiques | IC de Wilson, test de McNemar exact (métriques appariées), test des rangs signés de Wilcoxon, corrélation de Spearman, IC bootstrap |
| Ex æquo | rang pessimiste (le bon composé est classé après les ex æquo incorrects) |

**Validation des implémentations** :
- cosine et modified cosine comparés à **matchms** ;
- entropie comparée à **ms_entropy** (implémentation officielle de Li et al.) ;
- moteur vectorisé comparé aux fonctions de référence paire par paire.

## 4. Architecture

```
src/msannot/
├── io/            # MSP, MGF, téléchargement MassBank, licences, écriture
├── chem/          # RDKit : structures, masses, empreintes, dessins
├── processing/    # nettoyage des spectres, curation du jeu de données
├── similarity/    # cosine, modified cosine, entropie ; index vectorisé
├── search/        # bibliothèque et moteur de recherche (identification, analogues)
├── evaluation/    # identification, analogues, spectre/structure, réseau, statistiques
├── plotting/      # figures (spectres miroir, courbes, réseau)
├── report/        # rapport HTML autonome
├── pipeline.py    # préparation du jeu de données, benchmark
└── cli.py         # fetch | prepare | validate | search | benchmark | report | demo
app/               # dashboard Streamlit (interface uniquement)
```

## 5. Étapes

- [x] Choix du sujet ; exploration de MassBank (licences, volumes, qualité)
- [x] Vérification préalable : formule d'entropie décomposée == ms_entropy
- [ ] Squelette, configuration, E/S (MSP, MGF, licences)
- [ ] Chimie (RDKit) et curation
- [ ] Similarités de référence + index vectorisé + tests de parité matchms/ms_entropy
- [ ] Moteur de recherche
- [ ] Évaluations (identification, analogues, spectre/structure, réseau) + statistiques
- [ ] Figures, rapport, pipeline, CLI
- [ ] Jeu de démonstration (CC BY/CC0) + attribution
- [ ] Tests (unitaires, intégration, CLI, dashboard)
- [ ] Dashboard Streamlit
- [ ] Docker, Makefile, CI, pre-commit
- [ ] Documentation et README avec résultats réels
- [ ] Vérification finale, publication GitHub (nouveau dépôt public)

## 6. Règles scientifiques

- Aucun chiffre écrit à la main : tout vient d'une exécution réelle, avec version et date.
- Séparation stricte des niveaux technique, statistique et biologique/chimique.
- Paramètres par défaut fixés a priori (littérature) ; toute comparaison sur les données de
  démonstration est présentée comme telle.
