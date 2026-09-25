# Méthodologie

Ce document décrit ce que fait msannot, pourquoi, et avec quelles hypothèses. Les chiffres
viennent de `msannot benchmark config/full.yaml` et de `scripts/validate_similarity.py`
(msannot 0.1.0, MassBank 2026.03, exécution du 25/09/2026). Ils sont reproductibles avec les
fichiers du dépôt.

## 1. Question

> Avec quelle fiabilité peut-on reconnaître une petite molécule à partir de son spectre de
> fragmentation (MS/MS) en le comparant à une bibliothèque ? Et quand la molécule est
> absente de la bibliothèque, les spectres les plus proches désignent-ils des structures
> proches ?

En métabolomique non ciblée ou en analyse environnementale, la plupart des signaux mesurés
ne sont pas identifiés. La comparaison à une bibliothèque de spectres de référence est la
première étape d'annotation ; la recherche d'analogues sert quand la molécule n'y figure pas.

## 2. Données et curation (`processing/curation.py`)

Source : **MassBank Europe**, version 2026.03 : export NIST MSP, et export JSON-LD qui
contient la licence de chaque spectre (Horai et al., 2010). Les filtres ci-dessous sont
appliqués dans cet ordre et chacun est compté :

| Filtre | Raison | Spectres exclus |
|---|---|---|
| spectre MS2 | les MS1 et MSⁿ ne se comparent pas aux MS2 | 22 028 |
| mode positif | un seul mode, pour des spectres comparables | 38 390 |
| adduit [M+H]+ | le précurseur doit correspondre à la même espèce ionique | 9 605 |
| instrument haute résolution (QTOF, Orbitrap/FT, TOF) | la tolérance de 0,01 Da suppose une haute résolution | 5 168 |
| licence CC BY 4.0 ou CC0 | redistribution dans le dépôt | 29 969 |
| précurseur lisible | — | 287 |
| SMILES valide (RDKit) | la vérité de référence est la structure | 260 |
| InChIKey déclaré = InChIKey recalculé | cohérence de l'annotation | 0 |
| molécule neutre | [M+H]+ n'a pas de sens pour un cation permanent | 92 |
| précurseur à ±10 ppm (ou ±0,005 Da) de la masse théorique [M+H]+ | détecte une erreur de structure ou d'adduit | 39 |
| au moins 5 pics après nettoyage | spectre trop pauvre pour être comparé | 9 674 |

Sur 139 006 spectres, **23 494** sont conservés, soit 3 958 composés. Trois d'entre eux sont
mis de côté comme requêtes d'exemple, ce qui laisse une bibliothèque de **23 491 spectres**
issus de 24 contributeurs. L'erreur médiane sur la masse du précurseur est de 0,09 ppm.

**Composé.** Un composé est identifié par le premier bloc (14 caractères) de l'InChIKey
(Heller et al., 2015), qui encode la connectivité sans la stéréochimie. Des stéréo-isomères,
indiscernables en MS/MS classique, comptent donc comme un seul composé.

**Doublons.** 360 groupes de spectres strictement identiques (813 spectres) sont présents,
presque tous chez un même contributeur. Dans 55 groupes, les spectres identiques sont annotés
comme des composés différents. Un doublon exact n'est pas une mesure indépendante : il est
exclu des candidats de sa propre requête.

## 3. Nettoyage des spectres (`processing/cleaning.py`)

1. Retrait des pics de m/z supérieur à précurseur − 1,6 Da, comme Li et al. (2021). Le
   précurseur résiduel et ses isotopes ne sont pas des fragments, et gonfleraient
   artificiellement le modified cosine.
2. Retrait des pics inférieurs à 1 % du pic de base (bruit).
3. Conservation des 100 pics les plus intenses au maximum.
4. Rejet des spectres de moins de 5 pics.
5. Normalisation : pic de base = 1.

## 4. Mesures de similarité (`similarity/`)

**Appariement.** Deux pics correspondent si `|m/z_a − (m/z_b + décalage)| ≤ 0,01 Da`. Chaque
pic n'est apparié qu'une fois. Les paires sont choisies de façon **gloutonne**, par
contribution décroissante, avec les mêmes conventions et le même ordre en cas d'égalité que
matchms (Huber et al., 2020).

| Mesure | Définition | Référence |
|---|---|---|
| **Cosine** | Σ(√Ia·√Ib) sur les paires, normalisé par la norme de tous les pics (intensités^0,5) | Stein & Scott, 1994 |
| **Modified cosine** | idem, mais un pic peut aussi s'apparier après décalage de la différence de masse des précurseurs ; si les précurseurs sont confondus, c'est un cosine | Watrous et al., 2012 |
| **Entropie spectrale** | 1 − (2·S_AB − S_A − S_B) / ln 4, avec S l'entropie de Shannon des intensités normalisées ; spectres pauvres (S < 3) repondérés par I^(0,25 + 0,25·S) | Li et al., 2021 |

**Décomposition de l'entropie.** Pour des spectres normalisés (somme = 1), la similarité
d'entropie se réécrit :

```
similarité = Σ_(paires appariées) [ (a+b)·ln(a+b) − a·ln a − b·ln b ] / ln 4
```

Les pics non appariés ne contribuent pas. C'est la propriété exploitée par Flash entropy
(Li & Fiehn, 2023), et c'est elle qui permet à msannot de traiter les trois mesures avec le
même moteur.

## 5. Moteur de recherche vectorisé (`similarity/index.py`)

Comparer une requête à 23 491 spectres un par un serait lent. msannot procède donc ainsi :

1. Tous les pics de la bibliothèque sont concaténés et triés par m/z. Pour chaque pic de la
   requête, `searchsorted` donne en une opération tous les pics à moins de la tolérance.
2. Pour le modified cosine, la condition « fragment décalé de ΔP » équivaut à l'égalité des
   **pertes neutres** (précurseur − fragment). Un second index, trié par perte neutre, trouve
   ces paires alors que le décalage varie d'un spectre à l'autre.
3. Les contributions sont additionnées par spectre (`bincount`). L'appariement glouton exact
   n'est rejoué que pour les spectres où un même pic participe à plusieurs paires.
4. Les m/z publiés étant arrondis à 4 décimales, des écarts **exactement** égaux à la
   tolérance sont fréquents. La recherche utilise donc une fenêtre à peine élargie, puis chaque
   paire est revérifiée avec les mêmes opérations flottantes que la référence.

Temps mesurés : 2 ms (cosine) à 7 ms (modified cosine) par requête pour 6 000 spectres. Le
benchmark complet, soit plus de 22 000 requêtes d'identification et 4 000 d'analogues pour
trois mesures, prend environ 3 minutes.

## 6. Validation des implémentations

Commande : `python scripts/validate_similarity.py`. Paires de score non nul, requêtes tirées
au hasard (graine fixe) :

| Comparaison | Paires | Écart maximal | Pics appariés différents |
|---|---|---|---|
| moteur vectorisé / fonctions de référence (3 mesures) | 35 452 | 3,3 × 10⁻¹⁶ | 0 |
| cosine / matchms `CosineGreedy` | 6 000 | 3,3 × 10⁻¹⁶ | 0 |
| modified cosine / matchms `ModifiedCosineGreedy` | 6 000 | 4,4 × 10⁻¹⁶ | 0 |
| entropie / ms_entropy (Li et al.) | 6 000 | écart médian 6,6 × 10⁻⁹ | — |

Pour l'entropie, 34 paires diffèrent de plus de 10⁻³. Toutes ont une explication :

- **31** : un spectre contient deux pics distants de moins de deux tolérances. ms_entropy
  suppose ces pics fusionnés au préalable, et apparie dans l'ordre des m/z ; msannot choisit
  l'appariement qui contribue le plus, comme pour le cosine.
- **3** : un écart de m/z exactement égal à la tolérance, que la conversion en float32 de
  ms_entropy fait basculer.

Hors de ces deux cas, l'écart maximal est de 1,6 × 10⁻⁷. Les tests `test_reference_parity.py`
vérifient cette parité en continu. Ils passent avec matchms 0.31 et 0.33.

## 7. Protocole d'identification (`evaluation/identification.py`)

- **Requêtes** : tous les spectres dont le composé possède au moins un autre spectre, soit
  22 316 requêtes.
- **Candidats** : spectres dont le précurseur est à moins de 10 ppm de celui de la requête.
  La requête et ses doublons exacts sont exclus.
- **Scénarios** :
  - *tous laboratoires* ;
  - *autres laboratoires uniquement* : les candidats du même contributeur sont exclus. C'est
    le cas réaliste, car un spectre inconnu n'a pas de jumeau acquis sur le même instrument
    avec les mêmes réglages.
- **Rang pessimiste** : en cas d'égalité de score, le bon composé est classé après les
  candidats incorrects ex æquo.
- **Référence aléatoire** : probabilité de désigner le bon composé en tirant un candidat de
  même masse au hasard. C'est ce que donne la masse seule.
- **Requêtes avec concurrents** : requêtes dont la fenêtre de masse contient au moins un
  autre composé (isomère, isobare). Pour les autres, la masse suffit. Seules les premières
  mesurent vraiment le pouvoir discriminant du spectre.
- **Statistiques** :
  - intervalles de Wilson à 95 % (Wilson, 1927) ;
  - test de McNemar exact entre mesures, évaluées sur les mêmes requêtes (McNemar, 1947) ;
  - calibration : précision et couverture du premier résultat en fonction d'un seuil de score.

## 8. Protocole de recherche d'analogues (`evaluation/analogs.py`)

- **Requêtes** : une requête par composé, soit 3 958. Le spectre retenu est celui dont
  l'accession est la plus petite : un choix arbitraire, mais indépendant du contenu.
- **Molécule inconnue simulée** : tous les spectres du composé sont retirés des candidats.
- **Candidats** : écart de précurseur d'au plus 200 Da, soit une modification chimique
  plausible.
- **Mesure** : similarité de Tanimoto entre la requête et le premier résultat, sur empreintes
  de Morgan de rayon 2, 2048 bits, équivalentes aux ECFP4 (Rogers & Hahn, 2010), calculées
  avec RDKit.
- **Références** :
  - *hasard* : Tanimoto moyen de tous les candidats de la requête ;
  - *oracle* : meilleure structure disponible parmi les candidats, c'est-à-dire le mieux
    qu'une méthode puisse faire avec cette bibliothèque.
- **Critère « analogue proche »** : Tanimoto ≥ 0,7. Une requête sans aucun pic commun avec la
  bibliothèque compte comme un échec.
- **Statistiques** : médianes avec intervalle bootstrap (2 000 rééchantillonnages), test
  des rangs signés de Wilcoxon entre mesures et contre le hasard.

## 9. Relation spectre/structure et réseau moléculaire

- **Relation** : pour toutes les paires (requête, candidat) de score non nul, soit 16 à 25
  millions de paires selon la mesure, on calcule le Tanimoto moyen par classe de score, et
  une corrélation de Spearman.
- **Réseau** (Wang et al., 2016) :
  - un nœud par composé ;
  - une arête si le modified cosine est ≥ 0,7 avec au moins 6 pics appariés, et si chaque
    nœud figure parmi les 10 meilleurs voisins de l'autre ;
  - cohérence évaluée par le Tanimoto des arêtes, comparé à autant de paires tirées au
    hasard (test de Mann-Whitney).

## 10. Choix des paramètres

Tous les paramètres ont été fixés **a priori**, d'après la littérature et les usages. Aucun
n'a été optimisé sur les résultats.

| Paramètre | Valeur | Justification |
|---|---|---|
| tolérance fragments | 0,01 Da | instruments haute résolution |
| intensités cosine | ^0,5 | usage courant (GNPS), atténue la domination d'un pic |
| entropie | pondérée | recommandation de Li et al. (2021) |
| fenêtre de précurseur | 10 ppm | haute résolution, avec marge pour les écarts d'étalonnage |
| écart de masse (analogues) | 200 Da | au-delà, l'interprétation « même squelette modifié » devient peu plausible |
| réseau | 0,7 · 6 pics · top 10 | paramètres par défaut de GNPS |
| Tanimoto « proche » | 0,7 | seuil usuel pour des analogues sur ECFP4 |

Les requêtes de la démo sont échantillonnées (5 000 et 1 000, graine 0) pour qu'elle reste
rapide. `config/full.yaml` utilise toutes les requêtes.

## 11. Limites

- **Une seule bibliothèque.** Les résultats décrivent MassBank, restreint au sous-ensemble
  sous licence ouverte : laboratoires, classes chimiques (pesticides, médicaments,
  contaminants), énergies de collision. Ils ne se transposent pas tels quels à une autre
  bibliothèque ni à des spectres expérimentaux bruts.
- **Répétitions intra-laboratoire.** Plusieurs spectres d'un même composé viennent souvent du
  même laboratoire, à différentes énergies de collision. D'où le scénario « autres
  laboratoires ».
- **Annotations imparfaites.** Des spectres identiques peuvent porter deux annotations
  différentes (55 groupes). L'évaluation mesure l'accord avec les annotations, pas la vérité
  chimique.
- **Une seule mesure structurale.** Le Tanimoto sur empreintes de Morgan est une similarité
  parmi d'autres ; un seuil de 0,7 reste conventionnel.
- **Pas de données orthogonales.** Temps de rétention, mobilité ionique, isotopes ne sont pas
  utilisés : une annotation par bibliothèque correspond au niveau 2 de la Metabolomics
  Standards Initiative (Sumner et al., 2007) et doit être confirmée par un standard.
- **Pas de méthodes apprises.** Spec2Vec et MS2DeepScore ne sont pas comparés dans cette
  version (voir la roadmap).

## Références

- Heller S.R. et al. (2015). InChI, the IUPAC International Chemical Identifier. *Journal of
  Cheminformatics* 7:23.
- Horai H. et al. (2010). MassBank: a public repository for sharing mass spectral data for
  life sciences. *Journal of Mass Spectrometry* 45(7):703–714.
- Huber F. et al. (2020). matchms – processing and similarity evaluation of mass spectrometry
  data. *Journal of Open Source Software* 5(52):2411.
- Huber F. et al. (2021). Spec2Vec: Improved mass spectral similarity scoring through learning
  of structural relationships. *PLoS Computational Biology* 17(2):e1008724.
- Li Y. et al. (2021). Spectral entropy outperforms MS/MS dot product similarity for
  small-molecule compound identification. *Nature Methods* 18:1524–1531.
- Li Y., Fiehn O. (2023). Flash entropy search to query all mass spectral libraries in real
  time. *Nature Methods* 20:1475–1478.
- McNemar Q. (1947). Note on the sampling error of the difference between correlated
  proportions or percentages. *Psychometrika* 12(2):153–157.
- Rogers D., Hahn M. (2010). Extended-connectivity fingerprints. *Journal of Chemical
  Information and Modeling* 50(5):742–754.
- Stein S.E., Scott D.R. (1994). Optimization and testing of mass spectral library search
  algorithms for compound identification. *Journal of the American Society for Mass
  Spectrometry* 5(9):859–866.
- Sumner L.W. et al. (2007). Proposed minimum reporting standards for chemical analysis.
  *Metabolomics* 3:211–221.
- Wang M. et al. (2016). Sharing and community curation of mass spectrometry data with Global
  Natural Products Social Molecular Networking. *Nature Biotechnology* 34(8):828–837.
- Watrous J. et al. (2012). Mass spectral molecular networking of living microbial colonies.
  *PNAS* 109(26):E1743–E1752.
- Wilson E.B. (1927). Probable inference, the law of succession, and statistical inference.
  *Journal of the American Statistical Association* 22(158):209–212.
