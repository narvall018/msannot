# Dépannage

### « Python >= 3.11 requis »

Le `python3` par défaut est trop ancien, par exemple Anaconda en 3.9. Indiquez un
interpréteur explicitement :

```bash
make install PYTHON=/usr/bin/python3.12
```

### `ImportError: libXrender.so.1` (RDKit, dessin des molécules)

Sur une image Linux minimale, le rendu Cairo de RDKit a besoin de bibliothèques X11 :

```bash
sudo apt-get install libxrender1 libxext6
```

Le Dockerfile les installe déjà.

### `msannot demo` : « config/demo.yaml introuvable »

Lancez la commande depuis le dépôt (ou un de ses sous-dossiers), ou utilisez
`msannot benchmark <config.yaml>`.

### « trop peu de pics après nettoyage »

Un spectre doit garder au moins 5 pics, une fois retirés ceux proches du précurseur et ceux
sous 1 % du pic de base. Les spectres très pauvres, souvent de petites molécules, ne peuvent
pas être comparés de façon fiable.

### Aucun résultat en mode `identity`

Aucun spectre de la bibliothèque n'a un précurseur à moins de `--ppm` de celui de la requête.
Vérifiez l'adduit (les bibliothèques préparées ne contiennent que des [M+H]+) et l'étalonnage
de masse, ou élargissez `--ppm`.

### `msannot fetch` échoue

Le téléchargement (environ 670 Mo) passe par GitHub. En cas de coupure, relancez : les
fichiers complets déjà présents ne sont pas retéléchargés, et un fichier partiel
(`.part`) n'est jamais laissé comme fichier final.

### Tests de parité ignorés (`skipped`)

matchms ou ms_entropy ne sont pas installés : `pip install -e ".[ref]"`.

### Docker : `permission denied … docker.sock`

L'utilisateur n'appartient pas au groupe `docker`. Utilisez `sudo docker …`, ou
`sudo usermod -aG docker $USER` puis reconnectez-vous. Ce second choix donne des droits
équivalents à root ; réservez-le à une machine personnelle. La CI construit et teste l'image
à chaque push.

### Le dashboard est lent au premier affichage

Le chargement et l'indexation de 23 000 spectres prennent environ 5 s, puis sont mis en cache.
