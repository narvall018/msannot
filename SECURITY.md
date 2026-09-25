# Politique de sécurité

## Versions prises en charge

| Version | Correctifs |
|---|---|
| 0.1.x | oui |

## Signaler une vulnérabilité

Merci de **ne pas ouvrir d'issue publique**. Utilisez le signalement privé de GitHub
(onglet **Security**, puis **Report a vulnerability**). Une première réponse est donnée sous
7 jours.

## Périmètre

- **Fichiers d'entrée.** Les fichiers MSP/MGF et les SMILES sont analysés par msannot et
  RDKit. N'analysez pas de fichiers d'origine inconnue sur une machine sensible.
- **Réseau.** Seule la commande `fetch` accède au réseau, en HTTPS vers GitHub (versions de
  MassBank).
- **Dashboard.** Streamlit n'a pas d'authentification et accepte des fichiers importés. Il
  est prévu pour un usage local : ne l'exposez pas sur Internet sans protection.
- **Docker.** L'image s'exécute avec un utilisateur non root.
