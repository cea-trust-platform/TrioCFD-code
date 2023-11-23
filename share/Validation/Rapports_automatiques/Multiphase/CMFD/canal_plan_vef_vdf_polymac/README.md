# Fiche validation canal plan

Cette fiche permet lancer des simulation d'un écoulement turbulent dans un canal plan dans les configurations suivantes :
"VEF_k-epsilon", "VDF_k-epsilon", "VDF_k-tau", "VDF_k-omega", "PolyMAC_k-tau", "PolyMAC_k-omega"

2 jdd de base sont mis en place (un pour Pb_hydraulique_turbulent et l'autre pour Pb_Multiphase)

un script functions.py contient plusieurs fonctions utiles au pré et post traitement des calculs

### Pour ajouter une configuration il faut :
  - Aller dans src/ et ouvrir le jdd correspondant à la configuration à ajouter
  - Voir s'il y a un mot clé à modifier en plus de ceux qui sont déjà modifié par le code (les mots qui commence par $...)
  - Si oui ajoutez une nouvelle variable ($var) à la place du mot clé
  - Aller dans functions.py
  - mettre à jour les listes : available_config, available_config_pb_hydr et available_config_pb_multi 
  - modifier si nécessaire le dictionnaire param_config (qui définit comme écrire les jdd pour chaque config)
  - Aller dans la fiche
  - Ajouter la nouvelle configuration dans config
  
### Remarque : 
toujours appeler les configurations de la façon suivante : "{method}_{model}"
où method peut etre "VEF", "VDF" ou "PolyMAC" et model peut etre "k-epsilon", "k-omega" ou "k-tau"

