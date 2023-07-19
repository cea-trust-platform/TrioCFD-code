# Fiche validation conduite circulaire

Cette fiche permet lancer des simulation d'un écoulement turbulent dans une conduite circulaire dans les configurations suivantes :
"VEF_k-epsilon", "VDF_k-epsilon", "VDF_k-tau", "VDF_k-omega", "PolyMAC_k-tau", "PolyMAC_k-omega"

4 jdd de base sont mis en place 

un script functions.py contient plusieurs fonctions utiles au pré et post traitement des calculs

### Maillages

Avant de lancer la fiche assurez-vous que les fichiers de maillage ".med" sont présents dans src/

Pour créer les maillage :
- 1- charger les modules de salome dans un terminal
- 2- executer dans src/ : "salome maillage_tube_cut_vef.py"
- 3- executer dans src/ : "salome mesh_0_polymac.py"
- 4- executer dans src/ : "salome mesh_1_polymac.py"

### Pour ajouter une configuration il faut :
  - Aller dans src/ et ouvrir le jdd correspondant à la configuration à ajouter
  - Voir s'il y a un mot clé à modifier en plus de ceux qui sont déjà modifié par le code (les mots qui commence par $...)
  - Si oui ajoutez une nouvelle variable ($var) à la place du mot clé
  - Aller dans functions.py
  - mettre à jour les listes : available_config, available_config_pb_hydr et available_config_pb_multi 
  - mettre à jour la fonction file_name si nécessaire
  - modifier si nécessaire le dictionnaire param_config (qui définit comme écrire les jdd pour chaque config)
  - Aller dans la fiche
  - Ajouter la nouvelle configuration dans config
  
### Remarque : 
toujours appeler les configurations de la façon suivante : "{method}_{model}"
où method peut etre "VEF", "VDF" ou "PolyMAC" et model peut etre "k-epsilon", "k-omega" ou "k-tau"
