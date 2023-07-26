~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
git branch : gbZ/Vy_init_nulle/dev
    commit : a54a6ecfe3644d444d248595537ed71b4e59ded4
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
La redaction de la fiche n'est pas terminee mais elle permet deja de valider separement deux fonctioninalites de fort interet pour le code.

I - Ce que la fiche valide
----------------------------
    
1) Superposition vitesse reprise et vitesse connue
 Lors d'une simulation en reprise il y a desormais la possibilite de supperposer le champs de vitesse de reprise a un champ dont on connait l'expression.
 Pour appliquer cette superposition, il faut mettre le mot clef " ajout_init_a_reprise " (si possible le placer soit proche du mot "nom_reprise", soit proche des mots "expression_v[i]_init"). Le chmap donne par "expression_v[i]_init" correspond au champ qui sera superpose au champ de reprise.

2) Repetabilite du tirage aleatoire
 Lors d'une simulation avec le forcage spectral utilise de sorte a superposer un ecoulement de THI a l'ecoulement simule (bloc forcage, type 3) un tirage de nombre aleatoire intervient.
 En ajoutant le mot clef "random_fixed" la suite des tirages aleatoires dans la simulation devient peterministe (ou pseudo-deterministe). C'est a dire que l'iteration N de toute simulation genere systematiquement le meme champ de THI des lors que le mot clef "random_fixed" est present.



II - Ce qu'il reste a faire 
-----------------------------

1) Rediger ces remarques dans le prm, et produire les post-traitement qui permettent de les justifier tres clairement.

2) Developpemet IJK pour realiser une reprise du champ de force spectral. Dans l'utilisation de type 3, on doit aussi faire une reprise du generateur de nombres aleatoire. Un MWE pour la reprise de nombre aleatoires est disponible ici : gr262753@is243146:/volatile/Premier_Projet_CPP/TEST_RANDOM_GENERATOR/version_avec_mt19937 . 
