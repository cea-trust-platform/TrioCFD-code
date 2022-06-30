#!/bin/bash

#############################################################################
# Appelle le python-visit pour faire la suite d'images pouis appelle fait_film
# pour le rendre en mp4
# mode d'emploi : ./fait_film_total.sh fichier_lata    nom_champ first_last_state chemin_logos liste_logos
#                                         1              2           3          4
# 1 : chemin jusqu'au lata, en mettant .lata
# 2 : Nom du champ tel qu'il apparait dans visit - dnas les latas
# 3 : A donner comme ca : [first_state,last_state]
# 4 : chemin jusqu'au dossier des logos a apposer
# 5 : A donner comme ca : [logo1,logo2,logo3] ou [logo1,logo2,0] ou [loo1,0,0]
#############################################################################


chemin_ijk=/volatile/IJK/env_IJK.sh
chemin_lata=$1
nom_champ=$2
nom_fichier_movie=movie_$nom_champ
first_last_state=$3

chemin_images=./MOVIE/
chemin_logos=$4
liste_logos=$5
nom_film=$nom_fichier_movie

source $chemin_ijk
mkdir -p $chemin_images

visit -cli -nowin -s $project_directory/share/PyTools/SOME_PYTHON/vs_makemovie.py $chemin_lata $nom_champ $chemin_images $nom_fichier_movie $first_last_state
#visit -s $project_directory/share/PyTools/SOME_PYTHON/vs_makemovie.py $chemin_lata $nom_champ $chemin_images $nom_fichier_movie $first_last_state

./fait_film $chemin_images $chemin_logos $liste_logos $nom_film.mp4
