#!/bin/bash

########################################################
# mode d'emploi : .:place_sondes.sh chemin_du_data chemin_ou_poser_le_data_avec_sondes retirer/ajouter
#                                        1                          2                         3 
#                                   [axial,plan,splan,mode]
#                                     4    5    6     7
# 1 : nom du jeu de donnees, sans l'extension .data
# 4 : Nombre de plans normaux a Y contenant une sonde X
#     Et nombre de plans normaux a Z contenant une sonde X
# 5 : Nombre de plans normaux a X contenant des segments Y et des segments Z
# 6 : Nombre de plans normaux à Z contenant des segments Y et des segments Z
# 7 : . -> segments isolés
#     v -> faisceaux de trois segments, disposes en coin
#     + -> faisceaux de cinq segments, disposes en croix
########################################################

if [ "$project_directory" == "" ]; then 
	echo "IJK environment not set! You need to source IJK to launch this script"
	exit 1
fi

chemin_python=${project_directory}/share/PyTools/probes/genSondes_GR.py
jdd=$1
new_jdd=$2
ajout=$3
a=$(echo "$4"| tr -d '\[' | tr -d '\]' | cut -d"," -f1)
p=$(echo "$4"| tr -d '\[' | tr -d '\]' | cut -d"," -f2)
s=$(echo "$4"| tr -d '\[' | tr -d '\]' | cut -d"," -f3)
m=$(echo "$4"| tr -d '\[' | tr -d '\]' | cut -d"," -f4)
echo "$a"
echo "$p"
echo "$s"

#  = On a deja un bloc de sondes. Il nous plait pas. On en veut un nouveau
if [ "$ajout" == "remplace" ]; then
	num_ligne_Sondes=$(grep -n -e "  Sondes" $jdd.data | cut -d":" -f1)
	echo "$num_ligne_Fin, $num_ligne_Sondes"
	echo "$(($num_ligne_Sondes +1))"
	# Ecriture jusqu'au bloc sondes dans le nouveau jdd
	head -n$(($num_ligne_Sondes - 1)) $jdd.data > $new_jdd.data
	# Ecriture du bloc sondes dans un fichier tampon
	python $chemin_python $jdd.data -a$a -p$p -n$s -f$new_jdd -m$m
	# Ajout concatenation du fichier tampon au nouveau jdd
	cat $new_jdd >> $new_jdd.data
	# On termine tranquillement le nouveau jdd
	echo "}" >> $new_jdd.data
	echo "  Fin" >> $new_jdd.data
fi

# = On a pas de bloc de sondes. On est jaloux de ceux qu'en ont. On en fait un
if [ "$ajout" == "cree" ]; then
	num_ligne_Fin=$(grep -n -e "Fin" $jdd.data | cut -d":" -f1)
        # Ecriture jusqu'à presque la fin sondes dans le nouveau jdd
        head -n$(($num_ligne_Fin - 2)) $jdd.data > $new_jdd.data
        # Ecriture du bloc sondes dans un fichier tampon
        python $chemin_python $jdd.data -a$a -p$p -n$s -f$new_jdd -m$m
        # Ajout concatenation du fichier tampon au nouveau jdd
        cat $new_jdd >> $new_jdd.data
        # On termine tranquillement le nouveau jdd
        echo "}" >> $new_jdd.data
        echo "  Fin" >> $new_jdd.data
fi

# = On a un bloc de sondes. Il est pas mal. On veut plus. On ajoute des sondes a celles existantes
if [ "$ajout" == "ajout" ]; then
        num_ligne_Fin=$(grep -n -e "Fin" $jdd.data | cut -d":" -f1)
        # Ecriture jusqu'à presque la fin sondes dans le nouveau jdd
        head -n$(($num_ligne_Fin - 3)) $jdd.data > $new_jdd.data
	# Ecriture du bloc sondes dans un fichier tampon
        python $chemin_python $jdd.data -a$a -p$p -n$s -f$new_jdd -m$m
        sed -in '1,2'd $new_jdd
	sed -i '/}/'d $new_jdd
	# Ajout concatenation du fichier tampon au nouveau jdd
        cat $new_jdd >> $new_jdd.data
        # On termine tranquillement le nouveau jdd
        echo "}" >> $new_jdd.data
        echo "}" >> $new_jdd.data
        echo "  Fin" >> $new_jdd.data
fi
