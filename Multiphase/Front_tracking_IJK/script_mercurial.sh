#!/bin/bash

# Script cree par B.Mathieu pour J.Dumas pour gerer ses bases mercurial
# sur differentes machines de travail.
# Principe propose par Eli Laucoin:
#  Une base centralisee sur une machine accessible de partout
#   (choix de stocker dans un espace sauvegarde). Pas de workdir dans
#   cette base.
#  Des bases sur chaque machine de travail.
#  Deux operations: push vers la base centrale, pull pour mettre a jour.
#  Le script verifie et propose les operations necessaires
#   (add des fichiers oublies, commit, merge lors du pull)

# La machine contenant la base centralisee de Jonathan est:
MACHINE_BASE_JD=gre058509
# Le repertoire contenant la base centralisee sur cette machine est:
# (sur le homedir pour que ca soit sauvegarde)
DIR_BASE_JD=/homelinux/jd235850/Maquette_ifs_central_repository
# Path mercurial sur cette machine:
MERCURIAL_PATH_BASE_JD=/work/jd235850/usr/mercurial-2.1.1/hg

# Pour chaque machine de travail de JD, nom de la machine et nom du repertoire contenant
# la copie locale de la base, path pour mercurial:

MACHINES1=gre058509
HGLOCAL1=/work/jd235850/usr/mercurial-2.1.1/hg

MACHINES2=is221705
HGLOCAL2=hg

MACHINES3=eris-ib.cluster
HGLOCAL3=hg

# Affichage d'une aide si script lance sans parametres:
if test -z "$1"
then
    cat <<EOF
options:
 pull   : met a jour la copie locale a partir de la base centrale
          (realise commit, pull, update et merge)
 push   : enregistre tous les changements dans la base centrale
          (realise commit et push)
EOF
fi

# Avant de commencer a travailler, mettre a jour la version locale avec la base centrale:
if test "$HOSTNAME" == "$MACHINES1"
    then
    LOCALHG=$HGLOCAL1
elif test "$HOSTNAME" == "$MACHINES2"
    then
    LOCALHG=$HGLOCAL2
elif test "$HOSTNAME" == "$MACHINES3"
    then
    LOCALHG=$HGLOCAL3
else
    echo Erreur: machine $HOSTNAME inconnue
    exit
fi

DIRLOCAL=`$LOCALHG root`
if test -z "$DIRLOCAL"
then
    echo Erreur: le repertoire courant n\'est pas dans une base mercurial
    exit
fi
cd $DIRLOCAL
if test "$PWD" == "$DIR_BASE_JD"
then
    echo Erreur: vous etes dans la base centrale.
    exit
fi

function test_if_commit_needed {
    echo Test s\'il y a des modifications non enregistrees dans la base locale:
    n=`$LOCALHG status -d | wc -l`
    if test $n -gt 0
	then
	echo Il y a des fichiers effaces qui sont encore dans la base\:
	echo 
	$LOCALHG status -d
	echo
	echo Il faut les retirer de la base avec \"$LOCALHG remove FICHIER\"
        echo Voulez-vous les retirer tous automatiquement maintenant \(o/n\) ?
	read reponse
	if test "$reponse" == "o"
	    then
	    $LOCALHG remove `$LOCALHG status -dn`
	    echo Commande executee: $LOCALHG remove 
	else
	    echo Non. ok je sors.
	    exit
	fi
    fi
    n=`$LOCALHG status -mar | wc -l`
    if test $n -gt 0
	then
	echo Il y a des modifications non enregistrees dans la base\:
	echo 
	$LOCALHG status -mard
	echo 
        echo Pour voir les modifications effectuees, utiliser ceci:
	echo "   "  $LOCALHG kdiff
	echo Voulez-vous que je lance la commande \"hg commit -u `whoami`\" \(o/n\) ?
	read reponse
	if test "$reponse" == "o"
	    then 
	    # entre anti-cotes, remplace le texte par le resultat de la commande
	    $LOCALHG commit -u `whoami`
	    echo Commande executee: $LOCALHG 
	else
	    echo Non. ok je sors.
	    exit
	fi
    else
	echo Non. C\'est bien.
    fi    
}

if test "$1" == "pull"
    then
    echo UPDATE_LOCAL

    test_if_commit_needed

    echo mise a jour de la copie $DIRLOCAL
    echo a partir de $DIR_BASE_JD sur $MACHINE_BASE_JD

    if test "$HOSTNAME" == "$MACHINE_BASE_JD"
	then
	$LOCALHG pull $DIR_BASE_JD 
    else
	$LOCALHG pull ssh://jd235850@$MACHINE_BASE_JD/$DIR_BASE_JD --remotecmd $MERCURIAL_PATH_BASE_JD
    fi
    echo Lancement de la commande update pour mettre a jour les fichiers non modifies en local
    $LOCALHG update
    echo Lancement de la commande merge pour resoudre les conflits entre les fichiers modifies des deux cotes
    $LOCALHG merge

elif test "$1" == "push"
    then
    echo COMMIT
    echo Liste des fichiers non enregistres:
    $LOCALHG status -u
    echo Reste-t-il des fichiers a enregistrer dans cette liste \(o/n\)?
    read reponse
    if test "$reponse" != "n"
	then
	echo Oui \? Et bien il faut  faire d\'abord  \"hg add FICHIER\".
	exit
    fi
    echo C\'est bien. 

    test_if_commit_needed

    echo On va maintenant faire push sur la base centralisee.

    if test "$HOSTNAME" == "$MACHINE_BASE_JD"
	then
	$LOCALHG push $DIR_BASE_JD 
    else
	$LOCALHG push ssh://jd235850@$MACHINE_BASE_JD/$DIR_BASE_JD --remotecmd $MERCURIAL_PATH_BASE_JD
    fi
    
else
    exit
fi
