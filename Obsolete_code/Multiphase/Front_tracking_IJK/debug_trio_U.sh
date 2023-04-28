#!/bin/bash
# Ce script sert a debugger Trio_U en parallele (ou en sequentiel).
# Il est concu pour etre execute comme une commande systeme par Trio_U, 
# dans MAIN.cpp
# Pour l'utiliser, il faut passer le chemin du script comme troisieme
# parametre de Trio_U:
#  $Mpirun -np 4 /home/myaccount/Trio_U/Trio_U_mpich datafile -mpi -debugscript=Path_to_this_script.sh

# Le script recoit le numero du processeur me() comme parametre.

#--------------------------------------------------------
# utiliser le test ci-dessous si on ne veut ouvrir le debugger que sur certains processeurs:
if test \( $1 = 0 \)
then
#	exit
echo 
fi

rm -f gdb_ready$1.gdb
# Recherche du numero de processus Trio_U:
ppid=$PPID
# Si echec de la methode ci-dessus, on peut essayer ceci:
#  pid=$$
#  ppid=`ps -l $pid | awk 'FNR==2{print $5}'`

echo Running gdb attach $ppid
# Construction d'un fichier de commandes pour gdb
cat <<END >cmd$1.gdb
 # Nom de l'executable Trio_U:
 file $exec
 # On brache gdb sur l'executable Trio_U
 attach $ppid
 # Ne pas s'arreter sur SIGUSR1 (utilise par MPI ?)
 handle SIGUSR1 pass nostop

 # Definitions pratiques:
 set print static-members 0

 # On peut placer des instructions ici (breakpoints, etc)...
 break 'Process::exit(Nom const&, int)' 

 # gdb a pris de controle de Trio_U, le script peur rendre la main a Trio_U.
 shell touch gdb_ready$1.gdb
 # Continuer l'execution de Trio_U (une fois que gdb a pris le controle,
 # il stoppe l'execution de Trio_U).
 # break Navier_Stokes_FT_Disc.cpp:622
 # break IJK_Interfaces::supprimer_duplicata_bulles()
 # break IJK_Interfaces.cpp:1076
 
cont
END
# Ouverture d'un xterm et lancement de gdb avec le fichier de commande:
xterm -T "debug PE $1" -e /usr/bin/gdb -x cmd$1.gdb &

# Il faut attendre que gdb ait pris le controle de Trio_U, sinon
# l'execution continue et on risque de rater l'evenement recherche.
# On attend l'ecriture du fichier gdb_ready par le debugger.
echo Waiting for gdb 
until [ -f gdb_ready$1.gdb ]
do
sleep 1
done
rm -f gdb_ready$1.gdb
echo Resume execution
