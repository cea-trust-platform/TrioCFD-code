Fiche basic_temperature
-----------------------

DESCRIPTION : Cette fiche n'est qu'une copie de la fiche monte. La différence est l'ajout d'un champ de température.
ETAT_RUN : KO ?
ETAT_PRM : KO
PB : probleme lancement \*_fine, voir prepare.sh
ACTION : Suppression ?

Fiche basic_temeprature_test
----------------------------

DESCRIPTION : Cette fiche reprend les différents cas de conditions aux limites dans le cas de la reprise et du calcul parallèle. 
ETAT_RUN : OK
ETAT_PRM : KO
PB : 
           
   Post traitement pour EUL 
    Processing T EUL for 8 at point 2.62500000e-03, 2.62500000e-03.  Tini=-1.43440336e-02
    Processing T EUL for 16 at point 2.81250000e-03, 2.81250000e-03.  Tini=-1.43440336e-02
    Processing T EUL for 32 at point 2.90625000e-03, 2.90625000e-03.  Tini=-1.43440336e-02
    Processing T EUL for 64 at point 2.95312500e-03, 2.95312500e-03.  Tini=-1.43440336e-02
   line 0: warning: iconv failed to convert degree sign
   
   set terminal png size 640,480 enhanced font "Helvetica,12"
                ^
   "plot.gplot" line 2: unknown or ambiguous terminal type; type just 'set terminal' for a list
   
   ==============================
   -> Building PDF report file...
   ==============================
   Traceback (most recent call last):
     File "/home/asonolet/Trio/TRUST-1.8.0/Validation/Outils/Genere_courbe/src/genererCourbes.py", line 521, in <module>
       app.genererRapport(debug_figure)
     File "/home/asonolet/Trio/TRUST-1.8.0/Validation/Outils/Genere_courbe/src/genererCourbes.py", line 367, in genererRa
   pport
       indice = chapitre.genererGraphes(destTMP, indice,debug_figure, self.novisit)
     File "/home/asonolet/Trio/TRUST-1.8.0/Validation/Outils/Genere_courbe/src/Chapitre.py", line 115, in genererGraphes
       figure.genererGraphe(dest, indice,debug_figure)
     File "/home/asonolet/Trio/TRUST-1.8.0/Validation/Outils/Genere_courbe/src/Visu.py", line 230, in genererGraphe
       ficPlot.write('ok=AddPlot(\"Mesh\",\"%s\")\n'%plot[2])
   IndexError: list index out of range
   Error, see Run.log

Fiche stat_temperature
----------------------

DESCRIPTION : reprise d'un fichier lata avec bulles
ETAT_RUN : OK
ETAT_PRM : OK
PB : Segmentation fault lors du calcul... ? Résolu avec le nouvel exec !



Fiche temperature_conv
----------------------

DESCRIPTION : Etude de la convergence des schémas en temps, du schéma Quick de diffusion et du schéma de convection par défaut
ETAT_RUN : OK
ETAT_PRM : KO
PB : set terminal ... c'est probablement un problème machine. NON, c'est probablement un problème de PRM
- [X] 2020-06-19 régler le PRM

Fiche temperature_conv_canal
----------------------------

DESCRIPTION : Cette fiche est une copie de temperature conv, à la différence près que la convergence ne se fait pas pour un SWARM (Perio k) mais pour un canal avec des conditions de température aux limites (en k)
- [X] 2020-06-19 régler le PRM

Fiche temperature_diphasique
----------------------------

DESCRIPTION : Cette fiche étudie le comportement du solveur de température en présence de bulles, dans les cas où le calcul est séquentiel, parallèle ou repris.
ETAT_RUN : OK
ETAT_PRM : OK

Fiche temperature_monop
-----------------------

DESCRIPTION : Une copie de une copie de temperature_diphasique mais sans les bulles. En canal, similaire à temperature_conv_canal mais là on étudie les reprise/parallèle et non les opérateurs.

ETAT_RUN : OK
ETAT_PRM : OK

Fiche temperature_stat
----------------------

DESCRIPTION : Dans cette fiche, on lance un calcul simple avec les statistiques de la thermique. Pour le moment il n'y a aucun plot des statistiques, tout est vérifié à la main. Les cas de reprise et de parallèle sont traités.
ETAT_RUN : OK
PB : Segmentation fault
- [X] 2020-06-18 vérifier que ça marche avec le nouvel exec
