# FTD_2D_Axi
echo Chute d\'une goutte sous gravite g=[0. 0. -10.]
cat <<END | gnuplot
set output "compare_sonde.eps"
set terminal postscript color eps
plot \
-10.*x title "vitesse theorique",\
"prepare_V_GOUTTE_reference.son" using 1:4 title "calcul de reference" ps 2,\
"prepare_V_GOUTTE.son" using 1:4 title "resultat test prepare.data",\
"FTD_reprise_xyz_vef_3d_V_GOUTTE.son" using 1:4 title "resultat test seq",\
"FTD_reprise_xyz_vef_3d_V_GOUTTE.son" using 1:4 title "resultat test par"
END
echo Pour voir le fichier resultat: gv compare_sonde.eps
gv compare_sonde.eps
