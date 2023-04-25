5/6/2013:

Le calcul anisotherme sur muc a un peu avance, les resultats sont franchements differents du point de depart:
profil de vitesse nettement plus asymetrique, profil de temperature plus raide aux bords, flux de chaleur
plus important.

# scp faulery@jade.cines.fr:/scratch/faulery/DNS-180/T_2/rep32/Spatiale_10.089011 ./les_jade_180_Spatiale_10.089011.txt
# scp pr3da072@supermuc.lrz.de:/gpfs/scratch/pr86bo/pr3da072/DNS_QC/valid_anisotherme/180_2/reprise0/stat*.txt .
# scp pr3da072@supermuc.lrz.de:/gpfs/scratch/pr86bo/pr3da072/DNS_QC/valid_anisotherme/180_2/tmp/statistiques_9.990006.txt .
# ../../../diff_stat.sh 10.061199 10.083895

gnuplot

set terminal postscript eps color
set output "profil_u_t.eps"
set y2tics format "%5.3f"
set grid
set key outside right
plot [][-0.5:3][][250:600]\
"les_jade_180_Spatiale_10.089011.txt" using 1:2  title "profil U LES" w p lc rgbcolor "black",\
"les_jade_180_Spatiale_10.089011.txt" using 1:9 axes x1y2 title "profil T LES" w p lc rgbcolor "black",\
"statistiques_9.990006.txt" using 1:2 w l lc rgbcolor "red"  title "profil U initial DNS",\
"statistiques_9.990006.txt" using 1:18 axes x1y2 w l lc rgbcolor "red"  title "profil T initial DNS",\
"diff_stat_0.0711988_0.0938946.txt" using 1:2 w l lc rgbcolor "blue"  title "profil U apres 0.09s DNS",\
"diff_stat_0.0711988_0.0938946.txt" using 1:18 axes x1y2 w l lc rgbcolor "blue"  title "profil T apres 0.09s DNS"

