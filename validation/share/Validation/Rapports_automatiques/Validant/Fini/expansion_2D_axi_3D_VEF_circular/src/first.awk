#!/bin/bash
#Postraitement pour expension brusque "
tps_init_moyen=0
LC_NUMERIC=C

#Recuperation des donnees necessaires dans jdd
#---------------------------------------------
mu=`grep "mu" $1.data| awk '{print $4}' | head -1` 
rho=`grep "rho" $1.data| awk '{print $4}' | head -1`
#echo $mu > mu.dat
#echo $rho > rho.dat
#echo $mu
#echo $rho

#Recuperation des rayons necessaires : en dur
#---------------------------------------------
R1=0.257
R2=0.5
L1=0.25
L2=5.75
D1=`awk "BEGIN {print 2*$R1}" ` 
D2=`awk "BEGIN {print 2*$R2}" ` 
V1=`grep "Entree frontiere_ouverte_vitesse_imposee Champ_Front_Uniforme" $1.data| awk '{print $6}' | head -1`
V2=`awk "BEGIN {print $V1*$D1*$D1/$D2/$D2}" ` 

echo "V1=" $V1
echo "V2=" $V2



#Nombre de reynolds 
#---------------------------------------------
Re2=`awk "BEGIN { print $rho*$V2*$D2/$mu }" `
Re1=`awk "BEGIN { print $rho*$V1*$D1/$mu }" `


# coefficients de frottement et perte de charge theoriques
#---------------------------------------------
f1_lam=`awk "BEGIN { print 64/$Re1 }" `
f2_lam=`awk "BEGIN { print 64/$Re2 }" `
f1_tur=`awk "BEGIN {print (0.316/$Re1^0.25)*$L1/$D1}" `
f2_tur=`awk "BEGIN {print (0.316/$Re2^0.25)*$L2/$D2}" `
dpfrot1=`awk "BEGIN {print $f1_tur*$rho*$V1*$V1/2}" `
dpfrot2=`awk "BEGIN {print $f2_tur*$rho*$V2*$V2/2}" `

k=`awk "BEGIN {print (1-$D1*$D1/($D2*$D2))^2}" `
dpfrot_exp=`awk "BEGIN {print $k*$rho*$V1*$V1/2}" `

dp_tot_theo=`awk "BEGIN {print $dpfrot1+$dpfrot2+$dpfrot_exp}" `


# Perte de charge calculee
#---------------------------------------------

head -2  $1_SONDE_PRESSION.coupe > delta_p1
tail -1  delta_p1 > delta_p
P1_sur_rho_trio=` awk '{print $2}' delta_p | head -1`
P1=`awk "BEGIN {print $P1_sur_rho_trio*$rho}" `
P2=0
dp_trio=`awk "BEGIN {print $P1+(-1*$P2)+$rho*($V1*$V1-$V2*$V2)/2}" `

# Pourcentage erreur
#---------------------------------------------
erreur=`awk "BEGIN {print 100*($dp_trio+(-1)*$dp_tot_theo)/$dp_trio}" `


#echo "Re1="$Re1 
#echo "Re2="$Re2 
#echo "dp_tot_theo=" $dp_tot_theo
#echo "P1_sur_rho_trio=" $P1_sur_rho_trio
#echo "P1" $P1

#echo "dp_trio=" $dp_trio
#echo "erreur=" $erreur


awk -v x1=$V1 -v x2=$Re1 -v x3=$dp_tot_theo -v x4=$dp_trio -v x5=$erreur  'BEGIN {printf ("%.1f %.1f %.2f %.2f %.2f\n",x1,x2,x3,x4,x5)}'  > delta_p.dat
cat delta_p.dat >> ../deltap_courbe2d.dat


exit
