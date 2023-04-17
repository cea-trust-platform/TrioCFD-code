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
#echo "mu"=$mu
#echo "rho"=$rho

#Recuperation des dimensions necessaires : en dur
#---------------------------------------------
h=0.05
H=0.1
L1=0.00
L2=1.2
S1=`awk "BEGIN {print $h*$h}" ` 
S2=`awk "BEGIN {print $H*$H}" ` 
#echo "S1"=$S1
#echo "S2"=$S2
V1=`grep "velocity" $1.data| awk '{print $3}' | head -1`
V2=`awk "BEGIN {print $V1*$S1/$S2}" ` 

echo "V1=" $V1
echo "V2=" $V2

# Diametre hydraulique
#---------------------
DH1=`awk "BEGIN {print 4*$h*$h/(4*$h)}" `
DH2=`awk "BEGIN {print 4*$H*$H/(4*$H)}" `

#echo "DH1=" $DH1
#echo "DH2=" $DH2


#Nombre de reynolds 
#---------------------------------------------
Re2=`awk "BEGIN { print $rho*$V2*$DH2/$mu }" `
Re1=`awk "BEGIN { print $rho*$V1*$DH1/$mu }" `

#echo "Re1=" $Re1
#echo "Re2=" $Re2

# coefficients de frottement et perte de charge theoriques
#---------------------------------------------
f1_lam=`awk "BEGIN { print 64/$Re1 }" `
f2_lam=`awk "BEGIN { print 64/$Re2 }" `
f1_tur=`awk "BEGIN {print (0.316/$Re1^0.25)*$L1/$DH1}" `
f2_tur=`awk "BEGIN {print (0.316/$Re2^0.25)*$L2/$DH2}" `
dpfrot1=`awk "BEGIN {print $f1_tur*$rho*$V1*$V1/2}" `
dpfrot2=`awk "BEGIN {print $f2_tur*$rho*$V2*$V2/2}" `

k=`awk "BEGIN {print (1-$S1/$S2)^2}" `
dpfrot_exp=`awk "BEGIN {print $k*$rho*$V1*$V1/2}" `

dp_tot_theo=`awk "BEGIN {print $dpfrot1+$dpfrot2+$dpfrot_exp}" `

#echo "f1tur="$f1_tur 
#echo "f2tur="$f2_tur 
#echo "dpfrot1="$dpfrot1
#echo "dpfrot2="$dpfrot2
#echo "k="$k
#echo "P1_sur_rho_trio=" $P1_sur_rho_trio
#echo "P1" $P1


# Perte de charge calculee
#---------------------------------------------

head -2  $1_SONDE_PRESSION.coupe > delta_p1
tail -1  delta_p1 > delta_p
P1_sur_rho_trio=`cat delta_p| awk '{print $2}' | head -1`
P1=`awk "BEGIN {print $P1_sur_rho_trio*$rho}" `
P2=0
dp_trio=`awk "BEGIN {print $P1+(-1*$P2)+$rho*($V1*$V1-$V2*$V2)/2}" `

# Pourcentage erreur
#---------------------------------------------
erreur=`awk "BEGIN {print 100*($dp_trio+(-1)*$dp_tot_theo)/$dp_trio}" `


echo "dp_tot_theo=" $dp_tot_theo
#echo "P1_sur_rho_trio=" $P1_sur_rho_trio
#echo "P1" $P1

echo "dp_trio=" $dp_trio
echo "erreur=" $erreur


awk -v x1=$V1 -v x2=$Re1 -v x3=$dp_tot_theo -v x4=$dp_trio -v x5=$erreur  'BEGIN {printf ("%.1f %.1f %.2f %.2f %.2f\n",x1,x2,x3,x4,x5)}'  > delta_p.dat
cat delta_p.dat >> ../deltap_courbe_vef.dat


exit
