!/bin/sh
#Pos traitement pour calcul conduction
LC_NUMERIC=C

# frequence de battements
f=2
pi=3.1416
# Nombre de points de la courbe trms
nb_trms=60.
# Temps de simulation pour la solution
temps_simu=10.
# Nombre de points de la courbe f(t)
nb_temps_simu=800
# Nombre de courbes f(t)
nb_courb=13.

#Recuperation des donnees necessaires dans jdd
#---------------------------------------------
rho=`grep "rho" $1.data| awk '{print $4}' | head -1`
lambda=`grep "lambda" $1.data| awk '{print $4}' | head -1` 
Cp=`grep "Cp" $1.data| awk '{print $4}' | head -1`
L=`grep "Longueurs" $1.data| awk '{print $2}' | head -1`

alpha=`awk "BEGIN { print $lambda/$rho/$Cp }" `

echo "rho=" $rho
echo "lambda=" $lambda
echo "Cp=" $Cp
echo "L=" $L
echo "alpha=" $alpha

echo $mu > ../mu.dat
echo $rho > ../rho.dat
echo $Cp > ../Cp.dat
echo $lambda > ../lambda.dat

# Amplitude
A=`grep "Gauche   paroi_temperature_imposee champ_front_fonc_txyz" $1.data| awk '{print substr($5,1,4)}'| head -1`
echo "amplitude=" $A

# Frequence de battement

omega=`awk "BEGIN { print 2*$pi*$f }" `
echo "omega=" $omega

# Coefficient de decroissance
coef=`awk "BEGIN { print -sqrt($omega/2./$alpha) }" `
echo "coef=" $coef

# Courbe Trms
#------------
# Premier point sur la condition limite
x=`awk "BEGIN { print 0. }" `
Trms=`awk "BEGIN { print $A/1.41421*exp($x*$coef) }" `
awk -v x=$x -v Trms=$Trms  'BEGIN {printf ("%.5f %.5f\n",x,Trms)}'  > Trms_theo.dat


for n in $(seq $nb_trms)
do
  x=`awk "BEGIN { print $n*$L/$nb_trms }" `
  Trms=`awk "BEGIN { print $A/1.41421*exp($x*$coef) }" `
#  echo "n="$n
  awk -v x=$x -v Trms=$Trms  'BEGIN {printf ("%.5f %.5f\n",x,Trms)}'  >>  Trms_theo.dat
done

# Courbes en fonction du temps.
#----------------------------


x=`awk "BEGIN { print 0. }" `

n=1

for n in $(seq $nb_courb)
do
  x=`awk "BEGIN { print (($n-1)*$L/($nb_courb-1)) }" `
  echo "courbe n="$n"  x="$x
  t=0.
  n_simu=0
  dt=`awk "BEGIN { print $temps_simu/$nb_temps_simu }" `

  T_theo=`awk "BEGIN { print $A*exp($x*$coef)*sin($omega*$t+$coef*$x) }" `
  echo "t="$t" T_theo="$T_theo
  fic_sor=`echo 'T_theo'$x'.dat'`
  awk -v t=$t -v T_theo=$T_theo  'BEGIN {printf ("%.8f %.8f\n",t,T_theo)}'  >  $fic_sor
  
  while [ $n_simu -le $nb_temps_simu ]
  do
    T_theo=`awk "BEGIN { print $A*exp($x*$coef)*sin($omega*$t+$coef*$x) }" `
#    echo "n_simu=" $n_simu "  t="$t"  T_theo="$T_theo
    t=`awk "BEGIN { print $t+$dt }" `
    n_simu=`awk "BEGIN { print ($n_simu+1) }" `
    awk -v t=$t -v T_theo=$T_theo  'BEGIN {printf ("%.8f %.8f\n",t,T_theo)}'  >>  $fic_sor
 done


done