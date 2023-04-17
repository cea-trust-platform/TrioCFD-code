#!/bin/bash
#Postraitement pour canal 2D plan : recuperation valeurs dans fichiers Nusselt et Tsortie"
tps_init_moyen=0
LC_NUMERIC=C

#Recuperation des donnees necessaires dans jdd
#---------------------------------------------
mu=`grep "mu" $1.data| awk '{print $4}' | head -1` 
rho=`grep "rho" $1.data| awk '{print $4}' | head -1`
lambda=`grep "lambda" $1.data| awk '{print $4}' | head -1` 
Cp=`grep "Cp" $1.data| awk '{print $4}' | head -1`

echo $mu > ../mu.dat
echo $rho > ../rho.dat
echo $Cp > ../Cp.dat
echo $lambda > ../lambda.dat

Prandtl=`awk "BEGIN { print $mu*$Cp/$lambda }" `

#Calcul de Tmoyen et DeltaTw theoriques 
#---------------------------------------------
X=`grep -i "longueur" $1.data| awk '{print $2}' | head -1`
Y=`grep -i "longueur" $1.data| awk '{print $3}' | head -1`
NX=`grep -i "Nombre_de_Noeuds" $1.data| awk '{print $2}' | head -1`
D=`awk "BEGIN { print $X/($NX-1)/2.} "`

phi=`grep "flux_impose" $1.data| awk '{print $5}' | head -1`

V=`grep "vitesse_imposee" $1.data| awk '{print $6}' | head -1`
Tmoyen=`awk "BEGIN {print 2*$phi*$Y/($rho*$Cp*$V*$X)}" `
tail -3 $1_pbf_Nusselt.face > fin_nusselt
# calcul du gradient de T theorique
# si Nu impose, calcul direct ; sinon on passe par Kader
logi=`grep -i "Loi_Paroi_Nu_Impose" $1.data | head -1`
if [ logi == "" ]
then
  # on recupere u* qui devrait etre different de 0, on calcule y+ et T+ avec Kader
  uetoi=`grep "u_tau" $1.data| awk '{print $7}' | head -1`
  yplus=`awk "BEGIN {print $D*$uetoi*$rho/$mu}" `
  beta=`awk "BEGIN {print (3.85*$Prandtl^(1/3)-1.3)^2+2.12*log($Prandtl)}" `
  gamma=`awk "BEGIN {print 0.01*($Prandtl*$yplus)^4/(1+5*$yplus*$Prandtl^3)}" `
  tplus=`awk "BEGIN {print $Prandtl*$yplus*exp(-$gamma)+(2.12*log(1+$yplus)+$beta)*exp(-1/$gamma)}" `
  # on en deduit DeltaTw=tplus*t_tau
  t_tau=`awk "BEGIN {print $phi/($uetoi*$Cp*$rho)}" `
  DeltaTw=`awk "BEGIN {print $t_tau*$tplus}" `
else
  Nu=`cat fin_nusselt| awk '{print $7}' | head -1`
  DeltaTw=`awk "BEGIN {print $phi*$D/($lambda*$Nu)}" `
fi
awk -v Tmoy=$Tmoyen -v dT=$DeltaTw 'BEGIN {printf ("%.4f %s %s %.4f %s %.4f\n",Tmoy,"-","-",dT,"-",dT)}'  > theoric.dat

#Lecture de Tmoyen Trio_U
#---------------------------------------------
TSmoyen=`tail -1 Tmoyen_sortie | awk '{print $2}' | head -1 `

#Lecture de Tfluide et Tparoi Trio_U
#---------------------------------------------
Tfluide=`cat fin_nusselt| awk '{print $11}' | head -1`
Tparoi=`cat fin_nusselt| awk '{print $13}' | head -1`
Tparoi_eq=`cat fin_nusselt| awk '{print $15}' | head -1`

#Calcul des deltaT Trio_U a la paroi
#---------------------------------------------
deltaT=`awk "BEGIN {print $Tparoi-$Tfluide}" `
deltaT_eq=`awk "BEGIN {print $Tparoi_eq-$Tfluide}" `

awk -v Tmoy=$TSmoyen -v Tf=$Tfluide -v Tp=$Tparoi -v Tp_eq=$Tparoi_eq -v dT=$deltaT -v dT_eq=$deltaT_eq 'BEGIN {printf ("%.4f %.4f %.4f %.4f %.4f %.4f\n",Tmoy,Tf,Tp,dT,Tp_eq,dT_eq)}'  > temp.dat

exit
