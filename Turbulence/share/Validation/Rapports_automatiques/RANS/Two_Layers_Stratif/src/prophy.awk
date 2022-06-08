#!/bin/bash
#Postraitement : recuperation valeurs des proprietes physiques
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

exit
