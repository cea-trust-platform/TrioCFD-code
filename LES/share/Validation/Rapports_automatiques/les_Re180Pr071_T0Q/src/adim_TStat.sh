#!/bin/bash

usage() {
	echo "UTILISATION : adim.sh nom.data fic_Moyenne u_tau t_tau t_w"
	exit -1
	}
	
if [ $# == 0 ] 
then
 usage
fi

mu=`grep "mu " $1 | awk '{print $4}'`
rho=`grep "rho" $1 | awk '{print $4}'`
Cp=`grep "Cp" $1 | awk '{print $4}'`

awk -v utau=$3 -v ttau=$4 -v tw=$5 -v mu=$mu -v rho=$rho -v Cp=$Cp 'BEGIN{nu=mu/rho;  ad=ttau ; } {if (($1!="")&&($1!="#")) print $1*utau/nu " " ($2-tw)/ad " " $3/ad " " $4/ad/utau " " $5/ad/utau " " $6/ad/utau " ";}' $2 > $2_adim.dat




