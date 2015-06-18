#!/bin/bash

usage() {
	echo "UTILISATION : adim.sh nom.data fic_Moyenne u_tau "
	exit -1
	}
	
if [ $# == 0 ] 
then
 usage
fi

mu=`grep "mu " $1 | tail -1 | awk '{print $4}'`
rho=`grep "rho" $1 | tail -1 | awk '{print $4}'`

awk -v utau=$3 -v mu=$mu -v rho=$rho 'BEGIN{nu=mu/rho;} {if (($1!="")&&($1!="#")) print $1*utau/nu " " $2/utau " " $3/utau " " $4/utau " " $5/utau " " $6/utau " " $7/utau " " $8/utau/utau " " $9/utau/utau " " $10/utau/utau " ";}' $2 > $2_adim.dat

col=""
block=""
tab=("-um" "-vm" "-wm" "-up" "-vp" "-wp" "-uv" "-uw" "-vw")


if [ $# -gt 3 ] 
then

log=0
for i in $*
do
ind=2
is_col=0
	for j in ${tab[*]}
	do
	if [ $i == $j ]
		then
		col=$col" -bxy 1:$ind"
		is_col=1
		if [ $ind = 2 ] 
		then 
			log=1
		fi
	fi
	let $[ind+=1]
	done
done
block=$block" -block $2_adim.dat $col"
if [ $log -eq 1 ]
then
block=$block" $HOME/local/bin/log.dat"
fi
echo $block
xmgrace -log x $block 

fi
