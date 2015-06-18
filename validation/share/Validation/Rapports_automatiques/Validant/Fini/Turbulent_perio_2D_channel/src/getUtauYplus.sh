#!/bin/bash
#recupere les valeurs de utau & yplus, en fin de calculs (PreRequis)

RACINE=`pwd`

echo "Running getUtauYplus.sh in $RACINE"

for dom in N0 N1 N2 N3 N4 N5  ; do
	cd $RACINE/$dom/
	python ..//getUTau.py  10.15 ./Amont ./Muscl ./EF_stab ./EF_stab02
	python ../getYPlus.py ./Amont ./Muscl ./EF_stab ./EF_stab02
	#mv yplus.dat yplus_$dom.dat
	#mv utau.dat utau_$dom.dat
done;

cd $RACINE/
rm -f frictionVelocity.dat
cat N0/utau.dat >  frictionVelocity.dat
cat N1/utau.dat >> frictionVelocity.dat
cat N2/utau.dat >> frictionVelocity.dat
cat N3/utau.dat >> frictionVelocity.dat
cat N4/utau.dat >> frictionVelocity.dat
cat N5/utau.dat >> frictionVelocity.dat

