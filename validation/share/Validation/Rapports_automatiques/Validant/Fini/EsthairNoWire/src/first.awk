#!/bin/bash
#Postraitement pour canal 2D plan : recuperation valeurs dans fichiers Nusselt et Tsortie"
tps_init_moyen=0
LC_NUMERIC=C

#Recuperation des donnes utiles dans jdd
#----------------------------------------
mu=`grep "mu" $1.data| awk '{print $4}' | head -1` 
rho=`grep "rho" $1.data| awk '{print $4}' | head -1`

echo $mu > mu.dat
echo $rho > rho.dat

Ubulk=`grep "vitesse Champ_Uniforme" $1.data| awk '{print $6}' | head -1`

# Valeurs figees liees a la geometrie Esthair
#---------------------------------------------
Dhy_ass=0.0112
Re_ass=`awk "BEGIN {print $Dhy_ass*$Ubulk*$rho/$mu}"`

# Recuperation de donnees dans le fichier err
#---------------------------------------------
Vtotal=`grep "sum(volume" *err| awk '{print $NF}' | head -1`
swall=`grep "Area of thi" *err| awk '{print $5}' | head -1`
sperio=`grep "Area of perio" *err| awk '{print $5}' | head -1`
sperio=`awk "BEGIN {print $sperio/2.} " `
srod=`grep "Area of aiguilles" *err| awk '{print $5}' | head -1`

DhyTrio=`awk "BEGIN { print $Vtotal/($swall+$srod)*4.} "`
 
# Recuperation de donnees dans le fichier Channel_Flow_Rate_repr_pb_perio
#------------------------------------------------------------------
Channel_Flow_Rate=`ls *Channel_Flow_Rate_repr_pb_perio`
DebiTrio=`cat $Channel_Flow_Rate| awk '{print $2}' | tail -1`
UTrio=`awk "BEGIN { print $DebiTrio/$sperio} "`
ReTrio=`awk "BEGIN { print $UTrio*$DhyTrio*$rho/$mu} "`

awk -v Dhy=$Dhy_ass -v DhyT=$DhyTrio -v Re=$Re_ass -v ReT=$ReTrio -v U=$UTrio 'BEGIN {printf ("%.4f %.4f %.1f %.1f %.1f %.1f %.2f\n",Dhy,DhyT,(DhyT-Dhy)/Dhy*100,Re,ReT,(ReT-Re)/Re*100,U)}'  > reyno.dat

exit
