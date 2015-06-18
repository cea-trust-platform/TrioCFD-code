#!/bin/bash
tps_init_moyen=0
LC_NUMERIC=C

# Recuperation des resultats des coeff de frottement
#---------------------------------------------
Yplus=`tail -n 3 EsthydSf_pb_Ustar.face | awk '{print $4}' | head -n 1`
PGTrio=`cat *Channel_Flow_Rate_repr_pb_perio| awk '{print $4}' |tail -1`
UTrio=`cat reyno.dat| awk '{print $7}' | head -n 1`
DhyTrio=`cat reyno.dat| awk '{print $2}' | head -n 1`
CfTrio=`awk "BEGIN { print $PGTrio*$DhyTrio/($UTrio*$UTrio)/2 } "`
ReAss=`cat reyno.dat| awk '{print $4}' | head -n 1`
CfBlasius=`awk "BEGIN { print $ReAss^(-0.25)*0.316/4 } "`

awk -v Yp=$Yplus -v Re=$ReAss -v PG=$PGTrio -v CfT=$CfTrio 'BEGIN {printf ("%.1f %.1f %.2f %.5f\n",Re,Yp,PG,CfT)}'  > press.dat
awk -v Re=$ReAss -v CfT=$CfTrio -v CfB=$CfBlasius 'BEGIN {printf ("%.1f %.5f %.5f %.1f\n",Re,CfT,CfB,(CfT-CfB)/CfB*100)}'  > frot.dat

exit
