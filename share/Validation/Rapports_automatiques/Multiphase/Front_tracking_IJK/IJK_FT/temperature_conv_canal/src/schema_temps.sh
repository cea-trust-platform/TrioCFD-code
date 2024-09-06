#!/bin/bash
# source ~triou/env_triou_1.6.8.sh
# exec=/export/home/gb218285/BALTIKS/IJK/v1.6.8/DNS_IJK/New_algo_qc_opt
jdd=schema_temps
# tar xvzf init.tgz

######################
# Space........
######################
rm -rf DX*

schx=SCHEMA_TEMPS

rm -rf $schx
mkdir $schx
echo "############### Cas : $schx  ###############"
cd $schx
cp -f ../${jdd}.data .
#ln -s ../init.lata* .

######################
# Time........
######################
rm -rf DT*
nx=5
sed -i -e "s/nbelem_i .*/nbelem_i $nx/g" -e "s/nbelem_j .*/nbelem_j $nx/g" ${jdd}.data
\cp -f ${jdd}.data ${jdd}.data.ori
for n in 1 2 4 8 16 32 64 128 256 512 
do
   mkdir -p DT_EUL_$n DT_RK_$n
   nit=`echo 2048/$n|bc`
   dt=`echo 0.001953125*$n|bc` # Soit un tfin de 4s
   # Euler :
   sed -i -e "s/nb_pas_dt_max .*/nb_pas_dt_max $nit/g" -e "s/timestep .*/timestep $dt/g" ${jdd}.data
   echo -n "Calculating DT_EUL_$n....."
   trust ${jdd} 1> out 2> err
   echo "Done!"
   \cp -f ${jdd}_acceleration.out DT_EUL_$n/acc.txt
   grep "ERROR FIELD" err | awk '{print $4, $5, $6, $7}' > DT_EUL_$n/L2.txt
   \cp -f ${jdd}_P.son DT_EUL_$n/
   #
   # Sauvegarde du lata pour une figure a la fin: 
   if [ $n == 64 ]; then 
      mkdir -p DT_EUL_$n/lata
      cp -rf ${jdd}.data $jdd.lata* DT_EUL_$n/lata/
      echo "   (Copying lata to DT_EUL_$n for later use in plots...OK)"
   fi
   #
   # RK3 :
   sed -e "/time_scheme/s/#//g" ${jdd}.data > ${jdd}_RK3.data 
   echo -n "Calculating DT_RK_$n....."
   trust ${jdd}_RK3 1> out 2> err
   echo "Done!"
   \cp -f ${jdd}_RK3_acceleration.out DT_RK_$n/acc.txt
   grep "ERROR FIELD" err | awk '{print $4, $5, $6, $7}' > DT_RK_$n/L2.txt
   \cp -f ${jdd}_RK3_P.son DT_RK_$n/${jdd}_P.son
done
\cp -f ${jdd}.data.ori ${jdd}.data

# Comparaison a la solution analytique : 
for sch in "EUL" "RK"
do
 : > cvgt_son_$sch.txt
 : > cvgt_L2_$sch.txt
 : > cvgt_acc_$sch.txt
 for n in 1 2 4 8 16 32 64 128 256 512
 do
   awk '{
	 print $1, $2, 1.-(cos($1)+sin($1)), $2-(1.-(cos($1)+sin($1)));
	 }' <  DT_${sch}_$n/${jdd}_P.son > DT_${sch}_$n/${jdd}_Pcomplet.son 
   awk 'END{print '$n', $1, ($4**2)**(0.5)}' < DT_${sch}_$n/${jdd}_Pcomplet.son >> cvgt_son_$sch.txt
   awk 'END{print '$n', $0}' DT_${sch}_$n/L2.txt >> cvgt_L2_$sch.txt
   # Calcul des erreurs sur la vitesse moyenne 
   # de l'Erreur sur rho*acc  et sur d(rho*acc)/dt (qu'on impose)
   # pas tous au meme temps... acceleration est a tnew...  
   awk 'END{print '$n', $2, 
   		   (($3-(1.-(cos($2)+sin($2))))**2)**(0.5),
		   (($4-500.*(cos($2)+sin($2)))**2)**(0.5), 
		   $5,
		   (($8-500.*(-cos($7)+sin($7)))**2)**(0.5)
		   }' DT_${sch}_$n/acc.txt >> cvgt_acc_$sch.txt
 done
done


##################
# Debut du post
##################
cat >> plot.gplot << EOF
#!/usr/bin/gnuplot
set terminal png large
set output './cvgt_son.png'
set log xy
set xlabel "DT/0.001953125"
set title "abs(err) finale pour la sonde P en fonction de DT sur maillage $nx"
plot   "./cvgt_son_EUL.txt" u 1:3 t "EUL", \
       "./cvgt_son_RK.txt" u 1:3 t "RK", \
       3e-3*x**(1) w l t 'o1', \
       3e-11*x**(3) w l t 'o3'

set output './cvgt_L2.png'
set log xy
set xlabel "DT/0.001953125"
set title "Norme L2 fonction de DT sur maillage $nx"
plot   "./cvgt_L2_EUL.txt" u 1:3 t "EUL", \
       "./cvgt_L2_RK.txt" u 1:3 t "RK", \
       3.9e-03*x**(1) w l t 'o1', \
       3.9e-11*x**(3) w l t 'o3'


set output './cvgt_acc.png'
set log xy
set xlabel "DT/0.001953125"
set key bottom right
set title "acc fonction de DT"
plot   "./cvgt_acc_EUL.txt" u 1:6 t "EUL - rho*acc", \
       "./cvgt_acc_RK.txt" u 1:6 t "RK - rho*acc", \
       "./cvgt_acc_EUL.txt" u 1:4 t "EUL - d(rho*acc)/dt", \
       "./cvgt_acc_RK.txt" u 1:4 t "RK - d(rho*acc)/dt", \
       0.198*x**(1) w l t 'o1', \
       1.5e-08*x**(3) w l t 'o3'
#       "./cvgt_acc_EUL.txt" u 1:3 t "EUL - vmoy", \
#       "./cvgt_acc_RK.txt" u 1:3 t "RK - vmoy", \


EOF

# gnuplot plot.gplot
# display ../*png
cd ..
