#!/bin/bash
# source ~triou/env_triou_1.6.8.sh
# exec=/export/home/gb218285/BALTIKS/IJK/v1.6.8/DNS_IJK/New_algo_qc_opt
jdd=diffusion
# tar xvzf init.tgz

cas=DIFFUSION

rm -rf $cas
mkdir $cas
echo "############### Cas : $cas  ###############"
cd $cas
cp ../${jdd}.data .
cp ../${jdd}_VDF.data .
#ln -s ../init.lata* .

######################
# Space........
######################
rm -rf DX*

for schx in OpDiffIJK
do
   mkdir $schx 
   cd $schx
   cp ../${jdd}.data .
   cp ../${jdd}_VDF.data .
   ln -s ../init.lata* .
   echo "Schema diffusion  :   $schx ....."
   for n in 8 16 32 64 128 # 256
   do
      mkdir -p DX_EUL_$n DX_RK_$n DX_VDF_$n
      # Euler :
      sed -i -e "s/nbelem_i .*/nbelem_i $n/g" -e "s/nbelem_j .*/nbelem_j $n/g" ${jdd}.data
      echo -n "    Calculating DX_EUL_$n....."
      trust ${jdd} 1> out 2> err
      echo "Done!"
      grep "ERROR FIELD" err | awk '{print $4, $5, $6, $7}' > DX_EUL_$n/L2.txt
      cp -f ${jdd}_PP_VX.son DX_EUL_$n/
      cp -f ${jdd}_PP_VY.son DX_EUL_$n/
      # Sauvegarde du lata pour une figure a la fin: 
      if [ $n == 64 ]; then 
         mkdir -p DX_EUL_$n/lata
         cp -rf $jdd.data $jdd.lata* DX_EUL_$n/lata/
         echo "       (Copying lata to DX_EUL_$n for later use in plots...OK)"
      fi
      #
      # RK3 :
      sed -e "/time_scheme/s/#//g" ${jdd}.data > ${jdd}_RK3.data 
      echo -n "    Calculating DX_RK_$n....."
      trust ${jdd}_RK3 1> out 2> err
      echo "Done!"
      grep "ERROR FIELD" err | awk '{print $4, $5, $6, $7}' > DX_RK_$n/L2.txt
      \cp -f ${jdd}_RK3_PP_VX.son DX_RK_$n/${jdd}_PP_VX.son
      \cp -f ${jdd}_RK3_PP_VY.son DX_RK_$n/${jdd}_PP_VY.son
      # VDF : 
      np=`echo $n+1|bc`
      sed -i -e "s/Nombre_de_Noeuds .*/Nombre_de_Noeuds $np $np/g" ${jdd}_VDF.data
      echo -n "    Calculating DX_VDF_$n....."
      trust ${jdd}_VDF 1> VDF.out 2> VDF.err
      echo "Done!"
      grep "ERROR FIELD" VDF.err | awk '{print $4, $5, $6, $7}' > DX_VDF_$n/L2.txt # Fichiers Vides
      \cp -f ${jdd}_VDF_PP_VX.son DX_VDF_$n/${jdd}_PP_VX.son
      \cp -f ${jdd}_VDF_PP_VY.son DX_VDF_$n/${jdd}_PP_VY.son
   done
   #
   #
   # Comparaison a la solution analytique : 
   for sch in "VDF" "EUL" "RK" 
   do
    : > cvgx_son_$sch.txt
    : > cvgx_son_VY_$sch.txt
    : > cvgx_son_VX_$sch.txt
    : > cvgx_L2_$sch.txt
    for n in 8 16 32 64 128 # 256
    do
      for compo in "VX" "VY" 
      do 
         fic=DX_${sch}_$n/${jdd}_PP_${compo}.son
	 ficout=DX_${sch}_$n/${jdd}_${compo}_complet.son 
         val=`awk 'NR==5{print $2}' $fic`
         valx=`awk 'NR==2{print $4}' $fic`
         valy=`awk 'NR==2{print $6}' $fic`
         echo " Processing $compo $sch for $n at point $valx, $valy.  uini=$val"
         if [ $compo == "VX" ] ; then fct="cos"; 
         else fct="sin";
         fi
         awk '{x='$valx';y='$valy';T=$1;
              uana=(2.*0.001*'$fct'(x*98.17477042468103+0.2617993877991496)*'$fct'(y*98.17477042468103-0.5235987755982987)*exp(-19276.571095877651*1.e-6*T));
              print T, uana,  $2-uana}' \
              <  $fic > $ficout
         awk 'END{print '$n', $1, ($3**2)**(0.5)}' < $ficout >> cvgx_son_${compo}_$sch.txt
      done
      awk 'END{print '$n', $0}' DX_${sch}_$n/L2.txt >> cvgx_L2_$sch.txt
    done
   done
   cd ..
done

##################
# Debut du post
##################
: > plot.gplot
cat >> plot.gplot << EOF
#!/usr/bin/gnuplot
set terminal png size 640,480 enhanced font "Helvetica,12"

set output './cvgx_son.png'
set log xy
set title "abs(err) finale pour la sonde P fonction de NX"
plot   "./$schx/cvgx_son_VX_EUL.txt" u 1:3 t "$schx - EUL", \
       "./$schx/cvgx_son_VX_RK.txt" u 1:3 t "$schx - RK", \
       "./$schx/cvgx_son_VX_VDF.txt" u 1:3 t "$schx - VDF", \
       "./$schx/cvgx_son_VY_VDF.txt" u 1:3 t "$schx - VDF - Composante VY", \
       1.07343e-05*(x/8)**(-2) w l ls 2 t 'o2', \
       2.32346e-05*(x/8)**(-2) w l ls 2 t 'o2' #, \
#       0.1*(x/8)**(-1) w l t 'o1', \
 #      0.1*(x/8)**(-2) w l t 'o2', \
  #     0.1*(x/8)**(-3) w l t 'o3', \
   #    0.1*(x/8)**(-4) w l t 'o4'

set output './cvgx_L2.png'
set log xy
set title "Norme L2 fonction de NX"
plot   "./$schx/cvgx_L2_EUL.txt" u 1:3 t "$schx - EUL", \
       "./$schx/cvgx_L2_RK.txt" u 1:3 t "$schx - RK", \
        2.20918e-05*(x/8)**(-2) w l ls 2 t 'o2' #, \
#       "./$schx/cvgx_L2_VDF.txt" u 1:3 t "$schx - VDF" #, \
#       2.20918e-05*(x/8)**(-1) w l ls 1 t 'o1', \
 #      2.39984e-06*(x/8)**(-1) w l ls 1 t 'o1', \
   #    2.39984e-06*(x/8)**(-2) w l ls 2 t 'o2', \
    #   2.20918e-05*(x/8)**(-3) w l ls 3 t 'o3', \
     #  2.39984e-06*(x/8)**(-3) w l ls 3 t 'o3'


EOF

# gnuplot plot.gplot
# display ./*png

cd ..
