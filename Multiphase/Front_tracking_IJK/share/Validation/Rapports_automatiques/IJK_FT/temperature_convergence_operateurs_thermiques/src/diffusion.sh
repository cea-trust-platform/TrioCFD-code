#!/bin/bash
jdd=basic_temperature
cas=BASIC_TEMPERATURE

# rm -rf $cas
# mkdir $cas
#mkdir CONVERGENCE
echo "############### Cas : $cas  ###############"
#ncd CONVERGENCE
# cp ../${jdd}.data .
#cp ../${jdd}_VDF.data .

######################
# Space........
######################
rm -rf DX*

#for schx in OpDiffIJK
#do
   # mkdir $schx 
   # cd $schx
   # cp ../${jdd}.data .
    #cp ../${jdd}_VDF.data .
    # ln -s ../init.lata* .
    echo "Schema diffusion  :   QUICK PAR DEFAUT ....."
    for n in 8 16 32 64 128 256
    do
        mkdir -p DX_EUL_$n # DX_RK_$n 
        # Euler :
        #nk=`echo $n/2 | bc`
        nk=$n
        sed -i -e "s/nbelem_i .*/nbelem_i $n/g" -e "s/nbelem_j .*/nbelem_j $n/g" \
            ${jdd}.data
        #       -e "s/nbelem_k .*/nbelem_k $nk/g" ${jdd}.data
        echo -n "    Calculating DX_EUL_$n....."
        triou ${jdd} 1> out 2> err
        [ $? != 0 ] && echo "Calculation DX_EUL_$n failed! Exiting..." && exit -1
        echo "Done!"
        cp ${jdd}_PP_T.son DX_EUL_$n/
        grep "ERROR T FIELD" err | awk '{print $5, $6}' > DX_EUL_$n/L2.txt
        # cp -f ${jdd}.data ${jdd}_PP_T.son DX_EUL_$n/
        # Sauvegarde du lata pour une figure a la fin: 
        if [ $n == 256 ]; then 
            mkdir -p DX_EUL_$n/lata
            cp -rf $jdd.data $jdd.lata* DX_EUL_$n/lata/
            echo "       (Copying lata to DX_EUL_$n for later use in plots...OK)"
        fi
        # calcul en 3D pour la stabilite : 
        if [ $n == 32 ]; then 
            mkdir -p DX_EUL_${n}_3D
            sed -e "s/nbelem_i .*/nbelem_i $n/g" -e "s/nbelem_j .*/nbelem_j $n/g" \
                -e "s#cos(y\*2\*Pi/0.06)#cos(y\*2\*Pi/0.06)*cos(z\*2\*Pi/0.06)#g" \
                -e "s#2\*(2\*Pi/0.06)^2)#3\*(2\*Pi/0.06)^2)#g" \
                -e "s/nbelem_k .*/nbelem_k $n/g" \
                ${jdd}.data > ${jdd}_3D.data
            echo -n "    Calculating DX_EUL_${n}_3D....."
            triou ${jdd}_3D 1> out 2> err
            [ $? != 0 ] && echo "Calculation DX_EUL_${n}_3D failed! Exiting..." && exit -1
            echo "Done!"
            grep "ERROR T FIELD" err | awk '{print $5, $6}' > DX_EUL_${n}_3D/L2.txt
            cp -f ${jdd}_3D.data ${jdd}_3D_PP_T.son DX_EUL_${n}_3D/ 
        fi
        #
        # RK3 :
        #sed -e "/time_scheme/s/#//g" ${jdd}.data > ${jdd}_RK3.data 
        #echo -n "    Calculating DX_RK_$n....."
        #triou ${jdd}_RK3 1> out 2> err
        #echo "Done!"
        #grep "ERROR T FIELD" err | awk '{print $4, $5, $6, $7}' > DX_RK_$n/L2.txt
        #\cp -f ${jdd}_RK3_PP_T.son DX_RK_$n/${jdd}_T_VX.son
    done
    #
    #
    # Comparaison a la solution analytique : 
    for sch in "EUL" # "RK"  "VDF" 
    do
        echo "Post traitement pour $sch "
        : > cvgx_son_$sch.txt
        : > cvgx_son_T_$sch.txt
        : > cvgx_L2_$sch.txt
        for n in 8 16 32 64 128  256
        do
            for compo in "T" 
            do 
                fic=DX_${sch}_$n/${jdd}_PP_${compo}.son
                ficout=DX_${sch}_$n/${jdd}_${compo}_complet.son 
                val=`awk 'NR==5{print $2}' $fic`
                valx=`awk 'NR==2{print $4}' $fic`
                valy=`awk 'NR==2{print $6}' $fic`
                echo " Processing $compo $sch for $n at point $valx, $valy.  Tini=$val"
                awk '{x='$valx';y='$valy';Time=$1;
                Tana=(0.5*exp(-(0.1/40.*2*(2*Pi/0.06)^2)*Time)*cos(x*2*Pi/0.06)*cos(ys*2*Pi/0.06));
                print Time, Tana,  $2-Tana}' \
                    <  $fic > $ficout
                awk 'END{print '$n', $1, ($3**2)**(0.5)}' < $ficout >> cvgx_son_${compo}_$sch.txt
            done
            awk 'END{print '$n', $0}' DX_${sch}_$n/L2.txt >> cvgx_L2_$sch.txt
        done
    done
 #    cd ..
 #done

##################
# Debut du post
##################
echo "Début post-traitement"
echo "---------------------\n\n"
: > plot.gplot
cat >> plot.gplot << EOF
# #!/usr/bin/gnuplot
set terminal png size 640,480 enhanced font "Helvetica,12"

#set output './cvgx_son.png'
#set log xy
#set title "abs(err) finale pour la sonde T fonction de NX"
#plot   "./$schx/cvgx_son_T_EUL.txt" u 1:3 t "$schx - EUL", \
    #       1.07343e-05*(x/8)**(-2) w l ls 2 t 'o2', \
    #       2.32346e-05*(x/8)**(-2) w l ls 2 t 'o2' #, \
    #       0.1*(x/8)**(-1) w l t 'o1', \
    #      0.1*(x/8)**(-2) w l t 'o2', \
    #     0.1*(x/8)**(-3) w l t 'o3', \
    #    0.1*(x/8)**(-4) w l t 'o4'
#       "./$schx/cvgx_son_VX_VDF.txt" u 1:3 t "$schx - VDF", \
    #       "./$schx/cvgx_son_VY_VDF.txt" u 1:3 t "$schx - VDF - Composante VY", \

set output './cvgx_L2.png'
set log xy
set title "Norme L2 fonction de NX"
plot   "./cvgx_L2_EUL.txt" u 1:3 t "PAR_DEF - EUL", \
    0.00440539*10**(-1)*(x/8)**(-2) w l ls 2 t 'o2t '
    # 0.00440539*(x/8)**(-1) w l ls 3 t 'o1' #, \
    #       "./$schx/cvgx_L2_RK.txt" u 1:3 t "$schx - RK", \
    #       "./$schx/cvgx_L2_VDF.txt" u 1:3 t "$schx - VDF" #, \
    #       2.20918e-05*(x/8)**(-1) w l ls 1 t 'o1', \
    #      2.39984e-06*(x/8)**(-1) w l ls 1 t 'o1', \
    #    2.39984e-06*(x/8)**(-2) w l ls 2 t 'o2', \
    #   2.20918e-05*(x/8)**(-3) w l ls 3 t 'o3', \
    #  2.39984e-06*(x/8)**(-3) w l ls 3 t 'o3'


EOF

export LC_ALL="en_US.UTF-8"
gnuplot plot.gplot
# display ./*png
rm *.lata*
# cd ..
