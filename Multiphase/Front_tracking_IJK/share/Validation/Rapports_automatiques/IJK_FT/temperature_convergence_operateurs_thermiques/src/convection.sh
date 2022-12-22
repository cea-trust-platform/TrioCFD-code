#!/bin/bash

jdd=CONVECTION 
cas=CONVECTION

echo "############### Cas : $cas  ###############"

######################
# Space........
######################
rm -rf DX*

for n in 8 16 32 64 128 256
do
    mkdir -p DX_EUL_$n DX_RK_$n
    # Euler :
    sed -i -e "s/nbelem_i .*/nbelem_i $n/g" -e "s/nbelem_j .*/nbelem_j $n/g" ${jdd}.data
    echo -n "    Calculating DX_EUL_$n....."
    triou ${jdd} 1> out 2> err
    echo "Done!"
    grep "ERROR T FIELD" err | awk '{print $5, $6}' > DX_EUL_$n/L2.txt
    cp ${jdd}_PP_T.son DX_EUL_$n/

    # Sauvegarde du lata pour une figure a la fin: 
    if [ $n == 256 ]; then 
        mkdir -p DX_EUL_$n/lata
        cp -rf $jdd.data $jdd.lata* DX_EUL_$n/lata/
        echo "       (Copying lata to DX_EUL_$n for later use in plots...OK)"
    fi

    # RK3 :
    sed -e "/time_scheme/s/#//g" ${jdd}.data > ${jdd}_RK3.data 
    echo -n "    Calculating DX_RK_$n....."
    triou ${jdd}_RK3 1> out 2> err
    [ $? != 0 ] && echo "Calculation DX_RK_${n} failed! Exiting..." && exit -1
    echo "Done!"
    grep "ERROR T FIELD" err | awk '{print $5, $6}' > DX_RK_$n/L2.txt
    \cp -f ${jdd}_RK3_PP_T.son DX_RK_$n/${jdd}_PP_T.son
done

# Comparaison a la solution analytique : 
for sch in "EUL" "RK"
do
    echo "Post traitement pour $sch "
    : > cvgx_son_$sch.txt
    : > cvgx_son_T_$sch.txt
    : > cvgx_L2_$sch.txt
    for n in 8 16 32 64 128 256
    do
        for compo in "T" 
        do 
            fic=DX_${sch}_$n/${jdd}_PP_${compo}.son
            ficout=DX_${sch}_$n/${jdd}_${compo}_complet.son 
            val=`awk 'NR==5{print $2}' $fic`
            valx=`awk 'NR==2{print $4}' $fic`
            valy=`awk 'NR==2{print $6}' $fic`
            echo " Processing $compo $sch for $n at point $valx, $valy.  Tini=$val"
            awk '{Time=$1;x='$valx'-0.01*Time;y='$valy'-0.01*Time;
            Tana=(0.5*cos(x*2*Pi/0.006)*cos(y*2*Pi/0.006));
            print Time, Tana,  $2-Tana}' \
                <  $fic > $ficout
            awk 'END{print '$n', $1, ($3**2)**(0.5)}' < $ficout >> cvgx_son_${compo}_$sch.txt
        done
        awk 'END{print '$n', $0}' DX_${sch}_$n/L2.txt >> cvgx_L2_$sch.txt
    done
done

##################
# Debut du post
##################
# cd QUICK
echo "Début post-traitement"
echo "---------------------"
: > plot.gplot
cat >> plot.gplot << EOF
# #!/usr/bin/gnuplot
set terminal png size 640,480 enhanced font "Helvetica,12"

# set output '../cvgx_son.png'
# set log xy
# set title "abs(err) finale pour la sonde P fonction de NX"
# plot   "../QUICK/cvgx_son_EUL.txt" u 1:3 t "Quick - EUL", \
#     "../CENTRE4/cvgx_son_EUL.txt" u 1:3 t "Centre4 - EUL", \
#     "../QUICK/cvgx_son_RK.txt" u 1:3 t "Quick - RK", \
#     "../CENTRE4/cvgx_son_RK.txt" u 1:3 t "Centre4 - RK", \
#     "../QUICK/cvgx_son_VX_EUL.txt" u 1:3 t "Quick - EUL", \
#     "../CENTRE4/cvgx_son_VX_EUL.txt" u 1:3 t "Centre4 - EUL", \
#     "../QUICK/cvgx_son_VX_RK.txt" u 1:3 t "Quick - RK", \
#     "../CENTRE4/cvgx_son_VX_RK.txt" u 1:3 t "Centre4 - RK", \
#     "../QUICK/cvgx_son_VY_EUL.txt" u 1:3 t "Quick - EUL - VY", \
#     "../CENTRE4/cvgx_son_VY_EUL.txt" u 1:3 t "Centre4 - EUL - VY", \
#     "../QUICK/cvgx_son_VY_RK.txt" u 1:3 t "Quick - RK - VY", \
#     "../CENTRE4/cvgx_son_VY_RK.txt" u 1:3 t "Centre4 - RK - VY", \
#     2.6e-9*(x/128)**(-2) w l t 'o2', \
#     1.e-8*(x/128)**(-4) w l t 'o4'

set output './cvgx_L2.png'
set log xy
set title "Norme L2 fonction de NX"
plot   "./cvgx_L2_EUL.txt" u 1:3 t "EUL", \
       "./cvgx_L2_RK.txt" u 1:3 t "RK", \
    0.0406266*(x/8)**(-1.5) w l ls 3 t 'o1.5', \
    0.0406266*(x/8)**(-2) w l ls 4 t 'o2', \
    0.0406266*(x/8)**(-3) w l ls 5 t 'o3'

EOF

export LC_ALL="en_US.UTF-8"
gnuplot plot.gplot
echo "You can display:", display $PWD/cvgx_L2.png
rm  *.lata*

#cd ..
