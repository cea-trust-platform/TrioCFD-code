#!/bin/bash
# source ~triou/env_triou_1.6.8.sh
# exec=/export/home/gb218285/BALTIKS/IJK/v1.6.8/DNS_IJK/New_algo_qc_opt
jdd=CONVECTION 
# tar xvzf init.tgz

cas=CONVECTION

# rm -rf $cas
# mkdir $cas
echo "############### Cas : $cas  ###############"
# cd $cas
# cp ../${jdd}.data .
#ln -s ../init.lata* .

######################
# Space........
######################
rm -rf DX*

declare -a MYARRY
MYARRY[8]=0.01200 
MYARRY[16]=0.01400 
MYARRY[32]=0.01300 
MYARRY[64]=0.01350 
MYARRY[128]=0.01325 
MYARRY[256]=0.013375 

for n in 8 16 32 64 128 256
do
    mkdir -p DX_EUL_$n # DX_RK_$n
    # Euler :
    sed -i -e "s/nbelem_i .*/nbelem_i $n/g" -e "s/nbelem_j .*/nbelem_j $n/g" ${jdd}.data
    echo -n "    Calculating DX_EUL_$n....."
    triou ${jdd} 1> out 2> err
    echo "Done!"
    grep "ERROR T FIELD" err | awk '{print $5, $6}' > DX_EUL_$n/L2.txt
    cp ${jdd}_PP_T.son DX_EUL_$n/
#    cp -f ${jdd}_P.son DX_EUL_$n/
#    cp -f ${jdd}_PP_VX.son DX_EUL_$n/
#    cp -f ${jdd}_PP_VY.son DX_EUL_$n/
    #
    # Sauvegarde du lata pour une figure a la fin: 
    if [ $n == 64 ]; then 
        mkdir -p DX_EUL_$n/lata
        cp -rf $jdd.data $jdd.lata* DX_EUL_$n/lata/
        echo "       (Copying lata to DX_EUL_$n for later use in plots...OK)"
    fi
    #
    # RK3 :
#     sed -e "/time_scheme/s/#//g" ${jdd}.data > ${jdd}_RK3.data 
#     echo -n "    Calculating DX_RK_$n....."
#     triou ${jdd}_RK3 1> out 2> err
#     echo "Done!"
#     grep "ERROR FIELD" err | awk '{print $4, $5, $6, $7}' > DX_RK_$n/L2.txt
#     \cp -f ${jdd}_RK3_P.son DX_RK_$n/${jdd}_P.son
#     \cp -f ${jdd}_RK3_PP_VX.son DX_RK_$n/${jdd}_PP_VX.son
#     \cp -f ${jdd}_RK3_PP_VY.son DX_RK_$n/${jdd}_PP_VY.son
done
#
#
# Comparaison a la solution analytique : 
for sch in "EUL"  # "RK"
do
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
            echo " Processing $compo $sch for $n at point $valx, $valy.  uini=$val"
           # if [ $compo == "VX" ] ; then fct="sin"; 
                # Je ne comprends pas ce que ce bloc fou la : 
                #      val=${MYARRY[$n]}
                #      echo "Processing $n at point $val."
                #      awk '{print $1, $2, 0.01*(1.+2.*sin(('$val'+0.01*$1)*98.17477042468103)), $2-0.01*(1.+2.*sin(('$val'+0.01*$1)*98.17477042468103))}' <  DX_${sch}_$n/${jdd}_P.son > DX_${sch}_$n/${jdd}_Pcomplet.son 
                #      awk 'END{print '$n', $1, ($4**2)**(0.5)}' < DX_${sch}_$n/${jdd}_Pcomplet.son >> cvgx_son_$sch.txt
                #	    
           # else fct="(-1.)+0.*";
           # fi
            awk '{x='$valx';y='$valy';T=$1;
            uana=0.01*(1.+2.*'$fct'(('$valy'+0.01*T)*98.17477042468103));
            print T, uana,  $2-uana}' \
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
: > plot.gplot
cat >> plot.gplot << EOF
#!/usr/bin/gnuplot
set terminal png large

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
    2.20918e-02*(x/8)**(-3) w l ls 3 t 'o3' \
#    2.39984e-06*(x/8)**(-3) w l ls 3 t 'o3' 
#    "../CENTRE4/cvgx_L2_EUL.txt" u 1:3 t "Centre4 - EUL", \
#    "../QUICK/cvgx_L2_RK.txt" u 1:3 t "Quick - RK", \
#    "../CENTRE4/cvgx_L2_RK.txt" u 1:3 t "Centre4 - RK", \
#    2.20918e-05*(x/8)**(-1) w l ls 1 t 'o1', \
#    2.39984e-06*(x/8)**(-1) w l ls 1 t 'o1', \
#    2.20918e-05*(x/8)**(-2) w l ls 2 t 'o2', \
#    2.39984e-06*(x/8)**(-2) w l ls 2 t 'o2', \
#    2.20918e-05*(x/8)**(-1.5) w l ls 4 t 'o1.5'

EOF

gnuplot plot.gplot

# cd ..
# display ./*png
rm  *.lata*

#cd ..
