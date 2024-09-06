#!/bin/bash

jdd=ijkft_ConvectionMultiSphere_seq 
cas=ijkft_ConvectionMultiSphere_seq

echo "############### Cas : $cas  ###############"

# Dans cette partie on fait une étude en maillage à vitesse constante

rm -rf DX*

for n in 8 12 16 24 32 40 50 64 80 100 128 # 256
do
    mkdir -p DX_EUL_$n # DX_RK_$n
    cp ${jdd}.data DX_EUL_$n
    cp init* DX_EUL_$n
    cd DX_EUL_$n
    n_splitting=`echo "scale=0;0.035 * $n / 0.2 + 5" | bc`
    sed -e "s/nbelem_i .*/nbelem_i $n/g" \
        -e "s/nbelem_j .*/nbelem_j $n/g" \
        -e "s/nbelem_k .*/nbelem_k $n/g" \
        -e "s/ijk_splitting_ft_extension .*/ijk_splitting_ft_extension $n_splitting/g" -i ${jdd}.data
    echo -n "    Calculating DX_EUL_$n....."
    # trust ${jdd} 8 1> $jdd.out 2> $jdd.err
    echo "Done!"
    cd ..
done

# Dans cette partie on fait une étude en influence de la vitesse à maillage fixe

n=32
rm -rf V_*

for v in 1 5 8 10 20 40 50 100 # 40 50 64 80 100 128 256
do
mkdir -p V_$v # DX_RK_$n
cp ${jdd}.data V_$v
cp init* V_$v
cd V_$v
n_splitting=`echo "scale=0;0.035 * $n / 0.2 + 5" | bc`
v_fluid=`echo "scale=3;0.01 * $v" | bc`
sed -e "s/nbelem_i .*/nbelem_i $n/g" \
    -e "s/nbelem_j .*/nbelem_j $n/g" \
    -e "s/nbelem_k .*/nbelem_k $n/g" \
    -e "s/expression_vx_init .*/expression_vx_init $v_fluid/g" \
    -e "s/ijk_splitting_ft_extension .*/ijk_splitting_ft_extension $n_splitting/g" -i ${jdd}.data
    echo -n "    Calculating V_$v....."
    echo "Done!"
    cd ..
done
