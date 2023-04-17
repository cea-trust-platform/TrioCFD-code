#!/bin/bash

PATH=$project_directory/share/bin:$PWD/../../../../../../../bin:$PATH
# cd N80/P111_U11
cree_reprise
cd R
sed -e "s/# LECTURE DES INTERFACES DANS UN FICHIER LATA #/nom_reprise ijkft_stat_diph_gradUP_par8_init.sauv # /g" \
    -e 's/nproc_i 1/nproc_i 2/;s/nproc_j 1/nproc_j 2/;s/nproc_k 1/nproc_k 2/' \
    -e "s/# check_divergence//g" ijkft_stat_diph_gradUP.data > ijkft_stat_diph_gradUP_par8.data

cd ..
