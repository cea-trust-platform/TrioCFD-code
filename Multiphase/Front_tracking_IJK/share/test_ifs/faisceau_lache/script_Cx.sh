#!/bin/bash

######################calcul h

 awk '{if ($1=="'z_cylindre_n'") print c, $2;}' von_karman.out > z_cylindre_n.txt
 awk '{if ($3=="'v_cylindre_n'") print c, $4;}' von_karman.out > v_cylindre_n.txt
 awk '{if ($5=="'a_cylindre_n'") print c, $6;}' von_karman.out > a_cylindre_n.txt
 awk '{if ($7=="'Az_n'") print c, $8;}' von_karman.out > Az_n.txt

 awk '{if ($1=="'F_IBC_precedent'") print c, $2;}' von_karman.out > F_IBC_precedent.txt
 awk '{if ($3=="'F_IBC_n'") print c, $4;}' von_karman.out > F_IBC_n.txt
 awk '{if ($5=="'F_IBC_extrapol'") print c, $6;}' von_karman.out > F_IBC_extrapol.txt
 
 
 awk '{if ($1=="'F_inertie'") print c, $2;}' von_karman.out > F_inertie.txt
 awk '{if ($3=="'F_IBC_n_bis'") print c, $4;}' von_karman.out > F_IBC_n_2.txt
 awk '{if ($5=="'F_f_s'") print c, $6;}' von_karman.out > F_f_s.txt

 awk '{if ($1=="'z_cylindre_n_1'") print c, $2;}' von_karman.out > z_cylindre_n_1.txt
 awk '{if ($5=="'v_cylindre_n_1'") print c, $6;}' von_karman.out > v_cylindre_n_1.txt
 awk '{if ($7=="'a_cylindre_n_1'") print c, $8;}' von_karman.out > a_cylindre_n_1.txt

     awk '{if ($1=="'Cx'") print c, $2;}' von_karman.out > Cx.txt
     awk '{if ($3=="'Cz'") print c, $4;}' von_karman.out > Cz.txt


     awk '{if (($1=="'F_cylindre'")&&($2=="'12'")) print c, $6;}' von_karman.out > Fx_12.txt
     awk '{if (($1=="'F_cylindre'")&&($2=="'12'")) print c, $8;}' von_karman.out > Fz_12.txt
