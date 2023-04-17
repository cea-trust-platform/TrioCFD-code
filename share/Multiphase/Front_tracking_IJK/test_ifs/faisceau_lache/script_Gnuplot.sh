#!/bin/bash


echo -n "plot 'z_cylindre_n.txt' w d " > script_1.gnuplot
echo -n ", 'v_cylindre_n.txt' w d " >> script_1.gnuplot; 
echo -n ", 'a_cylindre_n.txt' w d " >> script_1.gnuplot; 

echo -n "plot 'F_IBC_precedent.txt' w d " > script_2.gnuplot
echo -n ", 'F_IBC_n.txt' w d " >> script_2.gnuplot; 
echo -n ", 'F_IBC_extrapol.txt' w d " >> script_2.gnuplot; 

echo -n "plot 'F_inertie.txt' w d " > script_3.gnuplot
echo -n ", 'F_IBC_n_2.txt' w d " >> script_3.gnuplot; 
echo -n ", 'F_f_s.txt' w d " >> script_3.gnuplot; 


echo -n "plot 'z_cylindre_n_1.txt' w d " > script_4.gnuplot
echo -n ", 'v_cylindre_n_1.txt' w d " >> script_4.gnuplot; 
echo -n ", 'a_cylindre_n_1.txt' w d " >> script_4.gnuplot; 

echo -n "plot 'Cx.txt' w d " > script_5.gnuplot
echo -n ", 'Cz.txt' w d " >> script_5.gnuplot; 
echo -n ", 'Az_n.txt' w d " >> script_5.gnuplot; 
