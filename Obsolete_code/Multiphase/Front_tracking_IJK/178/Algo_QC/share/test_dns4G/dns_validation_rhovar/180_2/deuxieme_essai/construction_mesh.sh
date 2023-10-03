awk '{print $1-x;x=$1}' coord_k_level1.txt >taille_de_maille_niveau1.txt


seq 0 128 | awk '{printf "%16.15g\n",x;size=exp(-(($1*288.0/128-288)**2*2e-5))/81.6325441756531*0.0298460000000000/2;x=x+size}' >coord_k_level1_upper.txt
tac coord_k_level1_upper.txt  | awk 'NR>1{printf "%16.15g\n", 0.029846-$1}' >coord_k_level1_lower.txt
cat coord_k_level1_upper.txt coord_k_level1_lower.txt >coord_k_level1.txt

seq 0 96 | awk '{printf "%16.15g\n",x;size=exp(-(($1-96)**2*9e-5))/74.660491070367*0.0298460000000000/2;x=x+size}' >coord_k_level2_upper.txt
tac coord_k_level2_upper.txt  | awk 'NR>1{printf "%16.15g\n", 0.029846-$1}' >coord_k_level2_lower.txt
cat coord_k_level2_upper.txt coord_k_level2_lower.txt >coord_k_level2.txt
awk '{print $1-x;x=$1}' coord_k_level2.txt >taille_de_maille_niveau2.txt


seq 0 128 | awk '{printf "%16.15g\n",$1*0.029846/128}' >coord_k_level3.txt