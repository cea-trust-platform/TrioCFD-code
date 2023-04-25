Ub=0.01
nu=0.000001

extrait_coupe Cas VITESSE
extrait_coupe Cas TKE
extrait_coupe Cas EPSILON

for file in *coupe
do
	file_mod1=`echo "$file" | tr '[A-Z]' '[a-z]'`
	file_mod2=`echo "$file_mod1" | sed -e "s/cas_//g" | sed -e "s/.coupe/.dat/g"`
	mv "$file" "$file_mod2"
	sed -i 2d "$file_mod2"
done

rm -fr gnuplot *coupe
