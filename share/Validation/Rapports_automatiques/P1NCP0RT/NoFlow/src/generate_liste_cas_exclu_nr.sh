for REP in mesh_2 mesh_3
do
	   for f in $(find $REP -name "*.data" -type f)
	   do
	       echo $f >> liste_cas_exclu_nr
	   done
done
	   
