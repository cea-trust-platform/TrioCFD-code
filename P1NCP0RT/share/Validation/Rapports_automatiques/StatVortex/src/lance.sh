Parameter="Re1 Re100"
Discretization="P0 P1 P0P1 P0P1Pa P0RT"
for Discr in $Discretization ; do
  for Re in $Parameter ; do
    for num in `seq 1 3`; do
        fold="$project_directory/share/Validation/Rapports_automatiques/StatVortex/src/mesh_"$num"/"$Discr"/"$Re
        cd $fold
        $exec StatVortex.data 1>out 2>err&
    done
  done
done

