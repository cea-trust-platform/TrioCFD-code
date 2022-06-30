define printpoly
set $ptr=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_

set $i=((class Maillage_FT_IJK)maillage).convert_packed_to_ijk_cell(num_element).data_[0]
set $j=((class Maillage_FT_IJK)maillage).convert_packed_to_ijk_cell(num_element).data_[1]
set $k=((class Maillage_FT_IJK)maillage).convert_packed_to_ijk_cell(num_element).data_[2]
set $x1=$ptr->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[0]
set $x2=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_->get_coords_of_dof($i+1,$j,$k, IJK_Splitting::NODES).v[0]
set $y1=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[1]
set $y2=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_->get_coords_of_dof($i,$j+1,$k, IJK_Splitting::NODES).v[1]
set $z1=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[2]
set $z2=((class Maillage_FT_IJK)maillage).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k+1, IJK_Splitting::NODES).v[2]

printf "splot \"-\" w lp lt 1, \"-\" w lp lt 2, \"-\" w lp lt 3\n"
set $n=poly_.dimensions_[0]
set $i=0
while $i <= $n
 set $j = $i % $n
 set $x=coord_som[0][0]*poly_.operator($j,0) + coord_som[1][0]*poly_.operator($j,1) + coord_som[2][0]*poly_.operator($j,2)
 set $y=coord_som[0][1]*poly_.operator($j,0) + coord_som[1][1]*poly_.operator($j,1) + coord_som[2][1]*poly_.operator($j,2)
 set $z=coord_som[0][2]*poly_.operator($j,0) + coord_som[1][2]*poly_.operator($j,1) + coord_som[2][2]*poly_.operator($j,2)
 printf "%g %g %g\n", $x, $y, $z
 set $i = $i+1
end

printf "e\n"
set $i=0
while $i <= 2
 printf "%g %g %g\n", coord_som[$i][0], coord_som[$i][1], coord_som[$i][2]
 set $i = $i+1
end
set $i=0
printf "%g %g %g\ne\n", coord_som[$i][0], coord_som[$i][1], coord_som[$i][2]
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y1,$z1,$x2,$y1,$z1
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y2,$z1,$x2,$y2,$z1
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y1,$z2,$x2,$y1,$z2
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y2,$z2,$x2,$y2,$z2

printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y1,$z1,$x1,$y2,$z1
printf "%g %g %g\n%g %g %g\n\n\n",$x2,$y1,$z1,$x2,$y2,$z1
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y1,$z2,$x1,$y2,$z2
printf "%g %g %g\n%g %g %g\n\n\n",$x2,$y1,$z2,$x2,$y2,$z2

printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y1,$z1,$x1,$y1,$z2
printf "%g %g %g\n%g %g %g\n\n\n",$x2,$y1,$z1,$x2,$y1,$z2
printf "%g %g %g\n%g %g %g\n\n\n",$x1,$y2,$z1,$x1,$y2,$z2
printf "%g %g %g\n%g %g %g\n\n\n",$x2,$y2,$z1,$x2,$y2,$z2

printf "e\n"
end

# exempl: liste_intersections_element intersections num_elem
# affiche la liste des faceettes qui coupent l'element num_elem
define liste_intersections_element
set $intersect=$arg0
set $elem=$arg1
set $idx=$intersect.index_elem_facette_.data_[$elem]
while ($idx >= 0)
print $intersect.data[$idx].numero_facette_
set $idx=$intersect.data[$idx].index_facette_suivante_
end
end

# suppose qu'on a la zone dans refzone_dis_,
# imprime l'element eulerien elem dans le fichier dump_elem.txt pour gnuplot
define dumpelem
print dump_elem(*(class Zone_VDF *)&refzone_dis_.valeur().valeur(), $arg0, 0)
end

# exemple: dump_intersections_element intersections num_elem
# ajoute au fichier dump_facette.txt toutes les facettes qui coupe l'element num_elem
define dump_intersections_element
set $intersect=$arg0
set $elem=$arg1
set $idx=$intersect.index_elem_facette_.data_[$elem]
while ($idx >= 0)
print $intersect.data[$idx].numero_facette_
print dump_facette(mesh, $intersect.data[$idx].numero_facette_, 1)
set $idx=$intersect.data[$idx].index_facette_suivante_
end
end
