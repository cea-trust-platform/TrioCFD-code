# source triangle.gdb

define printpoly
set $ptr=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_

set $i=((class Maillage_FT_IJK)mesh).convert_packed_to_ijk_cell(elem).data_[0]
set $j=((class Maillage_FT_IJK)mesh).convert_packed_to_ijk_cell(elem).data_[1]
set $k=((class Maillage_FT_IJK)mesh).convert_packed_to_ijk_cell(elem).data_[2]
set $x1=$ptr->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[0]
set $x2=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_->get_coords_of_dof($i+1,$j,$k, IJK_Splitting::NODES).v[0]
set $y1=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[1]
set $y2=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_->get_coords_of_dof($i,$j+1,$k, IJK_Splitting::NODES).v[1]
set $z1=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k, IJK_Splitting::NODES).v[2]
set $z2=((class Maillage_FT_IJK)mesh).ref_splitting_.pointeur_->get_coords_of_dof($i,$j,$k+1, IJK_Splitting::NODES).v[2]

printf "splot \"-\" w lp lt 1, \"-\" w lp lt 2, \"-\" w lp lt 3\n"

printf "%g %g %g\n", coordA.v[0], coordA.v[1], coordA.v[2]
printf "%g %g %g\ne\n", coordB.v[0], coordB.v[1], coordB.v[2]

set $som=0
while $som <= 2
 set $idx=facettes.operator(fa7, $som)
 set $s0=sommets.operator($idx, 0)
 set $s1=sommets.operator($idx, 1)
 set $s2=sommets.operator($idx, 2)
 printf "%g %g %g\n", $s0, $s1, $s2
 set $som = $som+1
end 

set $idx=facettes.operator(fa7, 0)
set $s0=sommets.operator($idx, 0)
set $s1=sommets.operator($idx, 1)
set $s2=sommets.operator($idx, 2)
printf "%g %g %g\ne\n", $s0, $s1, $s2

set $i=0
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

