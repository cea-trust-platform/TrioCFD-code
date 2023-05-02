awk '$1=="Centre_gravite_phases"{print $4, $6,$7,$8}' FTD_reprise_xyz_vef_3d/prepare.err>FTD_reprise_xyz_vef_3d/position_goutte.txt 
echo >>FTD_reprise_xyz_vef_3d/position_goutte.txt
echo >>FTD_reprise_xyz_vef_3d/position_goutte.txt
awk '$1=="Centre_gravite_phases"{print $4, $6,$7,$8}' FTD_reprise_xyz_vef_3d/FTD_reprise_xyz_vef_3d.err >>FTD_reprise_xyz_vef_3d/position_goutte.txt 

