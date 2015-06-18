# Extraction des propiétés physiques, de la vitesse de frottement, et du flux de chaleur à la paroi
#--------------------------------------------------------------------------------------------------

mu=$(grep mu Cas.data | awk '{printf "%.10f",$4}')
rho=$(grep rho Cas.data | awk '{print $4}')
Cp=$(grep Cp Cas.data | awk '{print $4}')
u_tau=$(grep "<u\*>=" Cas_pb_Ustar.face | awk '{print $6}' | sed -n '$'p)
qw=$(cat Cas_pb_Diffusion_chaleur.out | awk '{print -$2/2}' | sed -n '$'p)


# Adimensionnement des sondes
#----------------------------

sed '/^$/d' Cas_PROFIL_VITESSE1.coupe > tmp1 # Elimination de la première ligne (vide) dans les fichier .coupe
sed '/^$/d' Cas_PROFIL_TEMPERATURE1.coupe > tmp2

sed -n 1,18p tmp1 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2/'$u_tau'}' > Coupe_vitesse_LP_standard__SCHEMADIFF__adimensionnee.dat # Adimensionnement (on ne garde que les lignes correspondant à la première moitié du canal)
sed -n 1,18p tmp2 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2*'$rho'*'$Cp'*'$u_tau'/'$qw'}' > Coupe_temperature_LP_standard__SCHEMADIFF__adimensionnee.dat

rm -f tmp*
