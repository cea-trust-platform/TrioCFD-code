# Extraction des propiétés physiques, de la vitesse de frottement, et du flux de chaleur à la paroi
#--------------------------------------------------------------------------------------------------

mu=$(grep mu Cas.data | awk '{printf "%.10f",$4}')
rho=$(grep rho Cas.data | awk '{print $4}')
Cp=$(grep Cp Cas.data | awk '{print $4}')
u_tau=$(grep "u\*=" tble_post_vitesse_N=_NP_.dat | awk '{print $10}' | sed -n '$'p)
qw=$(cat Cas_pb_Diffusion_chaleur.out | awk '{print -$2/2}' | sed -n '$'p)


# Adimensionnement des sondes
#----------------------------

sed '/^$/d' Cas_PROFIL_VITESSE1.coupe > tmp1 # Elimination de la première ligne (vide) dans les fichiers .coupe
sed '/^$/d' Cas_PROFIL_TEMPERATURE1.coupe > tmp2

sed -n 1,18p tmp1 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2/'$u_tau'}' > Coupe_vitesse_TBLE_N=_NP__adimensionnee.dat # Adimensionnement (on ne garde que les lignes correspondant à la première moitié du canal)
sed -n 1,18p tmp2 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2*'$rho'*'$Cp'*'$u_tau'/'$qw'}' > Coupe_temperature_TBLE_N=_NP__adimensionnee.dat


# Adimensionnement des profils TBLE
#----------------------------------

num1=$(sed -n "/d_paroi/=" tble_post_vitesse_N=_NP_.dat | awk '{max_val=($1<max_val)?max_val:$1;} END{print max_val;}') # Maximum des numéros des lignes du fichier .dat où se trouve la chaîne "d_paroi"
num2=$(sed -n "/deq/=" T_tble_post_temperature_N=_NP_.dat | awk '{max_val=($1<max_val)?max_val:$1;} END{print max_val;}')

sed -n $(($num1+1)),'$'p tble_post_vitesse_N=_NP_.dat | awk '{print $1 "\t" $2}' > tmp3 # Affichage des lignes du fichier .dat comprises entre $num+1 et la dernière ligne, puis extraction des colonnes 1 et 2
sed -n $(($num2+1)),'$'p T_tble_post_temperature_N=_NP_.dat | awk '{print $1 "\t" $2}' > tmp4
sed '$d' tmp4 > tmp5

cat tmp3 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2/'$u_tau'}' > Profil_vitesse_TBLE_N=_NP__adimensionne.dat
cat tmp5 | awk '{print $1*'$u_tau'*'$rho'/'$mu' "\t" $2*'$rho'*'$Cp'*'$u_tau'/'$qw'}' > Profil_temperature_TBLE_N=_NP__adimensionne.dat

rm -f tmp*
