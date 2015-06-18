#!/bin/bash
export LC_NUMERIC="C"

jdd=$1

# fluide
X=`grep "LONGUEURS" ./"$jdd".data| awk '{print $2}' | head -1`
Y=`grep "LONGUEURS" ./"$jdd".data| awk '{print $3}' | head -1`
Zf=`grep "LONGUEURS" ./"$jdd".data| awk '{print $4}' | head -1`

# solide
Zs=`grep "longueurs" ./*.data| awk '{print $4}' | head -1`


# fluid physical properties 
rho=`grep "RHO Champ_Uniforme" $jdd.data| awk '{print $4}' | head -1`
mu=`grep "MU Champ_Uniforme" $jdd.data| awk '{print $4}' | head -1`
lambda_f=`grep "LAMBDA Champ_Uniforme" $jdd.data| awk '{print $4}' | head -1`
cp_f=`grep "CP Champ_Uniforme" $jdd.data| awk '{print $4}' | head -1`

echo $rho > rho_f.dat
echo $mu > mu_f.dat
echo $lambda_f > lambda_f.dat
echo $cp_f > cp_f.dat

#solid properties
rho_s=`grep "rho champ_uniforme" $jdd.data| awk '{print $4}' | head -1`
cp_s=`grep "cp champ_uniforme" $jdd.data| awk '{print $4}' | head -1`
lambda_s=`grep "lambda champ_uniforme" $jdd.data| awk '{print $4}' | head -1`

echo $rho_s > rho_s.dat
echo $lambda_s > lambda_s.dat
echo $cp_s > cp_s.dat

# mean velocity
U=`grep "Vitesse champ_uniforme" ./*.data| awk '{print $4}' | head -1`
echo $U > mean_velocity.dat

# volume thermal power
q=`grep "puissance_thermique" $jdd.data| awk '{print $6}' | head -1`
echo $q > puissance_volumique.dat


##############################
### geometrical dimensions ###
##############################

h=$Zs
h_sol=$h

h_liq=`awk "BEGIN { print 2*$Zf }"`

S=`awk "BEGIN { print $X*$Y }"`

H=`awk "BEGIN { print $h_sol+$h_liq }"`

echo $X $Y $h_sol h=$h_sol > solide_dim.dat
echo $X $Y $h_liq H=$H > fluide_dim.dat

# parametre alpha du schema de convection THERMIQUE
alpha_th=`grep "EF_STAB" ./"$jdd".data| awk '{print $6}' | head -1`


####################
### hydraulique ####
####################

# diametre hydraulique Dh=h_liq
Dh_hyd=`awk "BEGIN { print $h_liq }"`
Re_hyd=`awk "BEGIN { print 2*$U*$rho*$Dh_hyd/($mu) }"`

yplus=`tail -3 ./"$jdd"_pb1_Ustar.face | awk '{print $4}' | head -1`

# perte de charge Dean
lambda=`awk "BEGIN { print 0.08681/$Re_hyd^0.25 }"`
gradP=`awk "BEGIN { print 2*$lambda*$U*$U/(2*$Dh_hyd) }"`
gradP_calc=`tail -1 $jdd"_Pressure_Gradient_pb1_PerioX" | awk '{print $2}' | head -1`
diff_gradp=`awk "BEGIN { print -100*($gradP-$gradP_calc)/$gradP }"`

echo $Re_hyd $yplus $gradP $gradP_calc $diff_gradp > pressure.dat

# friction velocity
ustar=`awk "BEGIN { print (0.5*$lambda)^0.5*$U }"`
ustar_calc=`tail -3 ./"$jdd"_pb1_Ustar.face | awk '{print $6}' | head -1`
diff_ustar=`awk "BEGIN { print -100*($ustar-$ustar_calc)/$ustar }"`

echo $Re_hyd $yplus $ustar $ustar_calc $diff_ustar > ustar.dat

# profil log
extrait_coupe "$jdd".data SEG_VITESSE_FLUIDE 
i=2
echo "#" > profil_velocity.dat
while test $i -lt 25
do
	yplus=`head -$i ./"$jdd"_SEG_VITESSE_FLUIDE.coupe | awk -v h=$h_sol -v ustar=$ustar -v rho=$rho -v mu=$mu '{print ($1-h)*ustar*rho/mu}'| tail -1 `
	Uplus_theo=`head -$i ./"$jdd"_SEG_VITESSE_FLUIDE.coupe | awk -v h=$h_sol -v ustar_calc=$ustar_calc -v rho=$rho -v mu=$mu '{print (1/0.41*(log(($1-h)*rho*ustar_calc/mu))+5.17)}'| tail -1 `
	Uplus_calc=`head -$i ./"$jdd"_SEG_VITESSE_FLUIDE.coupe | awk -v ustar=$ustar -v rho=$rho -v mu=$mu '{print $2/ustar}'| tail -1 `
	
i=$[$i+1]
echo $yplus $Uplus_theo $Uplus_calc >> profil_velocity.dat
done


##################
### thermique ####
##################

## theo ##
#diametre hydraulique Dh=2.h_liq
Dh_th=`awk "BEGIN { print 2*$h_liq }"`
Re_th=`awk "BEGIN { print $U*$rho*$Dh_th/($mu) }"`

slopeT=`awk "BEGIN { print 2*$q*$h_sol/($rho*$cp_f*$h_liq) }"`

# Prandtl
Pr=`awk "BEGIN { print $mu*$cp_f/$lambda_f }"`
echo $Pr > Prandtl.dat

# correlation Seban 
Nusselt=`awk "BEGIN { print 0.023*$Pr^0.3333*$Re_th^0.8 }"`
Tcontact=`awk "BEGIN { print $q*$h_sol*$Dh_th/($lambda_f*$Nusselt) }"`
Tmax=`awk "BEGIN { print $Tcontact+($q*$h_sol*$h_sol)/(2*$lambda_s) }"`

## calcul ##

# flux a la paroi de contact
flux_theo=`awk "BEGIN { print $q*$h_sol }"`
flux_calc=`tail -1 ./"$jdd"_pb1_Diffusion_chaleur.out | awk -v S=$S '{print $7/S}' | head -1`
diff_flux=`awk "BEGIN { print -100*($flux_theo-$flux_calc)/$flux_theo }"`

# T de melange
Tmel=`tail -1 ./Tmoyen_* | awk '{print $2}' | head -1`

# Tcontact
Tc_1=`tail -1 ./"$jdd"_TEMP_CONTACT_SOLIDE.son | awk '{print $2}' | head -1`
Tc_calc=`awk "BEGIN { print $Tc_1-$Tmel }"`
diff_tc=`awk "BEGIN { print -100*($Tcontact-$Tc_calc)/$Tcontact }"`

# T max
Tm_1=`tail -1 ./"$jdd"_TEMP_MAX.son | awk '{print $2}' | head -1`
Tmax_calc=`awk "BEGIN { print $Tm_1-$Tmel }"`
diff_tmax=`awk "BEGIN { print -100*($Tmax-$Tmax_calc)/$Tmax }"`

# pente temporelle
t1=`tail -1 ./"$jdd"_TEMP_CONTACT_SOLIDE.son | awk '{print $1}' | head -1`
T1=`tail -1 ./"$jdd"_TEMP_CONTACT_SOLIDE.son | awk '{print $2}' | head -1`
t2=`tail -4 ./"$jdd"_TEMP_CONTACT_SOLIDE.son | awk '{print $1}' | head -1`
T2=`tail -4 ./"$jdd"_TEMP_CONTACT_SOLIDE.son | awk '{print $2}' | head -1`
slopeT_calc=`awk "BEGIN { print ($T1-$T2)/($t1-$t2) }"`
diff_slope=`awk "BEGIN { print -100*($slopeT-$slopeT_calc)/$slopeT }"`

# deltaT
deltaT=`awk "BEGIN { print $q*$h_sol*$h_sol/(2*$lambda_s) }"`
deltaT_calc=`awk "BEGIN { print $Tmax_calc-$Tc_calc }"`
diff_delta=`awk "BEGIN { print -100*($deltaT-$deltaT_calc)/$deltaT }"`

# Nusselt calcule
Nu_calc=`awk "BEGIN { print $flux_calc*$Dh_th/($lambda_f*$Tc_calc) }"`
diff_nu=`awk "BEGIN { print -100*($Nusselt-$Nu_calc)/$Nusselt }"`

# results files
echo $Re_th $flux_theo $flux_calc $diff_flux > flux.dat

echo $Re_th $deltaT $deltaT_calc $diff_delta > delta.dat

echo $Re_th $slopeT $slopeT_calc $diff_slope > temp_slope.dat

echo $Re_th $Nusselt $Nu_calc $diff_nu > nusselt.dat

echo $Re_th $Tmax $Tmax_calc $diff_tmax > Tmax.dat

echo $Re_th $Tcontact $Tc_calc $diff_tc > Tc.dat



#### profils de temperature dans le solide
## profil THEORIQUE
echo 1 | awk -v q=$q -v lambda_s=$lambda_s -v lambda_f=$lambda_f -v h_sol=$h_sol -v Dh_th=$Dh_th -v Nusselt=$Nusselt '
	## temperature solide theorique
	function T_s(y)
		{
		toto=0.5*q/lambda_s*(h_sol*h_sol-y*y)+q*h_sol*Dh_th/(lambda_f*Nusselt);
		return toto;
		}
	{
	for(z=0;z<h_sol;z=z+0.05) {printf "%.4f %.4f\n",z,T_s(z) > "profil_temp_solide.dat"}
	}
	'

# profil CALCULE
extrait_coupe "$jdd".data SEG_TEMP_SOLIDE 
j=2
echo "#" > profil_calc_temp_solide.dat
while test $j -lt 25
do
	z=`head -$j ./"$jdd"_SEG_TEMP_SOLIDE.coupe | awk -v Tmel=$Tmel '{print ($1)}'| tail -1 `
	Ts=`head -$j ./"$jdd"_SEG_TEMP_SOLIDE.coupe | awk -v Tmel=$Tmel '{print ($2-Tmel)}'| tail -1 `

j=$[$j+1]
echo $z $Ts >> profil_calc_temp_solide.dat
done

exit





