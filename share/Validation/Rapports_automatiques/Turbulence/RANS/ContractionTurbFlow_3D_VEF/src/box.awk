#!/bin/bash
# sudden contraction 
# recirculation box

jdd=$1


# geometrical dimensions
D1=`grep "Diameter D1" $jdd.data| awk '{print $5}' | head -1`
L=`grep "Length L" $jdd.data| awk '{print $5}' | head -1`
	echo D1=$D1 m > info.txt
	echo Length=$L m >> info.txt

# physical properties 
rho=`grep "rho" $jdd.data| awk '{print $4}' | head -1`
mu=`grep "mu" $jdd.data| awk '{print $4}' | head -1`
	echo rho=$rho kg/m^3 >> info.txt
	echo mu=$mu kg/m/s >> info.txt


PI=3.1415926

#### correlation pressure loss

# inlet velocity
U=`grep "vitesse Champ_uniforme" $jdd.data| awk '{print $6}' | head -1`
	echo U=$U m/s >> info.txt

# total flow rate
Q=`awk "BEGIN { print $PI/4*$D1*$D1*$U }"`
	echo Qtotal=$Q m^3/s >> info.txt

# Reynolds numbers
Re=`awk "BEGIN { print $rho/$mu*$U*$D1 }"`
	echo Re=$Re >> info.txt

# regular pressure losses for turbulent flow in a pipe
F=`awk "BEGIN { print (0.316/$Re^0.25)*$L/$D1 }"`

PLR=`awk "BEGIN { print 0.5*$F*$rho*$U*$U }"`
	echo PLR=$PLR Pa >> info.txt



#### calculated pressure loss

# static pressure calculated
gradP=`tail -1 $jdd"_Pressure_Gradient_pb1_PERIO" | awk '{print $2}' | head -1`

PLcalc=`awk "BEGIN { print $gradP*($rho*$L) }"`
	echo PLcalc=$PLcalc Pa >> info.txt

#### difference between correlation and calculation
diff=`awk "
	function ABS(a) 
	{
		if(a<0) a=-a;
		return a;
	}
BEGIN { print 100*ABS(($PLcalc-$PLR)/$PLR) }"`
	echo diff_PL=$diff % >> info.txt

### filling of the results file res.dat
awk -v Re=$Re -v PLR=$PLR -v PLcalc=$PLcalc -v diff=$diff 'BEGIN { printf ("%.0f %.4f %.4f %.2f \n",Re,PLR,PLcalc,diff)}' > res_pl.dat


#### friction
# correlation tau 
tau_theo=`awk "BEGIN { print 0.0395*$rho*$U^(1.75)*($mu/$rho)^(0.25)*$D1^(-0.25) }"`
	echo tau_theo=$tau_theo >> info.txt

# calculated contrainte visqueuse
FX_calc=`tail -1 ./"$jdd"_pb1_Contrainte_visqueuse.out | awk '{print $2}' | head -1`

# surface of the wall
S=0.20943   

tau_calc=`awk "BEGIN { print $FX_calc/$S }"`
	echo tau_calc=$tau_calc >> info.txt

diff_2=`awk "
	function ABS(a) 
	{
		if(a<0) a=-a;
		return a;
	}
BEGIN { print 100*ABS(($tau_theo-$tau_calc)/$tau_theo) }"`
	echo diff_friction=$diff_2 % >> info.txt
	
	
### filling of the results file res.dat
awk -v Re=$Re -v tau_theo=$tau_theo -v tau_calc=$tau_calc -v diff_2=$diff_2 'BEGIN { printf ("%.0f %.4f %.4f %.2f \n",Re,tau_theo,tau_calc,diff_2)}' > res_friction.dat

### friction velocity
# theoritical friction velocity u* 
ustar_theo=`awk "BEGIN { print ((0.316/$Re^0.25)*$U*$U/8)^0.5 }"`
	echo ustar_theo=$ustar_theo m/s >> info.txt

# calculated y+
yplus_calc=`tail -3 ./"$jdd"_pb1_Ustar.face | awk '{print $4}' | head -1`

# calculated friction velocity
ustar_calc=`tail -3 ./"$jdd"_pb1_Ustar.face | awk '{print $6}' | head -1`
	echo ustar_calc=$ustar_calc m/s >> info.txt
	
diff_3=`awk "
	function ABS(a) 
	{
		if(a<0) a=-a;
		return a;
	}
BEGIN { print 100*ABS(($ustar_theo-$ustar_calc)/$ustar_theo) }"`
	echo diff_ustar=$diff_3 % >> info.txt

### filling of the results file res.dat
awk -v Re=$Re -v yplus=$yplus_calc -v ustar_theo=$ustar_theo -v ustar_calc=$ustar_calc -v diff_3=$diff_3 'BEGIN { printf ("%.0f %.0f %.4f %.4f %.2f \n",Re,yplus,ustar_theo,ustar_calc,diff_3)}' > res_ustar.dat

### preparation of the inlet profiles for the contration calculations
tmax=`grep "tmax" $jdd.data| awk '{print $2}' | head -1`

mv -f pb1_K_EPS_PERIO_"$tmax"*.dat pb1_K_EPS_PERIO.dat
mv -f pb1_VITESSE_PERIO_"$tmax"*.dat pb1_VITESSE_PERIO.dat

for dir in ../alpha_1 ../alpha_02
do
	cp -f pb1_K_EPS_PERIO.dat $dir
	cp -f pb1_VITESSE_PERIO.dat $dir
done

exit
