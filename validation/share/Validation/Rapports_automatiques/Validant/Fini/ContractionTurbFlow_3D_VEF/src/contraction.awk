#!/bin/bash
# sudden contraction 
# subscript 1: upstream, subscript 2: downstream
# for each calculation the file info.txt gives the calculated values

export LC_NUMERIC="C"

jdd=$1


# geometrical dimensions
D1=`grep "Diameter D1" $jdd.data| awk '{print $5}' | head -1`
D2=`grep "Diameter D2" $jdd.data| awk '{print $5}' | head -1`
L1=`grep "Length L1" $jdd.data| awk '{print $5}' | head -1`
L2=`grep "Length L2" $jdd.data| awk '{print $5}' | head -1`
	echo D1=$D1 m > info.txt
	echo D2=$D2 m >> info.txt
	echo L1=$L1 m >> info.txt
	echo L2=$L2 m >> info.txt

# physical properties 
rho=`grep "rho" $jdd.data| awk '{print $4}' | head -1`
mu=`grep "mu" $jdd.data| awk '{print $4}' | head -1`
	echo rho=$rho kg/m3 >> info.txt
	echo mu=$mu kg/m/s >> info.txt

PI=3.1415926



#### correlation pressure loss

# inlet velocity
U1=`grep "Vitesse champ_uniforme" $jdd.data| awk '{print $4}' | head -1`
	echo U1=$U1 m/s >> info.txt

# flow rate
Q=`awk "BEGIN { print $PI/4*$D1*$D1*$U1 }"`
	echo Qtotal=$Q m^3/s >> info.txt

# downstream velocity
U2=`awk "BEGIN { print ($D1*$D1)/($D2*$D2)*$U1 }"`
	echo U2=$U2 m/s >> info.txt

# Reynolds numbers
Re1=`awk "BEGIN { print $rho/$mu*$U1*$D1 }"`
Re2=`awk "BEGIN { print $rho/$mu*$U2*$D2 }"`
	echo Re1=$Re1 >> info.txt
	echo Re2=$Re2 >> info.txt

# regular pressure losses for turbulent flow
f1=`awk "BEGIN { print 0.316/$Re1^0.25 }"`
f2=`awk "BEGIN { print 0.316/$Re2^0.25 }"`

F1=`awk "BEGIN { print $f1*$L1/$D1 }"`
F2=`awk "BEGIN { print $f2*$L2/$D2 }"`

PLR1=`awk "BEGIN { print 0.5*$F1*$rho*$U1*$U1 }"`
PLR2=`awk "BEGIN { print 0.5*$F2*$rho*$U2*$U2 }"`
	echo PLR_1=$PLR1 Pa >> info.txt
	echo PLR_2=$PLR2 Pa >> info.txt

# singular pressure loss
ksi=`awk "BEGIN { print 0.5*(1-(($D2*$D2)/($D1*$D1))) }"`
	echo ksi=$ksi >> info.txt

PLS=`awk "BEGIN { print 0.5*$ksi*$rho*$U2*$U2 }"` 
	echo PLS=$PLS Pa >> info.txt

# total pressure loss
PLtot=`awk "BEGIN { print $PLS+$PLR1+$PLR2 }"`
	echo PL_tot=$PLtot Pa >> info.txt

#### calculated pressure loss

# static pressure calculated
P1=`tail -1 ./sudden_contraction_PRESSION_AMONT.son | awk '{print $2}' | head -1`
P2=`tail -1 ./sudden_contraction_PRESSION_AVAL.son | awk '{print $2}' | head -1`

# pressure loss calculated
PLcalc=`awk "BEGIN { print $P1-($P2)+0.5*$rho*($U1*$U1-$U2*$U2) }"`
	echo PLcalc=$PLcalc Pa >> info.txt


#### difference between correlation and calculation
diff=`awk "
	function ABS(a) 
	{
		if(a<0) a=-a;
		return a;
	}
BEGIN { print 100*ABS(($PLcalc-$PLtot)/$PLtot) }"`
	echo diff=$diff % >> info.txt

#### filling of the results file res.dat
awk -v Re2=$Re2 -v PLtot=$PLtot -v PLcalc=$PLcalc -v diff=$diff 'BEGIN { printf ("%.0f %.2f %.2f %.1f \n",Re2,PLtot,PLcalc,diff)}' > res.dat


exit
