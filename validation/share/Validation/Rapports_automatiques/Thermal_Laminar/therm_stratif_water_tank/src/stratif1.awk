#!/bin/bash
#Postraitement pour expension brusque "
tps_init_moyen=0
LC_NUMERIC=C
#---------------------------------------------------------------------------
# reads the last line,second variable in that line and puts it to t1

Readn=1.0
t1=`cat $1_SONDE_TEMP1.son| awk '{print $2}' | tail -1`
te1=`cat expttemp.dat| awk '{print $1}' | head -1`
error=`awk "BEGIN { print ($t1-$te1)/$te1*100.0 }" `

awk -v x1=$Readn -v x2=$t1 -v x3=$te1 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes1.dat

cat proberes1.dat >> result.dat

#---------------------------------------------------------------------------
Readn=2.0
t2=`cat $1_SONDE_TEMP2.son| awk '{print $2}' | tail -1`
te2=`cat expttemp.dat| awk '{print $2}' | head -1`
error=`awk "BEGIN { print ($t2-$te2)/$te2*100.0 }" `

awk -v x1=$Readn -v x2=$t2 -v x3=$te2 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes2.dat

cat proberes2.dat >> result.dat
#---------------------------------------------------------------------------
# Reading and processing 3rd reading
Readn=3.0
t3=`cat $1_SONDE_TEMP3.son| awk '{print $2}' | tail -1`
te3=`cat expttemp.dat| awk '{print $3}' | head -1`
error=`awk "BEGIN { print ($t3-$te3)/$te3*100.0 }" `

awk -v x1=$Readn -v x2=$t3 -v x3=$te3 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes3.dat

cat proberes3.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 4th reading
Readn=4.0
t4=`cat $1_SONDE_TEMP4.son| awk '{print $2}' | tail -1`
te4=`cat expttemp.dat| awk '{print $4}' | head -1`
error=`awk "BEGIN { print ($t4-$te4)/$te4*100.0 }" `

awk -v x1=$Readn -v x2=$t4 -v x3=$te4 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes4.dat

cat proberes4.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 5th reading
Readn=5.0
t5=`cat $1_SONDE_TEMP5.son| awk '{print $2}' | tail -1`
te5=`cat expttemp.dat| awk '{print $5}' | head -1`
error=`awk "BEGIN { print ($t5-$te5)/$te5*100.0 }" `

awk -v x1=$Readn -v x2=$t5 -v x3=$te5 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes5.dat

cat proberes5.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 6th reading
Readn=6.0
t6=`cat $1_SONDE_TEMP6.son| awk '{print $2}' | tail -1`
te6=`cat expttemp.dat| awk '{print $6}' | head -1`
error=`awk "BEGIN { print ($t6-$te6)/$te6*100.0 }" `

awk -v x1=$Readn -v x2=$t6 -v x3=$te6 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes6.dat
cat proberes6.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 7th reading
Readn=7.0
t7=`cat $1_SONDE_TEMP7.son| awk '{print $2}' | tail -1`
te7=`cat expttemp.dat| awk '{print $7}' | head -1`
error=`awk "BEGIN { print ($t7-$te7)/$te7*100.0 }" `

awk -v x1=$Readn -v x2=$t7 -v x3=$te7 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes7.dat

cat proberes7.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 8th reading
Readn=8.0

t8=`cat $1_SONDE_TEMP8.son| awk '{print $2}' | tail -1`
te8=`cat expttemp.dat| awk '{print $8}' | head -1`
error=`awk "BEGIN { print ($t8-$te8)/$te8*100.0 }" `

awk -v x1=$Readn -v x2=$t8 -v x3=$te8 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes8.dat

cat proberes8.dat >> result.dat

#---------------------------------------------------------------------------
# Reading and processing 9th reading

Readn=9.0

t9=`cat $1_SONDE_TEMP9.son| awk '{print $2}' | tail -1`
te9=`cat expttemp.dat| awk '{print $9}' | head -1`
error=`awk "BEGIN { print ($t9-$te9)/$te9*100.0 }" `

awk -v x1=$Readn -v x2=$t9 -v x3=$te9 -v x4=$error 'BEGIN {printf ("%.1f %.1f %.1f %.1f\n",x1,x2,x3,x4)}'  > proberes9.dat

cat proberes9.dat >> result.dat

exit
