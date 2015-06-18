#!/bin/bash
#Postraitement pour expension brusque "
tps_init_moyen=0
LC_NUMERIC=C
#---------------------------------------------------------------------------
#----- this for accumulating data for prev and this run
#----- strategy is read probres*.dat and probres*_*.dat
#----- and mix them into rescombine

Readn=1.0
t1=`cat proberes1_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes1.dat| awk '{print $2}' | head -1`
te=`cat proberes1.dat| awk '{print $3}' | head -1`
error1=`cat proberes1_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes1.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom1.dat


Readn=2.0
t1=`cat proberes2_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes2.dat| awk '{print $2}' | head -1`
te=`cat proberes2.dat| awk '{print $3}' | head -1`
error1=`cat proberes2_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes2.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom2.dat



Readn=3.0
t1=`cat proberes3_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes3.dat| awk '{print $2}' | head -1`
te=`cat proberes3.dat| awk '{print $3}' | head -1`
error1=`cat proberes3_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes3.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom3.dat


Readn=4.0
t1=`cat proberes4_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes4.dat| awk '{print $2}' | head -1`
te=`cat proberes4.dat| awk '{print $3}' | head -1`
error1=`cat proberes4_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes4.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom4.dat



Readn=5.0
t1=`cat proberes5_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes5.dat| awk '{print $2}' | head -1`
te=`cat proberes5.dat| awk '{print $3}' | head -1`
error1=`cat proberes5_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes5.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom5.dat



Readn=6.0
t1=`cat proberes6_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes6.dat| awk '{print $2}' | head -1`
te=`cat proberes6.dat| awk '{print $3}' | head -1`
error1=`cat proberes6_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes6.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom6.dat


Readn=7.0
t1=`cat proberes7_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes7.dat| awk '{print $2}' | head -1`
te=`cat proberes7.dat| awk '{print $3}' | head -1`
error1=`cat proberes7_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes7.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom7.dat


Readn=8.0
t1=`cat proberes8_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes8.dat| awk '{print $2}' | head -1`
te=`cat proberes8.dat| awk '{print $3}' | head -1`
error1=`cat proberes8_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes8.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom8.dat


Readn=9.0
t1=`cat proberes9_1.dat| awk '{print $2}' | head -1`
t2=`cat proberes9.dat| awk '{print $2}' | head -1`
te=`cat proberes9.dat| awk '{print $3}' | head -1`
error1=`cat proberes9_1.dat| awk '{print $4}' | head -1`
error2=`cat proberes9.dat| awk '{print $4}' | head -1`

awk -v x1=$Readn -v x2=$t1 -v x3=$t2 -v x4=$te -v x5=$error1 -v x6=$error2 'BEGIN {printf ("%.1f %.1f %.1f %.1f %.1f %.1f\n",x1,x2,x3,x4,x5,x6)}'  > recom9.dat

exit
