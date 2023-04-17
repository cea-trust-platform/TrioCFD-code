function extract_points {
    j=
    while read i
      do
      j="$j $i"
    done

    for k in $j
    do
      export k
      awk 'BEGIN{i=ENVIRON["k"]}$1==i{print $2}' $1
    done
}
function linear {
    export val1=`echo $1 | extract_points $3`
    export val2=`echo $2 | extract_points $3`
    val2m1=`echo $[$2-1] | extract_points $3`
    
    n=`echo $val1 $val2m1 $val2 | awk '{print int(($3-$1)/($3-$2))}'`
    if test -n "$4"
	then
	n=$[$n+$4]
    fi
    #echo val1=$val1 val2=$val2 val2m1=$val2m1 nvalues = $n
    export n
    seq 1 $n | awk 'BEGIN{a=ENVIRON["val1"];b=ENVIRON["val2"];n=ENVIRON["n"]}{printf "%20g\n",a+$1/n*(b-a)}'
}
function symetry {
    awk '{printf "%20g\n",$1;x[NR]=$1;n=n+1}END{for(i=1;i<n;i++){printf "%20g\n",2*x[n]-x[n-i]}}' $1
}

input=coord_k_level0.txt
awk '{printf "%4d %20g %20g %20g\n", NR, $1, $1-a, ($1-a)*3;a=$1}' $input >x
echo `seq 1 2 65` | extract_points x >y
symetry y >coord_k_level1.txt

input=coord_k_level1.txt
awk '{printf "%4d %20g %20g %20g\n", NR, $1, $1-a, ($1-a)*3;a=$1}' $input >x1

echo  `seq 1 2 33` | extract_points x1 >y1
symetry y1 >coord_k_level2.txt

input=coord_k_level2.txt
awk '{printf "%4d %20g %20g %20g\n", NR, $1, $1-a, ($1-a)*3;a=$1}' $input >x2
echo 1 | extract_points x2 >y2
linear 1 17 x2 2 >>y2
symetry y2 >coord_k_level3.txt
