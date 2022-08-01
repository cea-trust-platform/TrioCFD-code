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
    export val1=`echo $1 | extract_points $4`
    export val2=`echo $2 | extract_points $4`
    export n=$3
    export n
    seq 1 $n | awk 'BEGIN{a=ENVIRON["val1"];b=ENVIRON["val2"];n=ENVIRON["n"]}{printf "%20g\n",a+$1/n*(b-a)}'
}
function symetry {
    awk '{printf "%20g\n",$1;x[NR]=$1;n=n+1}END{for(i=1;i<n;i++){printf "%20g\n",2*x[n]-x[n-i]}}' $1
}
# this function takes as input a file with node coordinates 
# and generates a file with node indexes, positions and element size
function preprocess {
    input=$1
    output=$2
    awk '{printf "%4d %20g %20g\n", NR, $1, $1-a;a=$1}' $input >$output
}

input=coord_k_level0.txt
preprocess $input level0_data.txt

seq 1 2 69 | extract_points level0_data.txt >tmp
linear 69 101 21 level0_data.txt >>tmp
symetry tmp >coord_k_level1.txt

preprocess coord_k_level1.txt level1_data.txt
seq 1 2 23 | extract_points level1_data.txt >tmp
linear 23 36 10 level1_data.txt >> tmp
seq 37 56 | extract_points level1_data.txt >>tmp
symetry tmp >coord_k_level2.txt

preprocess coord_k_level2.txt level2_data.txt
seq 1 2 7 | extract_points level2_data.txt >tmp
linear 7 13 5 level2_data.txt >> tmp
seq 14 42 | extract_points level2_data.txt >>tmp
symetry tmp >coord_k_level3.txt

preprocess coord_k_level3.txt level3_data.txt
# premier essai: genere 72 mailles.
#echo 1 | extract_points level3_data.txt >tmp
#linear 1 4 2 level3_data.txt >> tmp
#seq 5 38 | extract_points level3_data.txt >>tmp
#symetry tmp >coord_k_level4.txt

#  on voudrait 64 mailles pour pouvoir deraffiner 3 fois et decouper sur 2 processeurs en y,
#  donc on simplifie:
echo 1 | extract_points level3_data.txt >tmp
linear 1 38 32 level3_data.txt >> tmp
symetry tmp >coord_k_level4.txt

preprocess coord_k_level4.txt level4_data.txt
exit

echo 1 3 5 7 | extract_points x >y
linear 7 17  >>y
seq 18 51 | extract_points x >>y
symetry y >coord_k_level1.txt

input=coord_k_level1.txt
awk '{printf "%4d %20g %20g %20g\n", NR, $1, $1-a;a=$1}' $input >x1

echo 1 3 5 | extract_points x1 >y1
linear 5 25 x1 >>y1
seq 26 45 | extract_points x1 >>y1
symetry y1 >coord_k_level2.txt

input=coord_k_level2.txt
awk '{printf "%4d %20g %20g %20g\n", NR, $1, $1-a, ($1-a)*3;a=$1}' $input >x2

echo 1 | extract_points x2 >y2
linear 1 18 x2 1 >>y2
seq 19 35 | extract_points x2 >>y2
symetry y2 >coord_k_level3.txt
