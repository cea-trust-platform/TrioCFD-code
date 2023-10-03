#!/bin/bash

# Extrait de deux moyennes statistiques  la valeur moyenne entre les deux instants

f1=statistiques_$1.txt
f2=statistiques_$2.txt
export t1=`head -n 1 $f1| awk '{print $3}'`
export t2=`head -n 1 $f2| awk '{print $3}'`
grep -v \# $f1 > $f1.tmp
grep -v \# $f2 > $f2.tmp
paste $f1.tmp $f2.tmp >pasted.tmp

export ncol=`head -n 1 $f1.tmp | wc -w`
echo t1=$t1
echo t2=$t2
echo ncol=$ncol

awk 'BEGIN{ncol=ENVIRON["ncol"];t1=ENVIRON["t1"];t2=ENVIRON["t2"]}
{ printf "%g ",$1;
for (i=2;i<=ncol;i++){printf"%g ",($(i+ncol)*t2-$(i)*t1)/(t2-t1)};printf"\n"}' pasted.tmp >diff_stat_${t1}_${t2}.txt

rm -f pasted.tmp $f1.tmp $f2.tmp
