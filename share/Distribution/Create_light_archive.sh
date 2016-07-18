
file_liste=$PWD/share/Distribution/Liste_fiche 



#  pour trouver les grosses fiches

#make distrib
# find . -name '*'.lml.gz -size +100k | grep _jdd | awk -F/ '{print $NF}' | awk -F_jdd '{print "s/"$1" STD/"$1" NNR/"}' | sort -u
mkdir Prov
cd Prov
#tar zxvf ../TrioCFD.tar.gz
git clone `pwd`/.. TrioCFD
cd TrioCFD
baltik_build_configure
./configure



find . -name '*'.prm |grep src |awk -F/src '{print $1}' | sort  > L1
awk '{print $1}' $file_liste | sort > L2

if [ "`diff L1 L2`" != "" ]
then

diff L1 L2
echo Pb $file_list not uptodate
exit
fi


for fiche in `cat $file_liste | grep NNR | awk '{print $1}'`
do
  ff=`basename $fiche`
  for nr in `find . -name $ff'*.lml.gz'`
    do
    cat /dev/null > $nr
    touch $(dirname $nr)/skipped
    # echo  $(dirname $nr)/skipped
  done
done

rm -rf archives_dependancies
rm -f share/Distribution/Liste_fiche
mkdir share 
git log -1 > share/.gittag
make distrib
cp TrioCFD.tar.gz ../../TrioCFD_light.tar.gz
cd ../..
rm -rf Prov
