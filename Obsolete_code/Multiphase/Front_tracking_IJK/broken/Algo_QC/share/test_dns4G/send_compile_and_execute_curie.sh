# Send the source code and data file to curie
# Compile the source code
# Execute the data file

datafile=$1
echo datafile    = $datafile
destdatadir=/ccc/scratch/cont002/den/mathieub/dns_4g/$2
echo destdatadir = $destdatadir
basenamedatafile=`basename $1`

rootdir=/work/mathieu/tmp/New_algo_qc
rootsrc=$rootdir/src
destsrc=/ccc/work/cont002/den/mathieub/atelier_dns4g
calcdir=/ccc/scratch/cont002/den/mathieub/dns_4g

nproc=`awk '$1=="splitting"{print $3*$4*$5}' $datafile`
echo nproc = $nproc
ssh mathieub@curie-ccrt.ccc.cea.fr "mkdir $destdatadir"
scp -p $datafile mathieub@curie-ccrt.ccc.cea.fr:$destdatadir

scp -p `find $rootsrc -type f`  mathieub@curie-ccrt.ccc.cea.fr:$destsrc

ssh mathieub@curie-ccrt.ccc.cea.fr "
source /etc/profile.d/modules.sh ;
source /ccc/scratch/cont002/den/triou/Trio_U-1.6.5_patch/Trio_U/bin/Init_Trio_U ; 
export rep_dev=$destsrc; 
cd destsrc ; 
Makeatelier.sccs 2 ;
export exec=$destsrc/exec_opt/Trio_U_mpi_opt
cd $destdatadir ;
triou $basenamedatafile $nproc"

