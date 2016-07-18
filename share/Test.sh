
ORG=`pwd`

#baltik_build_configure
#./configure
#time make debug dist_clean
prjs=`find . -mindepth 2  -name project.cfg -exec dirname {} \;`
for dir in $prjs
#for dir in  LES  Zoom validation
do
cd $dir
baltik_build_configure
./configure
env PAR_F=0 time make check_all_optim dist_clean
#time make optim dist_clean
cd $ORG
done
