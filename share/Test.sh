
ORG=`pwd`

baltik_build_configure
./configure
time make debug dist_clean

for dir in ALE Front_tracking_discontinu LES Phase_field  Rayonnement Rayonnement_semi_transp UtilitairesAssemblages Zoom validation K_Eps_non_std
#for dir in  LES  UtilitairesAssemblages Zoom validation
do
cd $dir
baltik_build_configure
./configure
time make check_optim dist_clean
cd $ORG
done
