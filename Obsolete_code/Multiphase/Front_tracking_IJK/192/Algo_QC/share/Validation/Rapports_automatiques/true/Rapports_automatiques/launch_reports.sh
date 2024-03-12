source ../../../env_New_algo_qc.sh

for file in *;
do
    if [ -d $file ];
    then
        cd $file
        Run_fiche
        cd -
    fi
done

echo "done"
