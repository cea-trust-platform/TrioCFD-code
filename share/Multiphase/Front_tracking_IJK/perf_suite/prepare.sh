#!/bin/bash
export prepro="python $TRUST_ROOT/bin/KSH/preprocessor.py"
export testcase_src_dir=$project_directory/share/perf_suite
export template_basename=mg_perf_test.template
export template_jdd=$testcase_src_dir/$template_basename

if [ "$project_directory" = "" ]
then
    echo "Baltik env not initialized"
    exit
fi

# Usage: lauch_case N resultfile
# Process template file tmp.template and runs the test case on N processors
# Extracts clock times for solvers except the two first and print into resultfile
# Total case execution time put on first line of resultfile, on first line with a comment
launch_case()
{
    resultdir=$2
    resultfile=${resultdir}.txt
    set $1
    solvtyp=$1
    ni=$2
    nj=$3
    nk=$4
    meshtype=$5
    nproc=$[$ni*$nj*$nk]

    tmpd=$resultdir
    mkdir $tmpd
    echo Running jdd on $nproc processors in $tmpd
ici=`pwd`
    cd $tmpd
    cp $template_jdd .

cat >tmp.template <<EOF
#Pinclude($template_basename)
#Pusemacro(DEF_MAIN)($solvtyp,$ni,$nj,$nk,$meshtype)
EOF
    $prepro tmp.template jdd.data
    if [ $nproc -eq 1 ]
	then
	datafile=jdd.data
    else
	make_PAR.data jdd.data
	datafile=PAR_jdd.data
    fi
    if [ "$BALTIKHOSTNAME" = "jade" ]
	then
	triou -create_sub_file PAR_jdd $nproc
	echo "$testcase_src_dir/extract_result.sh PAR_jdd.out >$resultfile"  >>sub_file
	qsub -lcluster=hpt sub_file
    elif [ "$BALTIKHOSTNAME" = "titane" ]
	then
	triou -create_sub_file PAR_jdd $nproc
	echo "$testcase_src_dir/extract_result.sh PAR_jdd.out >$resultfile"  >>sub_file
	ccc_msub -p gen6200 sub_file
    elif [ "$BALTIKHOSTNAME" = "occigen" ]
	then
	triou -create_sub_file PAR_jdd $nproc
	echo "$testcase_src_dir/extract_result.sh PAR_jdd.out >$resultfile"  >>sub_file
	sbatch sub_file
    else
	triou $datafile $nproc 1>PAR_jdd.out 2>PAR_jdd.err
	$testcase_src_dir/extract_result.sh PAR_jdd.out >$resultfile
    fi
    cd $ici
}

# testing gcp solver #
launch_case "GCP 1 1 1 MESH64x128x16" gcp_seq_64x128x16
launch_case "GCP 1 2 1 MESH64x128x16" gcp_par2_64x128x16
launch_case "GCP 2 2 1 MESH64x128x16" gcp_par4_64x128x16
launch_case "GCP 4 2 1 MESH64x128x16" gcp_par8_64x128x16

launch_case "MG 1 1 1 MESH64x128x16" mg_seq_64x128x16
launch_case "MG 1 2 1 MESH64x128x16" mg_par2_64x128x16
launch_case "MG 2 2 1 MESH64x128x16" mg_par4_64x128x16
launch_case "MG 4 2 1 MESH64x128x16" mg_par8_64x128x16

launch_case "MG 1 1 1 MESH128x128x128" mg_seq_128x128x128
launch_case "MG 1 2 1 MESH128x128x128" mg_par2_128x128x128
launch_case "MG 2 2 1 MESH128x128x128" mg_par4_128x128x128
launch_case "MG 2 2 2 MESH128x128x128" mg_par8_128x128x128

launch_case "GCP 1 1 1 MESH128x128x128" gcp_seq_128x128x128
launch_case "GCP 1 2 1 MESH128x128x128" gcp_par2_128x128x128
launch_case "GCP 2 2 1 MESH128x128x128" gcp_par4_128x128x128
launch_case "GCP 2 2 2 MESH128x128x128" gcp_par8_128x128x128
