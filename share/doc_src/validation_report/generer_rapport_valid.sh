#!/bin/bash

function usage()
{
  echo -e "Usage : ${0} path_to_validation_directory"
  echo -e "This script allows to generate the TrioCFD validation report."
  echo -e "The first argument must be the path to the validation directory"
  echo -e "which contains an 'archives' directory".
  echo -e "For example : ${0} /volatile/projets/trust-trio/validation/v1.9.1"
}

function performModif()
{
  local theTest=$1
  local theType=$2
  local theFile="${LOCAL_DIR}/${theTest}/.tmp/fic.tex"
  local theCorps="${LOCAL_DIR}/${theTest}/.tmp/corps.tex"

  sed -i '45d' ${theFile}
  sed -i "s@rhead{}@rhead{${theType}}@g" ${theFile}
  sed -i '/documentclass/d' ${theFile}
  sed -i '/usepackage/d' ${theFile}
  sed -i '/setlength/d' ${theFile}
  sed -i '/renewcommand/d' ${theFile}
  sed -i '/pagestyle/d' ${theFile}
  sed -i '/makeindex/d' ${theFile}
  sed -i '/begin{document}/d' ${theFile}
  sed -i '/huge/d' ${theFile}
  sed -i '/end{document}/d' ${theFile}
  sed -i '/\makesavenoteenv{tabular}/d' ${theFile}
  sed -i "s@orig{..}@orig{./archives/${theTest}}@g" ${theFile}
  sed -i "s@input{corps}@input{\\\orig/.tmp/corps}@g" ${theFile}

  # Specific action for "Channel_T1_T2_incompressible"
  # The name is too long and goes out the LaTeX table => it is cut !
  sed -ri "s|^Incompressible/.*/Canal.*_(M[01bistetraedrse_\\\]*) &|\1 \&|gmi" ${theCorps}
}

# ------------------------------------------------------------------------------
# Begining of the script
# ------------------------------------------------------------------------------
ARCHIVES="$1/archives"
LOCAL_DIR="archives"

if [[ $1 == "" ]]
then
  echo -e "Error: the first argument is mandatory !"
  usage
  exit 1
fi

if [[ ! -d  "${ARCHIVES}" ]]
then
  echo -e "Error: archives directory does not exist !"
  exit 1
fi 

# Creation of the '${ARCHIVES}' directory.
# All necessary .tgz will be copied here and will be de-archived
echo -e "Info: archives directory is ${ARCHIVES}"

if [[ -d ${LOCAL_DIR} ]]
then
  rm -rf ${LOCAL_DIR}
fi
mkdir -p ${LOCAL_DIR}

echo -e "Info: local directory is ${LOCAL_DIR}"

# Definition of the test list
declare -A testList

# Laminar reports
testList["Poiseuille_flow_2D_VDF_VEF"]="LAMINAR FLOW"
testList["Lid_driven_cavity"]="LAMINAR FLOW"
testList["Cir_Cyl_Re100"]="LAMINAR FLOW"

# Thermal and laminar reports
testList["Convection_Vahl_Davis"]="THERMAL LAMINAR FLOW"
testList["Oscillating_Flow"]="THERMAL LAMINAR FLOW"

# Turbulent reports
testList["OBI_diffuser_VEF_k_eps"]="TURBULENT FLOW"
testList["Mixing_length_VEF_WF"]="TURBULENT FLOW"

# Thermal and turbulent report
testList["Channel_T1_T2_incompressible"]="THERMAL TURBULENT FLOW"

# Front tracking reports
testList["FTD_oscillating_bubble"]="FRONT-TRACKING"
testList["FTD_hanging_drop"]="FRONT-TRACKING"

# ALE reports
testList["DivaALE"]="ALE"
testList["TwoCylindersALE"]="ALE"
testList["TwoOscillatingCylindersALE"]="ALE"
testList["3dOscillatingBeamALE"]="ALE"

for t in ${!testList[@]}
do
  echo -e "Copy in ${LOCAL_DIR}: ${t}.tgz"
  cp ${ARCHIVES}/${t}.tgz ${LOCAL_DIR}
  tar -xzf ${LOCAL_DIR}/${t}.tgz
  rm -rf ${LOCAL_DIR}/${t}.tgz
  mv preserve ${LOCAL_DIR}/${t}
  performModif "${t}" "${testList[${t}]}"
done

# Generation of the pdf report
DOC="TrioCFD_validation_report"
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex
pdflatex ${DOC}.tex

# Installation of the final pdf report
#cp ${DOC}.pdf ../../doc

# Cleaning
for ext in aux idx lof log lot toc
do
  rm "${DOC}.${ext}"
done

rm -rf ${LOCAL_DIR}

exit 0
