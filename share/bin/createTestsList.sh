#!/usr/bin/bash

# usage : ./share/bin/createTestLis.sh
# This script updates the test and validation lists or each module.
# All lists are put in ./share/testList
# These lists can be used with make check_optim or make validation. For example on Turbulence module :
# make check_optim TESTLIST="./share/testList/testList_Turbulence"
# make check_validation TESTLIST=""./share/testList/validatioList_Turbulence"

REPLIST="./share/testList"

MODULES+="Algo_QC "
MODULES+="Turbulence "
MODULES+="Rayonnement/Rayonnement_milieu_transparent "
MODULES+="Rayonnement/Rayonnement_semi_transp "
MODULES+="Schema_Euler_Implicite_Stationnaire "
MODULES+="Zoom "
MODULES+="Fluid_Structure_Interaction "
MODULES+="Optimisation/Aposteriori "
MODULES+="Sensitivity_analysis "
MODULES+="Critere_Entrainement_Gaz "
MODULES+="P1NCP0RT "
MODULES+="Multiphase/Phase_field "
MODULES+="Multiphase/CMFD "
MODULES+="Multiphase/Front_tracking_discontinu "
MODULES+="Multiphase/Front_tracking_IJK "
MODULES+="validation "

mkdir -p ${REPLIST}

for MOD in ${MODULES}
do
  # ----------------------------------------------------------------------------
  # Creation de la liste des tests de non-regression associes a un module donne
  # ----------------------------------------------------------------------------
  TESTLIST="testList_$(basename ${MOD})"
  echo -e "Creation de la liste ${TESTLIST}"
  
  [[ -f ${TESTLIST} ]] && rm ${TESTLIST}
  touch ./${TESTLIST}
  
  for f in $(find "./tests/Reference/${MOD}" -iname "*.data" | sort -duf)
  do
    if [[ "$(basename ${f%.*})" == "$(basename $(dirname ${f}))" &&
         ! -f "$(dirname ${f})/skipped" ]]
    then
      #echo "$f"
      echo "$(basename ${f%.*})" >> ${TESTLIST}
    fi
  done
  
  for f in $(find "./share/Validation/Rapports_automatiques" -iwholename "*$(basename ${MOD})*src" | sort -duf)
  do
    rapport=$(basename $(dirname $f))
    echo -e " Rapport de validation '${rapport}'"
  
    for g in $(find "./tests/Reference/Validation" -iname "${rapport}_jdd*.data" | sort -duf)
    do
      #echo -e ${g}
      echo -e "$(basename ${g%.*})" >> ${TESTLIST}
    done
  
  done

  mv ${TESTLIST} ${REPLIST}

  # ----------------------------------------------------------------------------
  # Creation de la liste des rapports de validation associes a un module donne
  # ----------------------------------------------------------------------------
  VALIDLIST="validationList_$(basename ${MOD})"
  
  if [[ -d "./share/Validation/Rapports_automatiques/${MOD}" ]]
  then
    echo -e "Creation de la liste de validation ${VALIDLIST}"
    
    [[ -f ${VALIDLIST} ]] && rm ${VALIDLIST}
    touch ./${VALIDLIST}
    
    for r in $(find "./share/Validation/Rapports_automatiques/${MOD}" -iname "*.prm"  -not -path "*/build/*" -o -iname "*.ipynb" -not -path "*/.ipynb_checkpoints*" | sort -duf)
    do
      echo -e "Rapport de validation : ${r}"
      REPVAL=$(dirname ${r})
   
      if [[ "$(basename ${REPVAL})" == "src" ]]
      then
        REPVAL=$(dirname ${REPVAL})
      fi
   
      if [[ ! -f "${REPVAL}/skip_prm" ]] && [[ ! -f "${REPVAL}/src/skip_prm" ]]
      then
        echo -e "Ajout du rapport de validation : ${REPVAL}"
        echo -e ${REPVAL} >> ${VALIDLIST}
      fi
    done

    mv ${VALIDLIST} ${REPLIST}
  fi

done
