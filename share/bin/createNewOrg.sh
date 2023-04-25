#!/bin/bash
#set -xv

function exitOnError()
{
  local msg=${1}

  echo -e "${msg}"
  echo -e "Error in the script : exit now !"
  exit 1
}

function commentInProject()
{
  # Le module a commenter est passe en argument
  local MOD=$1
  for f in $(find ./ -iname "project.cfg")
  do
    #echo -e "Traitement du fichier ${f} pour ${MOD}"
    sed -ri "s/(.*${MOD} *$)/#\1/gmi" ${f}
  done
}

function findSymLink()
{
  local MOD=$1
  RES="laListeDesLiens.txt"
  
  echo -e "\n Lien pour le module ${MOD}" >> ${RES}
  find ${MOD} -type l >> ${RES} 
}

function testSymLink()
{
  for f in $(find ./ -type l)
  do
    ls -L ${f} > /dev/null 2>&1
    if [[ $? -ne 0 ]]
    then
      echo "Le lien semble mort : ${f}"
    fi
  done
}

function performChange()
{
  local MOD=$1
  echo -e "Change for module : ${MOD}"

  git status

  find ${MOD} -iname "*check_sources*" | xargs rm
  find ${MOD} -iname ".cproject" | xargs rm
  find ${MOD} -iname ".project" | xargs rm
  
  commentInProject "$(basename ${MOD})"
  
  git rm -f ${MOD}/project.cfg
  git rm ${MOD}/src/Version_kernel

  if [ -f "${MOD}/${DOC_SRC}/TRAD2_ajout0" ]
  then
    git rm -f "${MOD}/${DOC_SRC}/TRAD2_ajout0" 
  fi

  if [[ $(ls ./${MOD}/src/*) ]]
  then
    mkdir -p ./src/${MOD}
    git mv ${MOD}/src/* ./src/${MOD}/
    [[ $? -ne 0 ]] && exitOnError "Problem during mv ${MOD}/src"
    rmdir ${MOD}/src
  fi

  # Traitement des Rapports_automatiques
  if [[ -d ${MOD}/${RAPPORTS} ]]
  then
    mkdir -p "./${RAPPORTS}/${MOD}"
    for r in $(ls "./${MOD}/${RAPPORTS}")
    do
      echo -e "r vaut $r"
      git mv ./${MOD}/${RAPPORTS}/${r} ./${RAPPORTS}/${MOD}/
      [[ $? -ne 0 ]] && exitOnError "Problem during mv ${MOD}/${RAPPORTS}"
    done
    rmdir "${MOD}/${RAPPORTS}"
    rmdir "${MOD}/share/Validation"
  fi
  
  if [[ $(ls ${MOD}/share/doc/*) ]]
  then
    mkdir -p share/doc/${MOD}
    git mv ${MOD}/share/doc/* share/doc/${MOD}
    [[ $? -ne 0 ]] && exitOnError "Problem during mv share/doc"
  fi
  [[ -d ${MOD}/share/doc ]] && rmdir ${MOD}/share/doc

  if [[ $(ls ${MOD}/share/*) ]]
  then
    mkdir -p share/${MOD}
    git mv  ${MOD}/share/* share/${MOD}
    [[ $? -ne 0 ]] && exitOnError "Problem during mv share/*"
  fi
  [[ -d ${MOD}/share ]] && rmdir ${MOD}/share

  if [[ -d ${MOD}/${VALID} ]]
  then
    # Mise a jour des liens vers fiches de validation
    for f in $(find "./${MOD}/${VALID}/" -name "lien_fiche_validation")
    do
      #echo -e "Lien trouve : $f"
      sed -ri "s|(Rapports_automatiques)|\1/${MOD}|gmi" ${f}
    done

    git mv ${MOD}/${VALID}/* ./${VALID}/
    [[ $? -ne 0 ]] && exitOnError "Problem during mv ${MOD}/${VALID}/*"
    rmdir ${MOD}/${VALID}
  fi
  
  if [[ $(ls ${MOD}/${REFE}/*) ]]
  then
    mkdir -p ${REFE}/${MOD}
    git mv ${MOD}/${REFE}/* ${REFE}/${MOD}
    [[ $? -ne 0 ]] && exitOnError "Problem during mv ${MOD}/${REFE}/*"
  fi
  [[ -d ${MOD}/${REFE} ]] && rmdir ${MOD}/${REFE}

  if [[ $(ls ${MOD}/tests/*) ]]
  then
    mkdir -p ./tests/${MOD}
    git mv ${MOD}/tests/* ./tests/${MOD}
    [[ $? -ne 0 ]] && exitOnError "Problem during mv ${MOD}/tests/*"
  fi
  [[ -d "${MOD}/tests" ]] && rmdir "${MOD}/tests"

  rmdir ${MOD}
}


# ------------------------------------------------------------------------------
#
# Debut du script
#
# ------------------------------------------------------------------------------

RAPPORTS="share/Validation/Rapports_automatiques"
DOC="share/doc"
DOC_SRC="share/doc_src"
REFE="tests/Reference"
VALID="${REFE}/Validation"


# On fait un peu de menage pour retrouver un src "clean"
find ./src/* -type d | xargs rm -rf
find ./Multiphase/ -iname "build" | xargs -n 1 rm -rf
rm -rf ./share/archives
git checkout ./

# On cree un repertoire "tests" qui n'existait pas...
mkdir -p ./share/Validation
mkdir -p ./${REFE}
mkdir -p ./${VALID}

# A utiliser si on doit mettre a jour une branche qui touche le
go=0
if [[ ${go} -eq 1 ]]
then
  # Fichiers quasiment en double
  #meld ./Rayonnement/Rayonnement_milieu_transparent/src/Problemes_rayo.h ./Rayonnement/Rayonnement_semi_transp/src/patch/Problemes_rayo.h
  
  #meld ./Rayonnement/Rayonnement_milieu_transparent/src/Problemes_rayo.cpp ./Rayonnement/Rayonnement_semi_transp/src/patch/Problemes_rayo.cpp
  
  git rm -f ./Rayonnement/Rayonnement_milieu_transparent/src/Problemes_rayo.h
  git rm -f ./Rayonnement/Rayonnement_milieu_transparent/src/Problemes_rayo.cpp
  
  # Fichier en double :
  # - ./Turbulence/src/Specializations/VDF/Lois_Paroi/Hydr/Paroi_std_hyd_VDF.h
  # - ./Multiphase/Front_tracking_discontinu/src/patch/Paroi_std_hyd_VDF.h
  #rm -rf ./Multiphase/Front_tracking_discontinu/src/patch/Paroi_std_hyd_VDF.h
  #meld Turbulence/src/Specializations/VDF/Lois_Paroi/Hydr/Paroi_std_hyd_VDF.h ./Multiphase/Front_tracking_discontinu/src/patch/Paroi_std_hyd_VDF.h
  
  #cp ./Multiphase/Front_tracking_discontinu/src/patch/Paroi_std_hyd_VDF.h Turbulence/src/Specializations/VDF/Lois_Paroi/Hydr/Paroi_std_hyd_VDF.h
  
  git rm -f ./Multiphase/Front_tracking_discontinu/src/patch/Paroi_std_hyd_VDF.h
  
  # Fichier surcharge dans FT
  # - ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.h
  # - ./Turbulence/src/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.h
  # Il faut supprimer le fichier dans Turbulence
  #rm -rf ./src/Turbulence/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.h
  #rm -rf ./src/Turbulence/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.cpp

  #meld ./Turbulence/src/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.h ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.h
  #meld ./Turbulence/src/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.cpp ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.cpp

  cp ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.h ./Turbulence/src/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.h
  git rm -f ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.h
  
  cp ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.cpp ./Turbulence/src/ThHyd/Modeles_Turbulence/LES/Hydr/Mod_turb_hyd_ss_maille.cpp
  git rm -f ./Multiphase/Front_tracking_discontinu/src/patch/Mod_turb_hyd_ss_maille.cpp

  git add -u

  #git commit -m "Solve source files conflicts"

  exit 0
fi

go=0
if [[ ${go} -eq 1 ]]
then
  ./share/bin/post_configure
  if [[ $(git diff "./share/doc_src/TRAD2_ajout0") ]]
  then
    echo -e "Ecart dans TRAD2_ajout0"
    exit 1
  fi
  sed -ri "s|^(updateTRAD2_ajout0)$|# \1|gmi" ./share/bin/post_configure

  git add -u
  git commit -m "Stop updating 'TRAD2_ajout0'"
fi

go=1
if [[ ${go} -eq 1 ]]
then
  git cherry-pick --no-commit e2ceebfaee37f500c539eb2553d0f77755467e35
  git add -u
  git commit -m "Solve source files conflicts"
fi

# A decommenter a la demande
#MOD+="Turbulence "
#MOD+="Rayonnement/Rayonnement_milieu_transparent "
#MOD+="Rayonnement/Rayonnement_semi_transp "
#MOD+="Rayonnement "
#MOD+="Schema_Euler_Implicite_Stationnaire "
#MOD+="Zoom "
#MOD+="Fluid_Structure_Interaction "
#MOD+="Optimisation/Aposteriori "
#MOD+="Optimisation "
#MOD+="Sensitivity_analysis "
#MOD+="Critere_Entrainement_Gaz "
#MOD+="P1NCP0RT "
#MOD+="Multiphase/Phase_field "
#MOD+="Multiphase/CMFD "
#MOD+="Multiphase/Front_tracking_discontinu "
#MOD+="Multiphase/Front_tracking_IJK/IJK_Kernel "
#MOD+="Multiphase/Front_tracking_IJK "
#MOD+="Multiphase "
#MOD+="validation "

for rep in ${MOD}
do
  echo -e "Modifications pour le module ${rep}"
  
  #findSymLink "${rep}"
  if [[ "${rep}" == "Multiphase/CMFD" ]]
  then
    mkdir -p ./share/archives
    git mv ${rep}/archives/test_ftsatp.cpp ./share/archives
    rmdir ${rep}/archives
  fi

  if [[ "${rep}" == "Multiphase/Front_tracking_IJK/IJK_Kernel" ]]
  then
    mkdir -p ./share/archives
    git mv ${rep}/share/archives/* ./share/archives
    [[ $? -ne 0 ]] && exitOnError "Problem during mv share/archives"
    rmdir ${rep}/share/archives
    monTexte="[prerequisite1]\nname : fftw-3.3.8\nprog_test : test_fftw3.cpp\nlibrairies_flag : -lfftw3_mpi -lfftw3\nconfigure_flag : --enable-mpi --enable-shared\n"
    sed -ri "s/(^\[dependencies\]$)/${monTexte}\1/gmi" ./project.cfg
  fi
  
  if [[ "${rep}" == "validation" ]]
  then
    git rm -rf ${rep}/src
    for f in $(find ./${rep}/share/Validation/Rapports_automatiques -type l)
    do
      if [[ ! $(ls ${f}/ > /dev/null 2>&1) ]]
      then
        #echo "Le lien semble mort : ${f}"
        git rm ${f}
      fi
    done
    git reset HEAD "${rep}/share/Validation/Rapports_automatiques/Laminar/ConvergenceTaylorGreenWithDiffusion/src/"
    git checkout "${rep}/share/Validation/Rapports_automatiques/Laminar/ConvergenceTaylorGreenWithDiffusion/src/"
  fi

  performChange "${rep}"
  
  git add -u
  git commit -m "'$(basename ${rep})' is in TrioCFD"
  [[ $? -ne 0 ]] && exitOnError "Problem during the commit"

  if [[ "${rep}" == "Turbulence" ]]
  then
    git cherry-pick --no-commit 49103846f68f9765f1d65e3a277c49a895bf033a
    git add -u
    git commit -m "Fix test 'MachineLearning'"
  fi
  
  if [[ "${rep}" == "P1NCP0RT" ]]
  then
    git cherry-pick --no-commit 5f0bfb9ac6301ba853e314bd6d0220d12e29ca33
    git add -u
    git commit -m "Fix 'prepare' scripts for 'P1NCP0RT'"
  fi

done
    
# Reparation des liens symboliques
#git cherry-pick --no-commit 60715d6d4848da69b9c0c422b273901017f9ed0b
#git add -u
#git commit -m "Update symbolic links"


# Obsolete Code creation
#mkdir -p Obsolete_code
#git mv Multiphase Obsolete_code
#git mv Sensitivity_analysis Obsolete_code
#git add -u
#git commit -m "Old code is moved in 'Obsolete_code'"


# Creation des listes de tests (non-regression et validation)
#./scripts/createTestsList.sh
#git add ./share/testList
#git commit -m "Adding the list of non-regression and validation tests"

exit 0
