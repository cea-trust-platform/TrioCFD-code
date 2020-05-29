#!/bin/bash
if [ ! -d triocfd-ihm ]
then
   git clone https://codev-tuleap.intra.cea.fr/plugins/git/triocfd/triocfd-ihm.git || exit -1
fi
cd triocfd-ihm

# Salome
[ "$SALOME_ROOT_DIR" = "" ] && echo "Le context de Salome 9.3.0 ou superieur doit etre charge." && exit -1
source $SALOME_ROOT_DIR/../../env_launch.sh # Car l'environnement Trust peut etre charge (omniOrb not found)

# Installation
[ ! -f _prefix/bin/triocfd.sh ] && echo "Installation..." && dev/install.sh

# Lancement
export TRIOCFD_ROOT_DIR=$(pwd)/_prefix
runSalome.py -mGEOM,SMESH,TRIOCFD
