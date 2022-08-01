#!/bin/bash
cd src
DIR=$PWD
cd IJK
$DIR/gen.sh
cd solveur_mg
$DIR/gen.sh
cd ../OpVDF
$DIR/gen.sh
cd ../..
cd ..

make $1
