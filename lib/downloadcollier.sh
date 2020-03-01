#!/bin/bash

set -euo pipefail

compiler=gfortran  #have to set the SAME compiler in the JHUGen makefile
declare -i ncores
ncores=0

if [[ "$#" -eq 1 ]] && [[ "$1" == *"clean"* ]]; then
  rm -rf COLLIER*
  rm -f *.tar.gz
  exit 0
elif [[ "$#" -eq 1 ]] && [[ "$1" == *"-j"* ]]; then
  ncores=-1
elif [[ "$#" -eq 3 ]] && [[ "$1" == *"-j"* ]] && [[ "$2" == *"compiler"* ]]; then
  ncores=-1
  compiler=$3
elif [[ "$#" -eq 3 ]] && [[ "$3" == *"-j"* ]] && [[ "$1" == *"compiler"* ]]; then
  ncores=-1
  compiler=$2
elif [[ "$#" -eq 0 ]]; then
  : ok
elif [[ "$#" -eq 2 ]] && [[ "$1" == *"-j"* ]]; then
  ncores=$2
elif [[ "$#" -eq 2 ]] && [[ "$1" == *"compiler"* ]]; then
  compiler=$2
elif [[ "$#" -eq 4 ]] && [[ "$1" == *"-j"* ]] && [[ "$3" == *"compiler"* ]]; then
  ncores=$2
  compiler=$4
elif [[ "$#" -eq 4 ]] && [[ "$3" == *"-j"* ]] && [[ "$1" == *"compiler"* ]]; then
  ncores=$4
  compiler=$2
fi


cd $(dirname $0)

wget http://www.hepforge.org/archive/collier/collier-1.2.0.tar.gz
tar -xvzf collier-1.2.0.tar.gz
rm collier-1.2.0.tar.gz
mv COLLIER-1.2 COLLIER

cd COLLIER/build
cmake -DCMAKE_Fortran_COMPILER=$compiler ..
if [[ $ncores -eq 0 ]];then
  make
elif [[ $ncores -eq -1 ]];then
  make -j
else
  make -j $ncores
fi
