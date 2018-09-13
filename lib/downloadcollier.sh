#!/bin/bash

set -euo pipefail

compiler=gfortran  #have to set the SAME compiler in the JHUGen makefile
if [ $# -ge 1 ]; then
  compiler = $1
fi

cd $(dirname $0)
wget http://www.hepforge.org/archive/collier/collier-1.2.tar.gz
tar -xvzf collier-1.2.tar.gz
rm collier-1.2.tar.gz
cd COLLIER-1.2/build
cmake -DCMAKE_Fortran_COMPILER=$compiler ..
make
