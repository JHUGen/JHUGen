#!/bin/bash

set -euo pipefail

cd $(dirname $0)
collierversion=1.1

printenv () {
    if [[ "$OSTYPE" == "darwin"* ]]; then
      varname=DYLD_LIBRARY_PATH
    else
      varname=LD_LIBRARY_PATH
    fi
    echo "export $varname=$(readlink -f COLLIER-$collierversion):$""$varname"
}

if [ $# -gt 1 ]; then
  echo "unknown command line arguments"
  exit 1
elif [ $# -eq 1 ] && [ $1 == clean ]; then
  rm -r COLLIER-$collierversion
  exit
elif [ $# -eq 1 ] && [ $1 == env ]; then
  printenv
  exit
fi


compiler=gfortran  #have to set the SAME compiler in the JHUGen makefile
if [ $# -eq 1 ]; then
  compiler = $1
fi

wget http://www.hepforge.org/archive/collier/collier-$collierversion.tar.gz
tar -xvzf collier-$collierversion.tar.gz
rm collier-$collierversion.tar.gz
cd COLLIER-$collierversion/build
cmake -DCMAKE_Fortran_COMPILER=$compiler ..
make
echo
echo
echo "Remember to:"
echo
printenv
echo
