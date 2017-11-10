#!/bin/bash

{

set -euo pipefail

scriptdir=$(dirname $0)
curdir=$(pwd)

cd $scriptdir

pkgname="collier-1.2"
pkgdir="COLLIER-1.2"
tarname=$pkgname".tar.gz"
tarweb="https://www.hepforge.org/archive/collier/"$tarname
libname="libcollier.so"
tmpdir="colliertmp"

if [[ $# > 0 ]] && [[ "$1" == *"clean"* ]];then

  rm -f *.so
  rm -f *.o
  rm -f *.mod
  rm -f *.f
  rm -f *.F
  rm -f *.F90
  for f in $(ls ./);do
    if [ -d $f ];then
      rm -rf $f
    fi
  done

  rm -f ../data/*/$libname

else

  if ! [ -f "../data/$SCRAM_ARCH/$libname" ]; then
    wget $tarweb
    mkdir $tmpdir
    tar -xvzf $tarname -C $tmpdir
    rm $tarname
    mv $tmpdir"/"$pkgdir"/src/"* ./
    rm -rf $tmpdir

    make
    cp $libname "../data/$SCRAM_ARCH/$libname"
  fi

fi

cd $curdir

}
