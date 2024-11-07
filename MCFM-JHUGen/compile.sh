#!/bin/bash

(
set -euo pipefail

cd $(dirname $0)

JHUGenbasedir=../
if [[ -z "${SCRAM_ARCH+x}" ]];then
  GCCVERSION=$(gcc -dumpversion)
  if [[ "$GCCVERSION" == "4.3"* ]] || [[ "$GCCVERSION" == "4.4"* ]] || [[ "$GCCVERSION" == "4.5"* ]]; then # v1 of MCFM library
    export SCRAM_ARCH="slc5_amd64_gcc434"
  elif [[ "$GCCVERSION" == "4"* ]] || [[ "$GCCVERSION" == "5"* ]] || [[ "$GCCVERSION" == "6"* ]]; then # v2 of MCFM library
    export SCRAM_ARCH="slc6_amd64_gcc630"
  elif [[ "$GCCVERSION" == "7"* ]]; then # v3 of MCFM library
    export SCRAM_ARCH="slc7_amd64_gcc700"
  #elif [[ "$GCCVERSION" == "8"* ]]; then # v4 of MCFM library
  else
    export SCRAM_ARCH="slc7_amd64_gcc820"
  fi
fi
LIB=libmcfm_711.so

for aDir in QCDLoop/ff QCDLoop/ql QCDLoop TensorReduction/ov TensorReduction/pv TensorReduction/ov TensorReduction/recur/smallY TensorReduction/recur/smallP TensorReduction/recur/smallG TensorReduction/recur/smallF TensorReduction/recur
do
   pushd $aDir; mkdir -p obj; make "$@"; popd
done

mkdir -p obj
make "$@"
cd obj
targetdir="../"$JHUGenbasedir"/JHUGenMELA/MELA/data/"$SCRAM_ARCH
g++ -Wl,-soname,$LIB -shared -o $LIB *.o ../QCDLoop/ql/obj/ql*.o ../QCDLoop/ql/obj/a*.o ../QCDLoop/ff/obj/*.o ../TensorReduction/ov/*.a ../TensorReduction/pv/*.a ../TensorReduction/recur/*.a
mkdir -p $targetdir
cp $LIB $targetdir"/"
cd ..
)
