#!/bin/bash

(
set -euo pipefail

cd $(dirname $0)

JHUGenbasedir=../
if [[ -z "${SCRAM_ARCH+x}" ]];then
  export SCRAM_ARCH="slc6_amd64_gcc530"
fi
LIB=libmcfm_705.so

for aDir in QCDLoop/ff QCDLoop/ql QCDLoop TensorReduction/ov TensorReduction/pv TensorReduction/ov TensorReduction/recur/smallY TensorReduction/recur/smallP TensorReduction/recur/smallG TensorReduction/recur/smallF TensorReduction/recur
do
   pushd $aDir; mkdir -p obj; make; popd
done

mkdir -p obj
make
cd obj
targetdir="../"$JHUGenbasedir"/JHUGenMELA/MELA/data/"$SCRAM_ARCH
g++ -Wl,-soname,$LIB -shared -o $LIB *.o ../QCDLoop/ql/obj/ql*.o ../QCDLoop/ql/obj/a*.o ../QCDLoop/ff/obj/*.o ../TensorReduction/ov/*.a ../TensorReduction/pv/*.a ../TensorReduction/recur/*.a
mkdir -p $targetdir
cp $LIB $targetdir"/"
cd ..
)
