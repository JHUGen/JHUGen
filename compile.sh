#!/bin/sh

cmsswdir=../
scramarch=$SCRAM_ARCH
LIB=libmcfm_704.so

for aDir in QCDLoop/ff QCDLoop/ql QCDLoop TensorReduction/ov TensorReduction/pv TensorReduction/ov TensorReduction/recur/smallY TensorReduction/recur/smallP TensorReduction/recur/smallG TensorReduction/recur/smallF TensorReduction/recur
do
   pushd $aDir; mkdir -p obj; make; popd
done

mkdir -p obj
make
cd obj
targetdir="../"$cmsswdir"/ZZMatrixElement/MELA/data/"$scramarch 
g++ -Wl,-soname,$LIB -shared -o $LIB *.o ../QCDLoop/ql/obj/ql*.o ../QCDLoop/ql/obj/a*.o ../QCDLoop/ff/obj/*.o ../TensorReduction/ov/*.a ../TensorReduction/pv/*.a ../TensorReduction/recur/*.a
mkdir -p $targetdir
cp $LIB $targetdir"/"
cd ..
