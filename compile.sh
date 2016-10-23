#!/bin/sh

cmsswdir=../
scramarch=slc6_amd64_gcc530

for aDir in QCDLoop/ff QCDLoop/ql QCDLoop TensorReduction/ov TensorReduction/pv TensorReduction/ov TensorReduction/recur/smallY TensorReduction/recur/smallP TensorReduction/recur/smallG TensorReduction/recur/smallF TensorReduction/recur
do
   pushd $aDir; mkdir -p obj; make; popd
done

mkdir -p obj
make
cd obj
g++ -Wl,-soname,libmcfm_701.so -shared -o libmcfm_701.so *.o ../QCDLoop/ql/obj/ql*.o ../QCDLoop/ql/obj/a*.o ../QCDLoop/ff/obj/*.o ../TensorReduction/ov/*.a ../TensorReduction/pv/*.a ../TensorReduction/recur/*.a
mkdir -p $cmsswdir"/ZZMatrixElement/MELA/data/"$scramarch
cp libmcfm_701.so $cmsswdir"/ZZMatrixElement/MELA/data/"$scramarch"/"
cd ..


