#!/bin/sh

cmsswdir=/work-zfs/lhc/usarica/CMS-related/MELAdevelopment/CMSSW_8_0_12/src
scramarch=slc6_amd64_gcc530

make
cd obj
g++ -Wl,-soname,libmcfm_701.so -shared -o libmcfm_701.so *.o ../QCDLoop/ql/obj/ql*.o ../QCDLoop/ql/obj/a*.o ../QCDLoop/ff/obj/*.o ../TensorReduction/ov/*.a ../TensorReduction/pv/*.a ../TensorReduction/recur/*.a
mkdir -p $cmsswdir"/ZZMatrixElement/MELA/data/"$scramarch
cp libmcfm_701.so $cmsswdir"/ZZMatrixElement/MELA/data/"$scramarch"/"
cd ..


