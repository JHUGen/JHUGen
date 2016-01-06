#!/bin/bash

INDIR="../../"


cd $INDIR/JHUGenerator
make || ( echo "compiling JHUGen failed!"; exit 1 )
make clean
cd - > /dev/null
cd $INDIR/JHUGenMELA
make || ( echo "compiling JHUGenMELA failed!"; exit 1 )
make clean
cd - > /dev/null



VERSION=$1
if [ ! $VERSION ]; then
    VERSION=$(grep JHUGen_Version ../../JHUGenerator/mod_Parameters.F90 | grep JHUGen_Version | sed -r 's/.*"(.*)".*/\1/')
fi

cp $INDIR"/Manual/manJHUGenerator.pdf" "../manJHUGenerator."$VERSION".pdf"
tar -czvf "JHUGenerator."$VERSION".tar.gz" $INDIR"/JHUGenerator" $INDIR"/JHUGenMELA" $INDIR"/AnalyticMELA" "../manJHUGenerator."$VERSION".pdf"
