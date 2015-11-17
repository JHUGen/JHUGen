#!/bin/bash

INDIR="../../"


cd $INDIR/JHUGenerator
make || ( echo "compiling JHUGen failed!"; exit 1 )
cd - > /dev/null
cd $INDIR/JHUGenMELA
make || ( echo "compiling JHUGenMELA failed!"; exit 1 )
cd - > /dev/null



VERSION=$1
if [ ! $VERSION ]; then
    VERSION=$(grep JHUGen_Version ../../JHUGenerator/mod_Parameters.F90 | grep JHUGen_Version | sed -r 's/.*"(.*)".*/\1/')
fi

tar -czvf "JHUGenerator."$VERSION".tar.gz" $INDIR"/JHUGenerator" $INDIR"/JHUGenMELA" $INDIR"/AnalyticMELA"
