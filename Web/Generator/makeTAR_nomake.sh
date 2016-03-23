#!/bin/bash

INDIR="../../"


cd $INDIR/JHUGenerator
cd - > /dev/null
cd $INDIR/JHUGenMELA
cd - > /dev/null



VERSION=$1
if [ ! $VERSION ]; then
    VERSION=$(grep JHUGen_Version ../../JHUGenerator/mod_Parameters.F90 | grep JHUGen_Version | sed -r 's/.*"(.*)".*/\1/')
fi

pushd $INDIR"/Manual/"
# Do this three times
pdflatex manJHUGenerator.tex
pdflatex manJHUGenerator.tex
pdflatex manJHUGenerator.tex
popd
cp $INDIR"/Manual/manJHUGenerator.pdf" "../manJHUGenerator."$VERSION".pdf"
tar -czvf "JHUGenerator."$VERSION".tar.gz" $INDIR"/JHUGenerator" $INDIR"/JHUGenMELA" $INDIR"/AnalyticMELA" "../manJHUGenerator."$VERSION".pdf"
