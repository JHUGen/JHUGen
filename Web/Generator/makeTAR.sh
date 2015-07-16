#!bin/bash/

INDIR="../../"
VERSION=$1

tar -czvf "JHUGenerator."$VERSION".tar.gz" $INDIR"/JHUGenerator" $INDIR"/JHUGenMELA" $INDIR"/AnalyticMELA"
