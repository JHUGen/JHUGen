#!/bin/bash

set -euo pipefail

(
cd $(dirname $0)

if [ -f NNPDF30_lo_as_0130.LHgrid ]; then  #easiest case: the file already exists
    exit 0
fi

if [ -L NNPDF30_lo_as_0130.LHgrid ]; then  #file is a broken link: get rid of it and then continue dealing with it
    rm NNPDF30_lo_as_0130.LHgrid
fi

if [ -f ../../JHUGenMELA/MELA/downloadNNPDF.sh ]; then  #if there is a MELA release in the standard place: use that
    ../../JHUGenMELA/MELA/downloadNNPDF.sh
fi

if [ -f ../../JHUGenMELA/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid ]; then  #...and link it
    ln -s ../../JHUGenMELA/MELA/data/Pdfdata/NNPDF30_lo_as_0130.LHgrid
    exit 0
fi

#no other choice but to download it
#originally from http://pcteserver.mi.infn.it/~nnpdf/nnpdf30/NNPDF30_lo_as_0130.LHgrid.tgz
#the link is now broken
wget http://spin.pha.jhu.edu/Generator/NNPDF30_lo_as_0130.LHgrid.tgz
tar -zxvf ./NNPDF30_lo_as_0130.LHgrid.tgz
rm ./NNPDF30_lo_as_0130.LHgrid.tgz

)
