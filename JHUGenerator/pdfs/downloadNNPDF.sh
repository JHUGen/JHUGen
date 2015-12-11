#!/bin/bash

if ! [ $(pwd | sed "s|.*/||") == pdfs ]; then
    echo "Please run downloadNNPDF from the pdfs directory!"
    exit 1
fi

if [ -f NNPDF30_lo_as_0130.LHgrid ]; then
    exit 0
fi

if [ -f ../../JHUGenMELA/pdfs/NNPDF30_lo_as_0130.LHgrid ]; then
    ln -s ../../JHUGenMELA/pdfs/NNPDF30_lo_as_0130.LHgrid .
else
    # wget http://nnpdf.hepforge.org/html/nnpdf23/PDFsets/NNPDF23_lo_as_0130.LHgrid.tgz
    wget http://pcteserver.mi.infn.it/~nnpdf/nnpdf30/NNPDF30_lo_as_0130.LHgrid.tgz

    # tar -zxvf ./NNPDF23_lo_as_0130.LHgrid.tgz
    tar -zxvf ./NNPDF30_lo_as_0130.LHgrid.tgz

    # rm ./NNPDF23_lo_as_0130.LHgrid.tgz
    rm ./NNPDF30_lo_as_0130.LHgrid.tgz
fi

