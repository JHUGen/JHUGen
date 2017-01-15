#!/bin/bash

pdf="NNPDF30_lo_as_0130.LHgrid"
pdfdir="nnpdf30"
rootdir="../../"

if ! [ $(pwd | sed "s|.*/||") == pdfs ]; then
    echo "Please run downloadNNPDF.sh from the pdfs directory!"
    exit 1
fi

if [ ! -f $pdf ]; then

pdfsinside=$(find $rootdir -name $pdf)
npdfsinside=${#pdfsinside[@]}
if [ $npdfsinside -gt 0 ];then
    ln -s ${pdfsinside[0]} ./
else
    wget "http://pcteserver.mi.infn.it/~nnpdf/"$pdfdir"/"$pdf
    tar -zxvf "./"$pdf".tgz"
    rm "./"$pdf".tgz"
fi

fi
