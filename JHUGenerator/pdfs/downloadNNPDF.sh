#!/bin/bash

pdf="NNPDF30_lo_as_0130.LHgrid"
pdfdir="nnpdf30"
rootdir="../../"

if ! [ $(pwd | sed "s|.*/||") == pdfs ]; then
    echo "Please run downloadNNPDF.sh from the pdfs directory!"
    exit 1
fi

if [ ! -f $pdf ]; then

pdfsinside="$(find $rootdir -name $pdf)"
let found=0
for p in ${pdfsinside[@]}
do
    echo $p
    if [[ "$p" == *"$pdf"* ]];then
        ln -s $p ./
        let found=1
        break
    fi
done
if [ $found -eq 0 ];then
    wget "http://pcteserver.mi.infn.it/~nnpdf/"$pdfdir"/"$pdf".tgz"
    tar -zxvf $pdf".tgz"
    rm $pdf".tgz"
fi

fi
