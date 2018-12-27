#!/bin/bash

(
set -euo pipefail

cd $(dirname $0)

MELADIR="."
DATA_LIB_DIR="slc6_amd64_gcc530"

export SCRAM_ARCH=$DATA_LIB_DIR

printenv () {
    if [ -z "${LD_LIBRARY_PATH+x}" ]; then
      end=''
    else
      end=':$LD_LIBRARY_PATH'
    fi
    echo "export LD_LIBRARY_PATH=$(readlink -f $MELADIR)/data/$DATA_LIB_DIR$end"
    if [ -z "${PYTHONPATH+x}" ]; then
      end=''
    else
      end=':$PYTHONPATH'
    fi
    echo "export PYTHONPATH=$(readlink -f $MELADIR)/python$end"
}

if [[ "$#" -ge 1 ]] && [[ "$1" == "env" ]]; then
    printenv
elif [[ "$#" -ge 1 ]] && [[ "$1" == *"clean"* ]]; then
    COLLIER/setup.sh "$@"
    make clean
    pushd $MELADIR"/fortran/"
    make clean
    rm -f "../data/"$DATA_LIB_DIR"/libjhugenmela.so"
    popd
    make clean
else
    COLLIER/setup.sh "$@"
    tcsh data/retrieve.csh $DATA_LIB_DIR mcfm_705
    ./downloadNNPDF.sh
    pushd $MELADIR"/fortran/"
    make all
    if mv libjhugenmela.so "../data/"$DATA_LIB_DIR"/"; then
        echo
        echo "...and you are running setup.sh, so this was just done."
        echo
        popd
        make
        echo
        echo "remember to:"
        echo
        printenv
        echo
    else
        echo
        echo "ERROR: something went wrong in mv, see ^ error message"
        echo
        popd
        exit 1
    fi
fi
)
