#!/bin/bash

(
set -euo pipefail

cd $(dirname $0)

MELADIR="."
MCFMVERSION=mcfm_705
if [[ -z "${SCRAM_ARCH+x}" ]];then
  export SCRAM_ARCH="slc6_amd64_gcc530"
fi

printenv () {
    if [ -z "${LD_LIBRARY_PATH+x}" ]; then
      end=''
    else
      end=':$LD_LIBRARY_PATH'
    fi
    echo "export LD_LIBRARY_PATH=$(readlink -f $MELADIR)/data/$SCRAM_ARCH$end"
    if [ -z "${PYTHONPATH+x}" ]; then
      end=''
    else
      end=':$PYTHONPATH'
    fi
    echo "export PYTHONPATH=$(readlink -f $MELADIR)/python$end"
}
doenv () {
  ldlibappend="$(readlink -f $MELADIR)/data/${SCRAM_ARCH}"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${LD_LIBRARY_PATH+x}" != *"$ldlibappend"* ]];then
    export LD_LIBRARY_PATH="${ldlibappend}${end}"
    echo "Temporarily using LD_LIBRARY_PATH as ${LD_LIBRARY_PATH}"
  fi

  pythonappend="$(readlink -f $MELADIR)/python"
  end=""
  if [[ ! -z "${PYTHONPATH+x}" ]]; then
    end=":${PYTHONPATH}"
  fi
  if [[ "${PYTHONPATH+x}" != *"$pythonappend"* ]];then
    export PYTHONPATH="${pythonappend}${end}"
    echo "Temporarily using PYTHONPATH as ${PYTHONPATH}"
  fi
}

if [[ "$#" -eq 1 ]] && [[ "$1" == "env" ]]; then
    printenv
    exit
elif [[ "$#" -eq 1 ]] && [[ "$1" == *"clean"* ]]; then
    #echo "Cleaning COLLIER"
    COLLIER/setup.sh "$@"

    pushd $MELADIR"/fortran/"
    #echo "Cleaning FORTRAN"
    make clean
    rm -f "../data/"$SCRAM_ARCH"/libjhugenmela.so"
    popd

    #echo "Cleaning C++"
    make clean

    exit
elif [[ "$#" -eq 1 ]] && [[ "$1" == *"-j"* ]]; then
    : ok
elif [[ "$#" -eq 0 ]]; then
    : ok
elif [[ "$#" -eq 2 ]] && [[ "$1" == *"-j"* ]]; then
    : ok
else
    echo "Unknown arguments:"
    echo "  $@"
    echo "Should be nothing, env, or clean"
    exit 1
fi

doenv
COLLIER/setup.sh "$@"
tcsh data/retrieve.csh $SCRAM_ARCH $MCFMVERSION
./downloadNNPDF.sh
pushd $MELADIR"/fortran/"
make all
if mv libjhugenmela.so "../data/"$SCRAM_ARCH"/"; then
    echo
    echo "...and you are running setup.sh, so this was just done."
    echo
    popd
    make "$@"
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
)
