#!/bin/bash

(
set -euo pipefail

cd $(dirname $0)

MELADIR="."
MCFMVERSION=mcfm_707

declare -i needSCRAM
needSCRAM=0
if [[ -z "${SCRAM_ARCH+x}" ]];then
  needSCRAM=1
  GCCVERSION=$(gcc -dumpversion)
  if [[ "$GCCVERSION" == "4.3"* ]] || [[ "$GCCVERSION" == "4.4"* ]] || [[ "$GCCVERSION" == "4.5"* ]]; then # v1 of MCFM library
    export SCRAM_ARCH="slc5_amd64_gcc434"
  elif [[ "$GCCVERSION" == "4"* ]] || [[ "$GCCVERSION" == "5"* ]] || [[ "$GCCVERSION" == "6"* ]]; then # v2 of MCFM library
    export SCRAM_ARCH="slc6_amd64_gcc630"
  elif [[ "$GCCVERSION" == "7"* ]]; then # v3 of MCFM library
    export SCRAM_ARCH="slc7_amd64_gcc700"
  #elif [[ "$GCCVERSION" == "8"* ]]; then # v4 of MCFM library
  else
    export SCRAM_ARCH="slc7_amd64_gcc820"
  fi
fi


printenv () {
  ldlibappend="$(readlink -f $MELADIR)/data/${SCRAM_ARCH}"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$ldlibappend"* ]];then
    echo "export LD_LIBRARY_PATH=${ldlibappend}${end}"
  fi

  pythonappend="$(readlink -f $MELADIR)/python"
  end=""
  if [[ ! -z "${PYTHONPATH+x}" ]]; then
    end=":${PYTHONPATH}"
  fi
  if [[ "${end}" != *"$pythonappend"* ]];then
    echo "export PYTHONPATH=${pythonappend}${end}"
  fi

  if [[ $needSCRAM -eq 1 ]];then
    echo "export SCRAM_ARCH=${SCRAM_ARCH}"
  fi
}
doenv () {
  ldlibappend="$(readlink -f $MELADIR)/data/${SCRAM_ARCH}"
  end=""
  if [[ ! -z "${LD_LIBRARY_PATH+x}" ]]; then
    end=":${LD_LIBRARY_PATH}"
  fi
  if [[ "${end}" != *"$ldlibappend"* ]];then
    export LD_LIBRARY_PATH="${ldlibappend}${end}"
    echo "Temporarily using LD_LIBRARY_PATH as ${LD_LIBRARY_PATH}"
  fi

  pythonappend="$(readlink -f $MELADIR)/python"
  end=""
  if [[ ! -z "${PYTHONPATH+x}" ]]; then
    end=":${PYTHONPATH}"
  fi
  if [[ "${end}" != *"$pythonappend"* ]];then
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
    echo 'eval $(./setup.sh env)'
    echo "or"
    echo 'eval `./setup.sh env`'
    echo
    echo "or do"
    echo './setup.sh env'
    echo "and change the commands according to your shell in order to do something equivalent to set up the environment variables."
    echo
else
    echo
    echo "ERROR: something went wrong in mv, see ^ error message"
    echo
    popd
    exit 1
fi
)
