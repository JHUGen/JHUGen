#!/bin/bash

LIBVERSION=2
LIBFOLDER=ggF-SMEFTsim-standalone_all
LIBNAME="libMG_SMEFTsim_v${LIBVERSION}.so"

bold=$(tput bold)
normal=$(tput sgr0)

if [ "$1" = "clean" ]; then
    find -name "*.o" -delete -o -name "*.a" -delete -o -name "*.so" -delete
    exit 0
fi

python3 madMela_mergeProcessing.py

echo "Done compiling!"
echo "Moving compiled library to the top directory..."
mv "libSMEFTSIM/libSMEFTsim.so" ${LIBNAME} || exit 1
if [[ -z "${MELA_LIB_PATH}" ]]; then
    echo "${bold}Please go to the MELA directory and run eval \$(./setup.sh env) to set your MELA paths"
    echo "Then place this library within the path specified in \$MELA_LIB_PATH ${normal}"
else
    echo "${bold}Copying compiled library to the path specified in \$MELA_LIB_PATH..."
    cp "${LIBNAME}" "${MELA_LIB_PATH}" || exit 1
    echo "Moved library to ${MELA_LIB_PATH} ${normal}"
fi
