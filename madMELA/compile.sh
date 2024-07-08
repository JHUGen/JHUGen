#!/bin/bash

LIBFOLDER=ggF-SMEFTsim-standalone_all
LIBNAME=libSMEFTsim.so

bold=$(tput bold)
normal=$(tput sgr0)

cd "${LIBFOLDER}"/Source/MODEL || exit 1
make clean || exit 1
make || exit 1

cd ../DHELAS || exit 1
make clean || exit 1
make || exit 1

cd ../../SubProcesses || exit 1
make cpp || exit 1
echo "Done compiling!"
echo "Moving compiled library to the top directory..."
mv "${LIBNAME}" ../../ || exit 1
cd ../../ || exit 1
if [[ -z "${MELA_LIB_PATH}" ]]; then
    echo "${bold}Please go to the MELA directory and run eval \$(./setup.sh env) to set your MELA paths"
    echo "Then place this library within the path specified in \$MELA_LIB_PATH ${normal}"
else
    echo "${bold}Moving compiled library to the path specified in \$MELA_LIB_PATH..."
    mv "${LIBNAME}" "${MELA_LIB_PATH}" || exit 1
    echo "Moved library to ${MELA_LIB_PATH} ${normal}"
fi
