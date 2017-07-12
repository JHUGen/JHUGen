#!/bin/sh

MELADIR="."
DATA_LIB_DIR="slc6_amd64_gcc530"

if [[ "$1" == *"clean"* ]];then
	make clean
	pushd $MELADIR"/fortran/"
	make clean
	rm -f "../data/"$DATA_LIB_DIR"/libjhugenmela.so"
	popd
	make clean
else
	tcsh data/retrieve.csh $DATA_LIB_DIR mcfm_703
	bash downloadNNPDF.sh
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
                echo "export LD_LIBRARY_PATH=$(readlink -f $MELADIR)/data/$DATA_LIB_DIR"':$LD_LIBRARY_PATH'
                echo
	else
		echo
		echo "ERROR: something went wrong in mv, see ^ error message"
		echo
		popd
		return 1 >& /dev/null || exit 1 #return only works when sourced, exit will exit your whole session if sourced
	fi
fi
