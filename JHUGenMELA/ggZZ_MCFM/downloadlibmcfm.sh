#!/bin/bash

if ! [ $(pwd | sed "s|.*/||") == ggZZ_MCFM ]; then
    echo "Please run downloadlibmcfm from the ggZZ_MCFM directory!"
    exit 1
fi

if [ -f libmcfm_7p0.so ]; then
    exit 0
fi

wget http://spin.pha.jhu.edu/otherdownloads/libmcfm_7p0.so
