#!/bin/bash

set -euo pipefail

cd $(dirname $0)
wget http://qcdloop.fnal.gov/QCDLoop-1.96.tar.gz
tar -xvzf QCDLoop-1.96.tar.gz
rm QCDLoop-1.96.tar.gz
cd QCDLoop-1.96 
make
