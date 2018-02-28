#!/bin/bash

set -euo pipefail

cd $(dirname $0)
wget http://www.hepforge.org/archive/collier/collier-1.1.tar.gz
tar -xvzf collier-1.1.tar.gz
rm collier-1.1.tar.gz
cd COLLIER-1.1/build
cmake ..
make
