#!/bin/bash

if [ -f FLAP/static/libflap.a ]; then
    echo "libflap.a found -- nothing to build."
else
    git clone https://github.com/szaghi/FLAP
    cd FLAP
    git submodule update --init
    FoBiS.py build -mode static-intel
    cd -
    test "$?BASH_VERSION" = "0" || eval 'setenv() { export "$1=$2"; }'
    setenv FLAP_DIR FLAP
fi
