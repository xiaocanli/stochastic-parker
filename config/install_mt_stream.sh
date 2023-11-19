#!/bin/bash

build_dir=$PWD

if [ -f $HOME/.local/mt_stream/mt_stream.o ]; then
    echo "mt_stream.o found -- nothing to build."
else
    mkdir $HOME/.local/mt_stream
    cd $HOME/.local/mt_stream
    wget http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_f90.tar.gz
    tar zxvf mt_stream_f90.tar.gz
    cd mt_stream_f90-1.11
    sed -i -e 's/ifort/gfortran/g' Makefile
    sed -i -e 's/icpc/g++/g' Makefile
    sed -i -e 's/FFLAGS := -O3 -g -traceback/FFLAGS := -O3 -g -fno-range-check/g' Makefile
    sed -i -e 's/-static//g' Makefile
    # fix compiling errors when using new versions of gfortran
    sed -i -e 's/Z'\''123'\'', Z'\''234'\'', Z'\''345'\'', Z'\''456'\''/int\(Z'\''123'\''\), int\(Z'\''234'\''\), int\(Z'\''345'\''\), int\(Z'\''456'\''\)/g' check_stream.F90
    sed -i -e 's/avec = Z'\''9908b0df'\''/avec = int\(Z'\''9908b0df'\'', kind=INT32\)/g' f_jump_ahead_coeff/f_jump_coeff.F90
    sed -i -e 's/ZE = Z'\''eeeeeeee'\''/ZE = int\(Z'\''eeeeeeee'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
    sed -i -e 's/ZC = Z'\''cccccccc'\''/ZC = int\(Z'\''cccccccc'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
    sed -i -e 's/Z8 = Z'\''88888888'\''/Z8 = int\(Z'\''88888888'\'', kind=INT32\)/g' f_jump_ahead_coeff/gf2xe.F90
    sed -i -e 's/dc = Z'\''0'\''/dc = int\(Z'\''0'\'', kind=INT64\)/g' f_jump_ahead_coeff/gf2xe.F90
    make
    cd $build_dir
fi
