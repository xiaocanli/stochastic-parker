#!/bin/bash

sed -i -e 's/ifort/gfortran/g' Makefile
sed -i -e 's/icpc/g++/g' Makefile
sed -i -e 's/FFLAGS := -O3 -g -traceback/FFLAGS := -O3 -g -fno-range-check/g' Makefile
sed -i -e 's/-static//g' Makefile
