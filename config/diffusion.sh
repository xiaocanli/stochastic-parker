#!/bin/bash

conf=conf.dat
mpi_size=16
nptl=10000
ts=0
te=108
dist_flag=1 # 0 for Maxwellian. 1 for delta function.
change_variable () {
    sed -i -e "s/\($1 = \).*/\1$2/" $conf
}

run_stochastic () {
    change_variable momentum_dependency $1
    change_variable pindex $2
    change_variable mag_dependency $3
    cat $conf
    cd ..
    # mpirun -np $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts -te $te -df $dist_flag
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
}

stochastic () {
    run_stochastic 0 0.0 0 $1
    run_stochastic 1 1.0 0 $1
    run_stochastic 1 1.3333333 0 $1
    run_stochastic 1 1.0 1 $1
    run_stochastic 1 1.3333333 1 $1
}

mhd_run_dir=/net/scratch3/xiaocanli/mhd/S1E5/beta001/
stochastic $mhd_run_dir
