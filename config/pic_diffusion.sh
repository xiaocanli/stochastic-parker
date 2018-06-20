#!/bin/bash

conf=conf.dat
mpi_size=32
nptl=40000
ts=0
te=158
dist_flag=1      # 0 for Maxwellian. 1 for delta function.
split_flag=0     # 0 for without particle split, 1 for with split
whole_mhd_data=0 # whether to read the whole MHD data

track_particle_flag=0 # whether to track particles
nptl_selected=100     # number of selected particles to track
nsteps_interval=100   # steps interval to track particles

change_variable () {
    sed -i -e "s/\($1 = \).*/\1$2/" $conf
}

run_stochastic () {
    change_variable momentum_dependency $1
    change_variable pindex $2
    change_variable mag_dependency $3
    change_variable kpara0 $5
    change_variable kret $6
    # cat $conf
    cd ..
    diagnostics_directory=$7
    mkdir -p $diagnostics_directory
    # mpirun -np $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts -te $te -df $dist_flag
    srun -n $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts -te $te -df $dist_flag -wm $whole_mhd_data -sf $split_flag -tf $track_particle_flag -ns $nptl_selected -ni $nsteps_interval -dd $diagnostics_directory
    # mkdir -p data/$5/$6
    # mv data/*dat data/$5/$6
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
    change_variable kpara0 0.01
    change_variable kret 0.03
}

stochastic () {
    # run_stochastic 1 1.0 0 $1 $2 p100_b000
    # run_stochastic 1 1.0 1 $1 $2 p100_b100
    # diagnostics_directory=data/$2/p000_b000_0001/
    # run_stochastic 1 0.0 0 $1 0.001 0.05 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_001/
    # run_stochastic 1 1.3333333 0 $1 0.01 0.05 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b100_001/
    # run_stochastic 1 1.3333333 1 $1 0.01 0.05 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_01/
    # run_stochastic 1 1.3333333 0 $1 0.1 0.05 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b100_01/
    # run_stochastic 1 1.3333333 1 $1 0.1 0.05 $diagnostics_directory
    diagnostics_directory=data/$2/p133_b100_01_000/
    run_stochastic 1 1.3333333 1 $1 0.1 0.00 $diagnostics_directory
}

run_name=mime25_beta002_guide00_frequent_dump
mhd_run_dir=/net/scratch3/xiaocanli/reconnection/frequent_dump/$run_name/bin_data/
stochastic $mhd_run_dir $run_name
