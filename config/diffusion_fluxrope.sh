#!/bin/bash

conf=conf_fluxrope.dat
mpi_size=512
nptl=2000000
ts=0
te=43
dist_flag=1      # 0 for Maxwellian. 1 for delta function.
split_flag=1     # 0 for without particle split, 1 for with split
whole_mhd_data=1 # whether to read the whole MHD data
local_dist=1     # whether diagnose local particle distribution

track_particle_flag=0 # whether to track particles
nptl_selected=100     # number of selected particles to track
nsteps_interval=100   # steps interval to track particles
nptl_max=10000000

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
    rm $diagnostics_directory/*
    srun -n $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts \
         -te $te -df $dist_flag -wm $whole_mhd_data -sf $split_flag \
         -tf $track_particle_flag -ns $nptl_selected -ni $nsteps_interval \
         -dd $diagnostics_directory -cf $conf -ld $local_dist \
         -nm $nptl_max
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
    change_variable kpara0 0.01
    change_variable kret 0.03
}

stochastic () {
    # diagnostics_directory=data/$2/p000_b000_00043_100/
    # run_stochastic 0 1.3333333 0 $1 0.0043 1.0 $diagnostics_directory
    # diagnostics_directory=data/$2/p100_b000_00043_100/
    # run_stochastic 1 1.3333333 0 $1 0.0043 1.0 $diagnostics_directory
    diagnostics_directory=data/$2/p000_b000_00043_002/
    run_stochastic 0 1.3333333 0 $1 0.0043 0.02 $diagnostics_directory
}

mhd_run_dir=/net/scratch3/xiaocanli/mhd/small_fluxrope_0629b/S4E4/bin_data/
run_name=fluxrope_S4E4
# stochastic $mhd_run_dir $run_name

mhd_run_dir=/net/scratch3/xiaocanli/mhd/small_fluxrope_0629b/S4E4_long/bin_data/
run_name=fluxrope_S4E4_long
stochastic $mhd_run_dir $run_name
