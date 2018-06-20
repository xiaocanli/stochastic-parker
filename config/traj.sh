#!/bin/bash

conf=conf.dat
mpi_size=16
nptl=20000
ts=50
te=236
dist_flag=1      # 0 for Maxwellian. 1 for delta function.
split_flag=0     # 0 for without particle split, 1 for with split
whole_mhd_data=1 # whether to read the whole MHD data

track_particle_flag=1 # whether to track particles
nptl_selected=100     # number of selected particles to track
nsteps_interval=100    # steps interval to track particles

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
    srun -n $mpi_size --ntasks-per-node=8 ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts \
         -te $te -df $dist_flag -wm $whole_mhd_data -sf $split_flag \
         -tf $track_particle_flag -ns $nptl_selected -ni $nsteps_interval \
         -dd $diagnostics_directory -cf $conf
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
    change_variable kpara0 0.01
    change_variable kret 0.03
}

stochastic () {
    # diagnostics_directory=data/traj_data/$2/p000_b000_00001_100/
    # run_stochastic 0 1.3333333 0 $1 0.0001 1.0 $diagnostics_directory
    ## diagnostics_directory=data/traj_data/$2/p133_b000_00001_100/
    ## run_stochastic 1 1.3333333 0 $1 0.0001 1.0 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p133_b000_00001_005/
    # run_stochastic 1 1.3333333 0 $1 0.0001 0.05 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p133_b000_00001_001/
    # run_stochastic 1 1.3333333 0 $1 0.0001 0.01 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p000_b000_0001_100/
    # run_stochastic 0 1.3333333 0 $1 0.001 1.0 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p133_b000_0001_100/
    # run_stochastic 1 1.3333333 0 $1 0.001 1.0 $diagnostics_directory
    ## diagnostics_directory=data/traj_data/$2/p133_b000_0001_005/
    ## run_stochastic 1 1.3333333 0 $1 0.001 0.05 $diagnostics_directory
    ## diagnostics_directory=data/traj_data/$2/p133_b000_0001_001/
    ## run_stochastic 1 1.3333333 0 $1 0.001 0.01 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p000_b000_001_100/
    # run_stochastic 0 1.3333333 0 $1 0.01 1.0 $diagnostics_directory
    # diagnostics_directory=data/traj_data/$2/p133_b000_001_100/
    # run_stochastic 1 1.3333333 0 $1 0.01 1.0 $diagnostics_directory
    ## diagnostics_directory=data/traj_data/$2/p133_b000_001_005/
    ## run_stochastic 1 1.3333333 0 $1 0.01 0.05 $diagnostics_directory
    ## diagnostics_directory=data/traj_data/$2/p133_b000_001_001/
    ## run_stochastic 1 1.3333333 0 $1 0.01 0.01 $diagnostics_directory
    diagnostics_directory=data/traj_data/$2/p000_b000_0003_100/
    run_stochastic 0 1.3333333 0 $1 0.003 1.0 $diagnostics_directory
}

mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/S1E5_beta01_bg00/bin_data/
run_name=S1E5_beta01_bg00
stochastic $mhd_run_dir $run_name
