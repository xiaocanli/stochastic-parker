#!/bin/bash

conf=conf21.dat
mpi_size=32
nptl=40000
ts=50
te=51
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
    srun -n $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts \
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
    # run_stochastic 1 1.0 0 $1 $2 p100_b000
    # run_stochastic 1 1.0 1 $1 $2 p100_b100
    diagnostics_directory=data/$2/p000_b000_0001/
    run_stochastic 1 0.0 0 $1 0.001 0.05 $diagnostics_directory
    # diagnostics_directory=data/$2/p000_b000_001/
    # run_stochastic 1 0.0 0 $1 0.01 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_001/
    # run_stochastic 1 1.3333333 0 $1 0.01  $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_001/
    # run_stochastic 1 1.3333333 0 $1 0.01 1.0 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b100_001/
    # run_stochastic 1 1.3333333 1 $1 0.01  $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_01/
    # run_stochastic 1 1.3333333 0 $1 0.1   $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b100_01/
    # run_stochastic 1 1.3333333 1 $1 0.1   $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b000_0001/
    # run_stochastic 1 1.3333333 0 $1 0.001 $diagnostics_directory
    # diagnostics_directory=data/$2/p133_b100_0001/
    # run_stochastic 1 1.3333333 1 $1 0.001 $diagnostics_directory
}

# mhd_run_dir=/net/scratch3/xiaocanli/mhd/S1E5/beta001/
# run_name=S16E5_beta001_bg0_new
# stochastic $mhd_run_dir $run_name
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/S1E5/beta01/
# run_name=S16E5_beta01_bg0
# stochastic $mhd_run_dir $run_name
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/S16E5/beta001_bg10/bin_data/
# run_name=S1E4_beta001_bg05
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01/bin_data/
# run_name=athena_S1E5_beta01
# stochastic $mhd_run_dir $run_name
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01_bg02/bin_data/
# run_name=athena_S1E5_beta01_bg02
# stochastic $mhd_run_dir $run_name
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01_bg10/bin_data/
# run_name=athena_S1E5_beta01_bg10
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/S9E4/beta001_bg00/bin_data/
# run_name=S9E4_beta001_bg00
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/S1E5_beta01_bg00/bin_data/
# run_name=S1E5_beta01_bg00
# stochastic $mhd_run_dir $run_name
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/S1E5_beta01_bg02/bin_data/
# run_name=S1E5_beta01_bg02
# stochastic $mhd_run_dir $run_name
mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/S1E5_beta01_bg05/bin_data/
run_name=S1E5_beta01_bg05
stochastic $mhd_run_dir $run_name
