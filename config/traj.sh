#!/bin/bash

conf=conf.dat
mpi_size=1
nptl=1000
ts=0
te=200
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
    # cat $conf
    cd ..
    # mpirun -np $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts -te $te -df $dist_flag
    srun -n $mpi_size ./stochastic-2dmhd.exec -dm $4 -np $nptl -ts $ts -te $te -df $dist_flag -wm $whole_mhd_data -sf 0 -tf $track_particle_flag -ns $nptl_selected -ni $nsteps_interval
    # mkdir -p data/$5/$6
    # mv data/*dat data/$5/$6
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
}

stochastic () {
    # run_stochastic 0 0.0 0 $1 $2 p000_b000
    # run_stochastic 1 1.0 0 $1 $2 p100_b000
    # run_stochastic 1 1.3333333 0 $1 $2 p133_b000
    # run_stochastic 1 1.0 1 $1 $2 p100_b100
    run_stochastic 1 1.3333333 1 $1 $2 p133_b100
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
mhd_run_dir=/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01_bg02/bin_data/
run_name=athena_S1E5_beta01_bg02
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01_bg10/bin_data/
# run_name=athena_S1E5_beta01_bg10
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/S9E4/beta001_bg00/bin_data/
# run_name=S9E4_beta001_bg00
stochastic $mhd_run_dir $run_name
