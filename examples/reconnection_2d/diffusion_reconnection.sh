#!/bin/bash

MACHINE=Cori
# MACHINE=Frontera

conf=conf.dat
mhd_config_filename=mhd_config.dat
mpi_size=32
ntasks_per_node=32
size_mpi_sub=1   # Size of a MPI sub-communicator
nptl=500
ts=0
te=200
single_time_frame=0
time_interp=1    # whether to interpolation in-between time frames
dist_flag=1      # 0 for Maxwellian. 1 for delta function. 2 for power-law
power_index=6.2  # power-law index for the initial distribution
split_flag=1     # 0 for without particle split, 1 for with split
split_ratio=2.0  # momentum increase ratio for particle splitting
local_dist=1     # whether diagnose local particle distribution

track_particle_flag=0 # whether to track particles
nptl_selected=20     # number of selected particles to track
nsteps_interval=1   # steps interval to track particles
nptl_max=1000000

inject_new_ptl=1 # whether to inject new particles at every step
inject_same_nptl=1 # whether to inject the same number of particles every step
inject_part_box=0 # inject in part of the box
ptl_xmin=0.245
ptl_xmax=0.255
ptl_ymin=0.04
ptl_ymax=0.15
ptl_zmin=0.00
ptl_zmax=0.02

inject_large_jz=0 # whether to inject particles where jz is large
jz_min=200 # The minimum jz for injection
ncells_large_jz_norm=20000 # Normalization for the number of cells with large jz

inject_large_db2=0 # whether to inject particles where db2 is large
db2_min=0.03 # The minimum db2 for injection
ncells_large_db2_norm=2000 # Normalization for the number of cells with large db2

inject_large_divv=0 # whether to inject particles where divv is negatively large
divv_min=10.0 # The minimum divv for injection
ncells_large_divv_norm=2000 # Normalization for the number of cells with large divv

dpp_wave=0  # momentum diffusion due to wave scattering?
dpp_shear=0 # momentum diffusion due to flow shear?
weak_scattering=1 # whether in weak-scattering regime

deltab_flag=0   # whether to have spatially dependent turbulence amplitude
correlation_flag=0 # whether to have spatially dependent turbulence correlation length
ndim_field=2 # The dimension of the field
drift_param1=3.315064e+08 # parameter #1 for particle drift in 3D
drift_param2=9.535508e+08 # parameter #2 for particle drift in 3D
charge=-1 # charge in unit charge
spherical_coord=0 # whether the grid is spherical
uniform_grid=1 # whether the grid is uniform
check_drift_2d=0 # whether to check particle drift in 2D
particle_data_dump=0 # whether to dump particle data

include_3rd_dim=0 # whether to include transport along the 3rd-dim in 2D simulations

acc_by_surface=0 # whether the acceleration region is separated by a surface
surface_filename1=ts_plane_bot # surface 1
surface_norm1=+y # the surface norm direction
surface2_existed=.false.
is_intersection=.true.
surface_filename2=ts_plane_top # surface 2
surface_norm2=-y

change_variable () {
    sed -i -e "s/\($1 = \).*/\1$2/" $conf
}

if [ "$MACHINE" = "Cori" ]; then
    export OMP_PROC_BIND=true
    # export OMP_PROC_BIND=spread
    export OMP_PLACES=threads
    export OMP_NUM_THREADS=1
else
    export IBRUN_TASKS_PER_NODE=$ntasks_per_node
    export OMP_NUM_THREADS=1
fi

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
    tau0_scattering=$8
    commands="./stochastic-mhd.exec -sm $size_mpi_sub \
        -dm $4 -mc $mhd_config_filename -np $nptl \
        -ti $time_interp -ts $ts -te $te -st $single_time_frame \
        -df $dist_flag -pi $power_index \
        -sf $split_flag -sr $split_ratio \
        -tf $track_particle_flag -ns $nptl_selected -ni $nsteps_interval \
        -dd $diagnostics_directory -cf $conf -ld $local_dist \
        -nm $nptl_max -in $inject_new_ptl -ij $inject_large_jz \
        -sn $inject_same_nptl -ip $inject_part_box -jz $jz_min \
        -nn $ncells_large_jz_norm -ib $inject_large_db2 \
        -db2 $db2_min -nb $ncells_large_db2_norm \
        -iv $inject_large_divv -dv $divv_min \
        -nv $ncells_large_divv_norm -xs $ptl_xmin -xe $ptl_xmax \
        -ys $ptl_ymin -ye $ptl_ymax -zs $ptl_zmin -ze $ptl_zmax \
        -dw $dpp_wave -ds $dpp_shear -t0 $tau0_scattering \
        -ws $weak_scattering -db $deltab_flag -co $correlation_flag \
        -nd $ndim_field -dp1 $drift_param1 -dp2 $drift_param2 \
        -ch $charge -sc $spherical_coord -ug $uniform_grid \
        -cd $check_drift_2d -pd $particle_data_dump \
        -i3 $include_3rd_dim -as $acc_by_surface \
        -sf1 $surface_filename1 -sn1 $surface_norm1 \
        -s2e $surface2_existed -ii $is_intersection \
        -sf2 $surface_filename2 -sn2 $surface_norm2"
if [ "$MACHINE" = "Cori" ]; then
    srun -n $mpi_size --ntasks-per-node $ntasks_per_node \
        -c $OMP_NUM_THREADS --cpu-bind=cores $commands
else
    ibrun -n $mpi_size $commands
fi
    cd config
    change_variable momentum_dependency 1
    change_variable pindex 1.3333333
    change_variable mag_dependency 1
    change_variable kpara0 0.01
    change_variable kret 0.03
}

stochastic () {
    tau0=7.548180e-05 # scattering time for initial particles
    kpara0=7.419363e-03
    diagnostics_directory=data/$2/transport_test_run/
    run_stochastic 1 1.3333333 1 $1 $kpara0 0.01 $diagnostics_directory $tau0
}

mhd_run_dir=/global/cscratch1/sd/xiaocan/athena_reconnection_test/test_periodic_bc/bin_data/
run_name=test_periodic_bc
stochastic $mhd_run_dir $run_name
