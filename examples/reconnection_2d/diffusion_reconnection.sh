#!/bin/bash

MACHINE=Perlmutter
# MACHINE=Frontera

quota_hour=12.0
restart_flag=.false.  # whether to restart a previous simulation

conf=conf_reconnection.dat
mhd_config_filename=mhd_config.dat

focused_transport=.false.

mpi_size=128
ntasks_per_node=128
size_mpi_sub=1   # Size of a MPI sub-communicator
nptl=500
nptl_max=1000000
ts=0
te=200
tmax_mhd=200         # the maximum time frame for the MHD fields
single_time_frame=0  # whether to use only a single MHD frame
time_interp=1        # whether to interpolation in-between time frames
dist_flag=1          # 0 for Maxwellian. 1 for delta function. 2 for power-law
power_index=6.2      # power-law index for the initial distribution
split_flag=1         # 0 for without particle split, 1 for with split
split_ratio=2.0      # (> 1.0 ) momentum increase ratio for particle splitting
pmin_split=2.0       # (> 1.0) minimum momentum (in p0) to start splitting particles
local_dist=.true.    # whether diagnose local particle distribution
dump_escaped_dist=.true.  # whether to dump escaped particle distributions
dump_escaped=.false.  # whether to dump escaped particles

track_particle_flag=.false. # whether to track particles
particle_tags_file=tags_selected_01.h5 # HDF5 file containing particle tags selected
nsteps_interval=100           # steps interval to track particles

inject_new_ptl=.true. # whether to inject new particles at every step
inject_same_nptl=.true. # whether to inject the same number of particles every step
tmax_to_inject=200  # the maximum time frame to inject particles
inject_part_box=.false. # inject in part of the box
ptl_xmin=-0.25
ptl_xmax=0.25
ptl_ymin=0.00
ptl_ymax=1.00
ptl_zmin=-0.25
ptl_zmax=0.25

inject_large_jz=.false. # whether to inject particles where jz is large
jz_min=500 # The minimum jz for injection
ncells_large_jz_norm=20000 # Normalization for the number of cells with large jz

inject_large_absj=.false. # whether to inject particles where rho is large
absj_min=0.2 # The minimum absj for injection
ncells_large_absj_norm=20000 # Normalization for the number of cells with large rho

inject_large_db2=.false. # whether to inject particles where db2 is large
db2_min=0.03 # The minimum db2 for injection
ncells_large_db2_norm=2000 # Normalization for the number of cells with large db2

inject_large_divv=.false. # whether to inject particles where divv is negatively large
divv_min=10.0 # The minimum divv for injection
ncells_large_divv_norm=2000 # Normalization for the number of cells with large divv

inject_large_rho=.false. # whether to inject particles where rho is large
rho_min=2.0 # The minimum rho for injection
ncells_large_rho_norm=2000 # Normalization for the number of cells with large rho

dpp_wave=0  # momentum diffusion due to wave scattering?
dpp_shear=0 # momentum diffusion due to flow shear?
weak_scattering=1 # whether in weak-scattering regime

deltab_flag=0   # whether to have spatially dependent turbulence amplitude
correlation_flag=0 # whether to have spatially dependent turbulence correlation length
ndim_field=2 # The dimension of the field
spherical_coord=0 # whether the grid is spherical
uniform_grid=1 # whether the grid is uniform
check_drift_2d=0 # whether to check particle drift in 2D
particle_data_dump=0 # whether to dump particle data

include_3rd_dim=0 # whether to include transport along the 3rd-dim in 2D simulations

acc_by_surface=0 # whether the acceleration region is separated by a surface
surface_filename1=ts_plane_top # surface 1
surface_norm1=+y # the surface norm direction
surface2_existed=.false.
is_intersection=.false.
surface_filename2=ts_plane_top # surface 2
surface_norm2=+y

varying_dt_mhd=.false. # whether the time interval for MHD fields is varying

change_variable () {
    sed -i -e "s/\($1 = \).*/\1$2/" $conf
}

if [ "$MACHINE" = "Perlmutter" ]; then
    module load cpu cray-hdf5-parallel
    export OMP_PROC_BIND=true
    export OMP_PLACES=cores
    export OMP_NUM_THREADS=2
else
    module load phdf5
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
    if [ "$restart_flag" = ".false." ] && [ "$track_particle_flag" = ".false." ]; then
        rm $diagnostics_directory/*
    fi
    mkdir -p $diagnostics_directory/restart
    tau0_scattering=$8
    duu0=$9
    particle_v0=${10}
    drift_param1=${11}
    drift_param2=${12}
    charge=${13}
    commands="./stochastic-mhd.exec \
        -qh $quota_hour -rf $restart_flag \
        -ft $focused_transport -pv $particle_v0 \
        -sm $size_mpi_sub \
        -dm $4 -mc $mhd_config_filename -np $nptl \
        -ti $time_interp -ts $ts -te $te -tm $tmax_to_inject \
        -st $single_time_frame \
        -df $dist_flag -pi $power_index \
        -sf $split_flag -sr $split_ratio -ps $pmin_split \
        -tf $track_particle_flag -ptf $particle_tags_file \
        -ni $nsteps_interval \
        -dd $diagnostics_directory -cf $conf \
        -ld $local_dist -ded $dump_escaped_dist -de $dump_escaped \
        -nm $nptl_max -in $inject_new_ptl -ij $inject_large_jz \
        -sn $inject_same_nptl -tti $tmax_to_inject -ip $inject_part_box \
        -jz $jz_min -nn $ncells_large_jz_norm -ib $inject_large_db2 \
        -db2 $db2_min -nb $ncells_large_db2_norm \
        -iv $inject_large_divv -dv $divv_min -nv $ncells_large_divv_norm \
        -ir $inject_large_rho -rm $rho_min -nr $ncells_large_rho_norm \
        -iaj $inject_large_absj -ajm $absj_min -naj $ncells_large_absj_norm \
        -xs $ptl_xmin -xe $ptl_xmax \
        -ys $ptl_ymin -ye $ptl_ymax -zs $ptl_zmin -ze $ptl_zmax \
        -dw $dpp_wave -ds $dpp_shear -t0 $tau0_scattering \
        -ws $weak_scattering -db $deltab_flag -co $correlation_flag \
        -nd $ndim_field -dp1 $drift_param1 -dp2 $drift_param2 \
        -ch $charge -sc $spherical_coord -ug $uniform_grid \
        -cd $check_drift_2d -pd $particle_data_dump \
        -i3 $include_3rd_dim -as $acc_by_surface \
        -sf1 $surface_filename1 -sn1 $surface_norm1 \
        -s2e $surface2_existed -ii $is_intersection \
        -sf2 $surface_filename2 -sn2 $surface_norm2 \
        -vdt $varying_dt_mhd -du $duu0"
if [ "$MACHINE" = "Perlmutter" ]; then
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
    pindex=1.3333333 # 3.0-turbulence spectral slope
    momentum_dependency=1
    mag_dependency=1
    dir_mhd_data=$1

    # 1keV electron, 50Gauss, slab1, Lc333km
    tau0=7.53877e-5 # scattering time for initial particles (only used for momentum diffusion)
    kpara0=0.00743592
    duu0=5578.445 # not used for Parker transport
    kret=0.01  # kperp/kpara
    particle_v0=17.20195  # particle speed / velocity normalization (not used for Parker transport)
    drift_param1=850964.408   # parameter #1 for particle drift in 3D
    drift_param2=13575468.975 # parameter #2 for particle drift in 3D
    charge=-1 # charge in unit charge
    diagnostics_directory=data/$2/transport_test_run/
    run_stochastic $momentum_dependency $pindex $mag_dependency \
        $dir_mhd_data $kpara0 $kret $diagnostics_directory $tau0 \
        $duu0 $particle_v0 $drift_param1 $drift_param2 $charge
}

mhd_run_dir=/pscratch/sd/x/xiaocan/test/transport_test/athena_reconnection_test/bin_data/
run_name=athena_reconnection_test
stochastic $mhd_run_dir $run_name
