#!/bin/sh

mhd_analysis_dir=mhd_data

mhd_code=Athena++
output_type=reconnection.prim
run_name=athena_reconnection_test
run_dir=/pscratch/sd/x/xiaocan/test/transport_test/$run_name/
config_name=athinput.reconnection
boundary=0 # 0 for periodic, 1 for open
tframe=0
tstart=0
tend=200
rx_min=0 # starting point along x (0-1)
rx_max=1 # ending point along x (0-1)
ry_min=0
ry_max=1

cd $mhd_analysis_dir

# Check if it is python 2+ or 3+
ver=$(python -V 2>&1 | sed 's/.* \([0-9]\).\([0-9]\).*/\1\2/')
if [ "$ver" -lt "27" ]; then
    pycomm=python
else
    pycomm=python3
fi

$pycomm reorganize_fields.py \
    --run_name $run_name \
    --run_dir $run_dir \
    --mhd_code $mhd_code \
    --output_type $output_type \
    --with_z_component \
    --config_name $config_name \
    --boundary $boundary \
    --tframe $tframe \
    --tstart $tstart \
    --tend $tend \
    --multi_frames \
    --rx_min $rx_min \
    --rx_max $rx_max \
    --ry_min $ry_min \
    --ry_max $ry_max \
    --organize_mhd_field
    # --xmirror \
    # --gaussian_filter \

cd -
