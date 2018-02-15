#!/bin/bash

run_stochastic () {
    python stochastic.py --run_name $1 --run_dir $2 --ene_bin $3 $4 $5
}

plot_edist () {
    python momentum_energy_distributions.py --mhd_run $1 --sde_run $2 $3
}

plot_sdist () {
    python spatial_distribution.py --mhd_run $1 --mhd_run_dir $2 --sde_run $3 \
                                   --ene_band $4 --multi_frames --tstart $5 \
                                   --tend $6 --multi_bands
}

plot_sdist_multi () {
    for ene_bin in 0 1 2 3 4
    do
        for kappa0 in 00001 0001 001
        do
            sde_run=p000_b000_${kappa0}_100
            plot_sdist $1 $2 $sde_run $ene_bin $3 $4
            for kret in 100 005 001
            do
                sde_run=p133_b000_${kappa0}_$kret
                plot_sdist $1 $2 $sde_run $ene_bin $3 $4
            done
        done
    done
}

# mhd_run_name=S1E5_beta01_bg05
# mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/$mhd_run_name/
# run_name=$mhd_run_name/p133_b000_0001_001
# ene_bin=4
# dist_flag=--edist
# run_stochastic $run_name $mhd_run_dir $ene_bin $dist_flag $1

# mhd_run_name=S1E5_beta01_bg10
# sde_run=p133_b000_001_001
# plot_edist $mhd_run_name $sde_run $1

mhd_run=S1E5_beta01_bg00
mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/$mhd_run/
tstart=50
tend=237
plot_sdist_multi $mhd_run $mhd_run_dir $tstart $tend

mhd_run=S1E5_beta01_bg02
mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/$mhd_run/
tstart=50
tend=245
plot_sdist_multi $mhd_run $mhd_run_dir $tstart $tend

mhd_run=S1E5_beta01_bg05
mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/$mhd_run/
tstart=50
tend=246
plot_sdist_multi $mhd_run $mhd_run_dir $tstart $tend

mhd_run=S1E5_beta01_bg10
mhd_run_dir=/net/scratch3/xiaocanli/mhd/guide_field_scaling/$mhd_run/
tstart=50
tend=265
plot_sdist_multi $mhd_run $mhd_run_dir $tstart $tend
