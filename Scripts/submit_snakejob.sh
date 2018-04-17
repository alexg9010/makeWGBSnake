#!/bin/bash

#$ -V
#$ -j y
#$ -o logs
#$ -r yes
#$ -cwd
#$ -S /bin/bash
#$ -P control

export TMPDIR=/fast/users/${USER}/scratch/tmp
export LOGDIR=logs/${JOB_ID}
mkdir -p $LOGDIR

set -x

snakemake \
    --drmaa " \
        -V \
        -cwd \
        -P all \
        -l h_vmem=1g \
        -l h_rt=10:00:00 \
        -pe smp 8 \
        -j yes \
        -o $LOGDIR/" \
    -j 2 \
    -k \
    -p \
    -s /fast/users/kwreczy_m/projects/WGBS/snakemake/Snakemake.py \
    --configfile ./Config_files/wgbs.json
