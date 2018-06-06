#!/bin/bash


#$ -V
#$ -j y
#$ -o logs
#$ -r yes
#$ -cwd
#$ -S /bin/bash
#$ -P control

export TMPDIR=/fast/users/${USER}/scratch/tmp
export LOGDIR=/fast/users/${USER}/scratch/tmp/logs/${JOB_ID}
mkdir -p $LOGDIR

set -x

snakemake \
    --drmaa " \
        -V \
        -cwd \
        -P medium \
        -l h_vmem=17g \
        -l h_rt=168:00:00 \ 
        -pe 7 \
        -j yes \
        -o $LOGDIR/" \
    -j 24 \
    -k \
    -p \
    -s /fast/users/kwreczy_m/projects/makeWGBSnake/Snakemake_postprocessing.py \
    --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/cluster_wgbs_hg38.json
#    --forceall
#    --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/cluster_wgbs_subset.json
#    --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/cluster_wgbs_hg38.json 
#    --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/cluster_wgbs_subset.json
#    --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/wgbs.json
