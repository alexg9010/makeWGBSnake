#!/bin/bash

tools=/fast/users/kwreczy_m/programs/bin/miniconda3/envs/mybase/bin/
input='/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/RDS/methylBaseDB.obj_filtered_hg19_destrandT.RDS'
scripts=/fast/users/kwreczy_m/projects/makeWGBSnake/Scripts/
output=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/hg19/

for t in {1..4}
do
 # your-unix-command-here
 echo $t
 qsub -V -P highmem -N diffmeth$t -l h_rt=168:00:00 -l h_vmem=60g -pe smp 8 -b y " $tools/Rscript $scripts/MethDiff.R $input $output'differential_methylation/diffmeth_'$t'.RDS' $t $output'differential_methylation/diffmeth_'$t'.log' 2> $output'differential_methylation/diffmeth_'$t'.log.err' " 
done

