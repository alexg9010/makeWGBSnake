

java -ea -Xmx200m -cp /fast/users/kwreczy_m/programs/bin/miniconda3/envs/mybase/opt/bbmap-37.90/current/ jgi.ReformatReads in1=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw//AS-196016-LR-31106_R1.fastq.gz in2=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw//AS-196016-LR-31106_R2.fastq.gz out1=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_random10Mlanespe//Sampled_AS-196016-LR-31106_R1.fastq.gz out2=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_random10Mlanespe//Sampled_AS-196016-LR-31106_R2.fastq.gz samplereads=1000000 sampleseed=1564 int=t
Executing jgi.ReformatReads [in1=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw//AS-196016-LR-31106_R1.fastq.gz, in2=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw//AS-196016-LR-31106_R2.fastq.gz, out1=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_random10Mlanespe//Sampled_AS-196016-LR-31106_R1.fastq.gz, out2=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_random10Mlanespe//Sampled_AS-196016-LR-31106_R2.fastq.gz, samplereads=1000000, sampleseed=1564, int=t]

in1= samplereads=1000000 sampleseed=1564 int=t
samplereads=

#GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6
#KJ678J-RUNID-0160-FLOWCELL-BHF2FHCCXY-LANE-5
  
# http://seqanswers.com/forums/showthread.php?t=58221
# Usage:  reformat.sh in=<file> in2=<file2> out=<outfile> out2=<outfile2>
  
indir="/fast_new/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane/02_trimming/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6/"
f1=$indir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_1_val_1.fq.gz
f2=$indir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_2_val_2.fq.gz

outdir="/fast_new/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6/"
o1=$outdir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_1_val_1.fq.gz
o2=$outdir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_2_val_2.fq.gz


reformat.sh reads=100000 sampleseed=1564 in=$f1 in2=$f2 out=$o1 out2=$o2 overwrite=true 2> $outdir/err.txt


#https://community.plot.ly/t/how-to-build-dendrograms-in-r/1762
# https://community.plot.ly/t/how-to-build-dendrograms-in-r/1762

  

snakemake -s ~/work/projects/makeWGBSnake/Snakemake_postprocessing.py --keep-going -j 4  --configfile ~/work/projects/makeWGBSnake/Config_files/cluster_subset_wgbs_hg38.yaml --printshellcmds --cluster "qsub -V -cwd -b y  -l h_vmem=8g -pe smp 22 -N '{rule}_{wildcards.sample}'"


