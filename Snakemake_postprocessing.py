#!/usr/bin/env python

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

import glob, os, re
from snakemake.utils import R

# ./per_chrom
# ./22X3H1_unmapped_*.sam
# ./*.fq.gz ###########
# 22X3H1.bam
# 22X3H1_1_val_1_bismark_bt2_pe.bam
# 22X3H1_chr*_merged.bam
# 22X3H1_merged.bam 

# Chce zeby zostaly 
# 22X3H1_sorted.bam
# 22X3H1_sorted_merged.bam
# 22X3H1_unmapped_2_bismark_bt2_SE_report.txt
# 22X3H1_unmapped_1_bismark_bt2_SE_report.txt
# 22X3H1_unmapped_1_bismark_bt2.bam
# 22X3H1_unmapped_2_bismark_bt2.bam

# for x in $(ls ./); do echo rm -r ./$x/$x_unmapped_*.sam ./$x/$x.bam ./$x/$x_merged.bam ./$x/$x_chr*_merged.bam ./$x/$x_1_val_1_bismark_bt2_pe.bam; done
# for x in $(ls ./); do echo rm -r ./$x/$x_chr*_merged.bam; done

#a='ZTFVN7'
#cd $a
#rm $a'_1_val_1_bismark_bt2_pe.bam'* $a'_merged.bam'* $a'_chr*_merged.bam' $a'_chr*_merged_sorted.bam'*
#rm -r ./*log ./5QCV6R_unmapped_*.sam 5QCV6R_1_val_1_bismark_bt2_pe.bam* 5QCV6R_merged.bam* 5QCV6R_chr*_merged.bam 5QCV6R_chr*_merged_sorted.bam* ./per_chrom/

#rm WYKWK3_merged.bam* WYKWK3_unmapped_1_bismark_bt2.bam WYKWK3_unmapped_2.bam WYKWK3_chr*_merged.bam *temp*  WYKWK3_chr*_merged_sorted.bam* *.log.err *.log

# qsub /data/akalin/Projects/BIH_Neuroblastoma/Project/Scripts/merge_3_lanes.sh
# /data/akalin/Projects/BIH_Neuroblastoma/Project/Data/Raw/AS-195434-LR-30888_R1.fastq.gz 
# /data/akalin/Projects/BIH_Neuroblastoma/Project/Data/Raw/AS-195434-LR-30889_R1.fastq.gz
# /data/akalin/Projects/BIH_Neuroblastoma/Project/Data/Raw/AS-195434-LR-30890_R1.fastq.gz 
# /scratch/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Data/Raw_merged_lanes/22X3H1_1.fastq.gz

#!/bin/bash
#$ -pe smp 5 
#$ -l h_vmem=5G

# dir="/fast/users/kwreczy_m/scratch/Altuna_secret_project/raw_data/H26/"
# cat $dir/H26_RDML00375_HMYFJCCXY_L8_1.fq.gz  $dir/H26_RDML00375_HNCVKCCXY_L8_1.fq.gz  $dir/H26_RDML00375_HT2FHCCXY_L8_1.fq.gz H26_RDML00375_HNCVKCCXY_L7_1.fq.gz  $dir/H26_RDML00375_HT2FHCCXY_L7_1.fq.gz > $dir/H26_1.fq.gz
# cat $dir/H26_RDML00375_HMYFJCCXY_L8_2.fq.gz  $dir/H26_RDML00375_HNCVKCCXY_L8_2.fq.gz  $dir/H26_RDML00375_HT2FHCCXY_L8_2.fq.gz H26_RDML00375_HNCVKCCXY_L7_2.fq.gz  $dir/H26_RDML00375_HT2FHCCXY_L7_2.fq.gz > $dir/H26_2.fq.gz
# 
# dir="/fast/users/kwreczy_m/scratch/Altuna_secret_project/raw_data/H28/"
# cat $dir/H28_RDML00377_HMYFJCCXY_L8_1.fq.gz  $dir/H28_RDML00377_HNCLHCCXY_L3_1.fq.gz  $dir/H28_RDML00377_HT2FHCCXY_L7_1.fq.gz $dir/H28_RDML00377_HNCLHCCXY_L2_1.fq.gz  $dir/H28_RDML00377_HNCLHCCXY_L4_1.fq.gz > $dir/H28_1.fq.gz
# cat $dir/H28_RDML00377_HMYFJCCXY_L8_2.fq.gz  $dir/H28_RDML00377_HNCLHCCXY_L3_2.fq.gz  $dir/H28_RDML00377_HT2FHCCXY_L7_2.fq.gz $dir/H28_RDML00377_HNCLHCCXY_L2_2.fq.gz  $dir/H28_RDML00377_HNCLHCCXY_L4_2.fq.gz > $dir/H28_2.fq.gz


#  snakemake -s ~/projects/makeWGBSnake/Snakemake_postprocessing.py --keep-going -j 24  --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/cluster_wgbs_subset.json --printshellcmds    --cluster "qsub -V -cwd -b y -P medium -l h_vmem=15g -l h_rt=168:00:00 -pe smp 12 -N '{rule}_{wildcards.sample}'" 
##########


# deduplication
# a='ZTFVN7'
# cd $a
# rm ./per_chrom/*.log* ./per_chrom/*.dedup.bam
# rm $a'_dedup.bam' $a'_dedup.err' $a'_deduplication.log' $a'_deduplication.log.err' $a'_dedup.log'
# cd ..

#id="LJ8KVQ"
#cd $id
#rm $id'_'chr*_merged.bam $id'_1_val_1_bismark_bt2_pe.bam' $id'.bam' $id'_merged.bam' $id'_unmapped_2.sam' $id'_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam' $id'_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam' $id'_unmapped_2.bam' $id'_unmapped_1.bam' $id'_merged.bam' $id'_merged.bam.bai' $id'_merged_sort.log' $id'_merged_sort.log.err'
#rm LJ8KVQ_chr*_merged.bam LJ8KVQ_chr*_merged_sorted.bam*
#rm -r ./sambamba* ./per_chrom/
#cd ..
# ==========================================================================================
# Input/output parameters
# ==========================================================================================

# Load params from a config file
inputdir = config["input"]
outputdir = config["output"]
genomedir = config["genome"]
chromsfile = config["chromsfile"]
chromcanonicalfile = config["chromcanonicalfile"]
envs = config["env"]
tools = config['tools']
ARGS = config['args']

try:
    SAMPLES = config["samples"]
except KeyError:
    SAMPLES = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]

try:
  TREATMENT = config['treatment']
  TREATMENT_UNIQUE = set(TREATMENT)
  SAMPLE_TREAT_DICT = dict(zip(SAMPLES, TREATMENT))
except KeyError:
  print("No treatment supplied")
  TREATMENT, TREATMENT_UNIQUE,AMPLE_TREAT_DICT = [],[],[]
  
CHROMS = [line.rstrip('\n').split("\t")[0] for line in open(chromsfile)] 
CHROMS_LENGTH = [line.rstrip('\n').split("\t")[1] for line in open(chromsfile)]
CHROMS_CANON = [line.rstrip('\n').split("\t")[0] for line in open(chromcanonicalfile)]

ASSEMBLY=config["assembly"]
MINCOV=ARGS["MINCOV"]
MINQUAL=ARGS["MINQUAL"]


# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = '/fast/users/kwreczy_m/projects/makeWGBSnake/Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'07_bigwig_files/'
DIR_methcall    = outputdir+'06_methyl_calls/'
DIR_methcall_tabix    = outputdir+'06_methyl_calls/Tabix/'
DIR_deduped_picard     = outputdir+'05_deduplication/'
DIR_mapped      = outputdir+'04_mapping/'
DIR_posttrim_QC = outputdir+'03_posttrimming_QC/'
DIR_trimmed     = outputdir+'02_trimming/'
DIR_rawqc       = outputdir+'01_raw_QC/'
DIR_bam_per_chrom = DIR_mapped+'bam_per_chr/' 
DIR_seg = outputdir+'08_segmentation/'
DIR_diffmeth    = outputdir+'differential_methylation/'
DIR_ucsc_hub = outputdir+"09_ucsc_hub/"
DIR_multiqc = outputdir+"multiqc/"
DIR_ucschub = outputdir+"ucsc_hub/"


# ==========================================================================================
# Output files

# No rule to produce sort_index_bam_mapped (if you use input functions make sure that they don't raise unexpected exceptions).

# Construct all the files we're eventually expecting to have.
FINAL_FILES = []


# # FASTQC
# FINAL_FILES.extend(
#   expand(DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )
# 
# # Fastqc afater trimming
# FINAL_FILES.extend(
#    expand(DIR_posttrim_QC+"{sample}/{sample}_{ext}_val_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )
# 
# # Alignment
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )

# # Split files per chromosome
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam", sample=SAMPLES, chrom=CHROMS_CANON)
# )
# 
# # Sorting
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )


# #Align unmapped reads as sinle-end
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}.bam",sample=SAMPLES, ext=["1", "2"])
# )
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}_sorted.bam",sample=SAMPLES, ext=["1", "2"])
# )
# ## Single-end
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_sorted.bam",sample=SAMPLES)
# )
#
# # Merge PE and SE reads
# FINAL_FILES.extend(
#   expand(DIR_mapped+"{sample}/{sample}_sorted_merged.bam", sample=SAMPLES)
# )
#
# # Sort merged PE and SE reads
# FINAL_FILES.extend(
#   expand(DIR_mapped+"{sample}/{sample}_{chrom}_merged_sorted.bam", sample=SAMPLES, chrom=CHROMS_CANON)
# )
#
# # Deduplicate
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.dedup.sorted.bam",sample=SAMPLES, chrom=CHROMS_CANON)
# )
# Deduplicate PE+SE
FINAL_FILES.extend(
  expand(DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}_merged.dedup.sorted.bam",sample=SAMPLES, chrom=CHROMS_CANON)
)

# # Merge deduplicated reads
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam",sample=SAMPLES)
# )
# # Methylation call. files
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_cpg.txt.bgz",sample=SAMPLES, chrom=CHROMS_CANON)
# )
# # Methylation call. files PE+SE
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_merged_cpg.txt.bgz",sample=SAMPLES, chrom=CHROMS_CANON)
# )
# 
# # Merge methylation calling
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz",sample=SAMPLES)
# )
# # Merge methylation calling PE+SE
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/{sample}_merged_cpg_filtered.txt.bgz",sample=SAMPLES)
# )
# 
# # Unite methylation calling PE+SE
# FINAL_FILES.extend(
#    [DIR_methcall+"methylBase_per_chrom/methylBase_merged_cpg_dF.txt.bgz"]
# )

# 
# # Segmentation
# ##FINAL_FILES.extend(
# ##  expand(DIR_seg+"{sample}/{sample}.deduped_meth_segments.bed", sample=SAMPLES)
# ##)
# FINAL_FILES.extend(
#   expand(DIR_seg+"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments.bed", sample=SAMPLES, chrom=CHROMS_CANON)
# )
# FINAL_FILES.extend(
#   expand(DIR_seg+"{sample}/per_chrom/{sample}_{chrom}_merged.deduped_meth_segments.bed", sample=SAMPLES, chrom=CHROMS_CANON)
# )
# FINAL_FILES.extend(
#   expand(os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments.bed"), sample=SAMPLES)
# )

# 
# # Differential methylation between subgroups, pairwise
# FINAL_FILES.extend(
#     expand(DIR_diffmeth+'diffmeth_{treat}.RDS', treat=TREATMENT_UNIQUE)
# )
# 
# # BigWig files
# FINAL_FILES.extend(
#     expand(DIR_bigwig+"{sample}/{sample}.bw", sample=SAMPLES)
# )
# ## # BigWig files
# ## FINAL_FILES.extend(
# ##     expand(DIR_bigwig+"{sample}/per_chrom/{sample}_{chrom}.bw", sample=SAMPLES, chrom=CHROMS_CANON)
# ## )
# # BigWig files
# FINAL_FILES.extend(
#     expand(DIR_bigwig+"{sample}/{sample}_merged.bw", sample=SAMPLES)
# )


rule target:
  input: FINAL_FILES


# ==========================================================================================
# Snakemake rules
# ==========================================================================================

# ==========================================================================================
# Create a UCSC hub:

# rule ucsc_hub:
#   input:
#     methylation
#     diffmeth
#     segmentation
#   output:
#     genomes = DIR_ucschub + 'genomes.txt',
#     hub = DIR_ucschub + 'hub.txt',
#     trackDB = DIR_ucschub + ASSEMBLY + 'trackDb.txt',
#   params:
#     genome_name = ASSEMBLY,
#     paths       = DIR_ucschub,
#   log:
#     logfile = os.path.join(PATH_LOG, 'UCSC_HUB.log')
#   message:"""
#                 Running: UCSC_HUB:
#                     output: {output.hub}
#             """
#   shell:
#     "{tools}/Rscript {DIR_scripts}/Make_UCSC_HUB.R {input} {output}"
#   
#   
  
  
# ==========================================================================================
# Segmentation:
#
# rules for PE and SE+PE
#include: "./Rules/Segmentation_rules.py"
include: "./Rules/Segmentation_merged_rules.py"

# ==========================================================================================
# Differential methylation in a pairwise manner
#
# rules for PE and SE+PE
include: "./Rules/DMC_pairwise.py"


# ==========================================================================================
# Export a bigwig file:
#
# rules for PE and SE+PE
include: "./Rules/Export_BW.py"

       
# ==========================================================================================
# Methylation calling:
# 
# rules for PE and SE+PE
#include: "./Rules/Meth_preprocessing_rules.py"
include: "./Rules/Meth_preprocessing_merged_rules.py"


# ==========================================================================================
# Split bam file per chromosome:

rule split_bam_per_chr:
  input:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    temp(DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam")
  params:
    sample = "{sample}",
    chrom = "{chrom}"
  log:
    DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.log"
  # run: # this way I get "Waiting for output files 5 secs" and I can't get rid of it
  #   import re, os
  #   m = re.search('(?<=_)(.*)(?=.bam)', os.path.basename(output[0]))
  #   chrom = m.group(0).split("_").pop()
  #   mycmd = "%s/sambamba slice --output-filename=%s %s %s > ~/log.txt" % (tools, output[0], input[0], chrom)
  #   shell(mycmd)
  shell: # it cab be propably done in parallel, but this is what works for now.
    #"for chrom in {CHROMS_CANON}; do {tools}/sambamba slice --output-filename={DIR_mapped}{params.sample}/per_chrom/{params.sample}'_'$chrom'.bam' {DIR_mapped}{params.sample}/{params.sample}_sorted.bam $chrom ; done"
    "{tools}/sambamba slice --output-filename={output} {input} {params.chrom}"

# ==========================================================================================
# Process unaligned reads + deduplication

# treat unaligned reads as single-end, map the into a genome and merge them to a bam
# file with aligned reads
include: "./Rules/Unaligned_rules.py"


# ==========================================================================================
# Mapping:

rule sort_index_bam_mapped:
  input:
    DIR_mapped+"{sample}/{sample}.bam"
  output:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_sort.log"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"

rule align_pe:
     input:
         refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
         fin1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
         #qc   = [ DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.html",
         #        DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.html"]
     output:
         bam = temp(DIR_mapped+"{sample}/{sample}.bam"),
         report = DIR_mapped+"{sample}/{sample}_report.txt",
         un1 = DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz",
         un2 = DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz",
         odir = DIR_mapped+"{sample}/"
     params:
        # Bismark parameters
         bismark_args = config['args']['bismark'],
         genomeFolder = "--genome_folder " + genomedir,
         outdir = "--output_dir  "+DIR_mapped+"{sample}/",
         #nucCov = "--nucleotide_coverage",
         pathToBowtie = "--path_to_bowtie " + config['tools'],
         useBowtie2  = "--bowtie2 ",
         samtools    = "--samtools_path "+ config['tools']+'samtools',
         tempdir     = "--temp_dir "+DIR_mapped+"/{sample}"
     log:
         DIR_mapped+"{sample}/{sample}_bismark_pe_mapping.log"
     message: "Mapping paired-end reads to genome."
     run:
         commands = [
	       '{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err',
         'mv '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
         'mv '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
         'mv '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
         'mv '+output.odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
         ]
         for c in commands:
            shell(c)
            
      
# ==========================================================================================
# Pre-mapping rules:

# Paired-end
include: "./Rules/Prealign_rules.py"

# Single-end
include: "./Rules/Prealign_rules_SE.py"






