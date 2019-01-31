#!/usr/bin/env python

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq

import glob, os, re
from snakemake.utils import R

# ==========================================================================================
# Input/output parameters
# ==========================================================================================

# Load params from a config file
inputdir = config["input"]
outputdir = config["output"]
genomedir = config["genome"]
chromsfile = config["chromsfile"]
chromcanonicalfile = config["chromcanonicalfile"]
tools = config['tools']
ARGS = config['args']

try:
    SAMPLES = config["samples"]
except KeyError:
    SAMPLES = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]

########################### TODO [START]
SAMPLES = [os.path.basename(x)[:-8] for x in glob.glob(inputdir+"*_1.fq.gz")][:10]
#SAMPLES = ["GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6"]
#SAMPLES =["YYGAV6-RUNID-0142-FLOWCELL-AHF3YTCCXY-LANE-5"]
# SAMPLES =["L1PP31-RUNID-0163-FLOWCELL-AHFMF5CCXY-LANE-7",
# "QMQHSB-RUNID-0163-FLOWCELL-AHFMF5CCXY-LANE-1",
# "KJ678J-RUNID-0160-FLOWCELL-BHF2FHCCXY-LANE-4",
# "YYGAV6-RUNID-0142-FLOWCELL-AHF3YTCCXY-LANE-5",
# "KJ678J-RUNID-0160-FLOWCELL-BHF2FHCCXY-LANE-5",
# "BZZHHV-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-2",
# "A45S2Q-RUNID-0129-FLOWCELL-BHFCTJCCXY-LANE-8",
# "GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6"]
##########################  TODO [END]


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

########################### TODO [START]
DIR_trimmed_subset=outputdir+'subset_reads/'
########################### TODO [END]

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
# )#


# # # Create genome bisulfite index
# FINAL_FILES.extend(
#    expand(genomedir+"Bisulfite_Genome/{din}_conversion/genome_mfa.{din}_conversion.fa",din=["CT","GA"])
# ) 

# # # Subset reads
# FINAL_FILES.extend(
#    expand(DIR_trimmed_subset+"{sample}/{sample}_{ext}_val_{ext}.fq.gz", sample=SAMPLES, ext=["1", "2"])
# )


# 
# Alignment
FINAL_FILES.extend(
   expand(DIR_mapped+"{sample}/{sample}.bam",sample=SAMPLES)
)
# # Sorting
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )

# # Split files per chromosome
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam", sample=SAMPLES, chrom=CHROMS_CANON)
# )


# # Align unmapped reads as sinle-end
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
# # Second mate in pbat mode
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_2_pbat_sorted.bam",sample=SAMPLES)
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
# # Deduplicate PE+SE
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}_merged.dedup.sorted.bam",sample=SAMPLES, chrom=CHROMS_CANON)
# )

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


print(FINAL_FILES)

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
  
  
# # ==========================================================================================
# # Segmentation:
# #
# # rules for PE and SE+PE
# #include: "./Rules/Segmentation_rules.py"
# include: "./Rules/Segmentation_merged_rules.py"
# 
# # ==========================================================================================
# # Differential methylation in a pairwise manner
# #
# # rules for PE and SE+PE
# include: "./Rules/DMC_pairwise.py"
# 
# 
# # ==========================================================================================
# # Export a bigwig file:
# #
# # rules for PE and SE+PE
# include: "./Rules/Export_BW.py"
# 
#        
# # ==========================================================================================
# # Methylation calling:
# # 
# # rules for PE and SE+PE
# #include: "./Rules/Meth_preprocessing_rules.py"
# include: "./Rules/Meth_preprocessing_merged_rules.py"
# 
# 
# # ==========================================================================================
# # Split bam file per chromosome:
# 
# rule split_bam_per_chr:
#   input:
#     DIR_mapped+"{sample}/{sample}_sorted.bam"
#   output:
#     temp(DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam")
#   params:
#     sample = "{sample}",
#     chrom = "{chrom}"
#   log:
#     DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.log"
#   # run: # this way I get "Waiting for output files 5 secs" and I can't get rid of it
#   #   import re, os
#   #   m = re.search('(?<=_)(.*)(?=.bam)', os.path.basename(output[0]))
#   #   chrom = m.group(0).split("_").pop()
#   #   mycmd = "%s/sambamba slice --output-filename=%s %s %s > ~/log.txt" % (tools, output[0], input[0], chrom)
#   #   shell(mycmd)
#   benchmark:
#     outputdir+"benchmarks/{sample}.split_bam_per_chr.benchmark.txt"
#   shell: # it cab be propably done in parallel, but this is what works for now.
#     #"for chrom in {CHROMS_CANON}; do {tools}/sambamba slice --output-filename={DIR_mapped}{params.sample}/per_chrom/{params.sample}'_'$chrom'.bam' {DIR_mapped}{params.sample}/{params.sample}_sorted.bam $chrom ; done"
#     "{tools}/sambamba slice --output-filename={output} {input} {params.chrom}"
# 
# # ==========================================================================================
# # Process unaligned reads + deduplication
# 
# # treat unaligned reads as single-end, map the into a genome and merge them to a bam
# # file with aligned reads
# include: "./Rules/Unaligned_rules.py"
# 
# 
# # ==========================================================================================
# # Mapping:
# 
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
  benchmark:
    outputdir+"benchmarks/{sample}.sort_index_bam_mapped.benchmark.txt"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"

rule align_pe:
     input:
         refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
         # fin1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
         # fin2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
         fin1 = DIR_trimmed_subset+"{sample}/{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed_subset+"{sample}/{sample}_2_val_2.fq.gz",
         #qc   = [ DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.html",
         #        DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.html"]
     output:
         bam = DIR_mapped+"{sample}/{sample}.bam",
         report = DIR_mapped+"{sample}/{sample}_bismark_bt2_PE_report.txt",
         un1 = DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz",
         un2 = DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
         #odir = DIR_mapped+"{sample}/"
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
     benchmark:
        outputdir+"benchmarks/{sample}.align_pe.benchmark.txt"
     run:
         sample_name = os.path.basename(input.fin1[:-14])
         output_odir = DIR_mapped+sample_name+"/"
         commands = [
	       '{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err;',
         'mv '+output_odir+sample_name+'_1_val_1_bismark_bt2_pe.bam {output.bam};',
         'mv '+output_odir+sample_name+'_1_val_1_bismark_bt2_PE_report.txt {output.report};',
         'mv '+output_odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1};',
         'mv '+output_odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2};'
         ]
         for c in commands:
           shell(c)

      
# 
# # ==========================================================================================
# # Generate methyl-converted version of the reference genome:
#            

# rule bismark_genome_preparation:
#     input:
#         ancient(genomedir)
#     output:
#         genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
#         genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
#     params:
#         pathToBowtie = "--path_to_bowtie "+ tools,
#         useBowtie2   = "--bowtie2 ",
#         verbose      = "--verbose "
#     log:
#         genomedir+'bismark_genome_preparation_'+ASSEMBLY+'.log'
#     message: "Converting {ASSEMBLY} Genome into Bisulfite analogue"
#     shell:
#         "bismark_genome_preparation {params} {input} > {log} 2> {log}.err"         
#  
 
           
# ==========================================================================================
# Subset reads:           
# indir="/fast_new/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane/02_trimming/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6/"
# f1=$indir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_1_val_1.fq.gz
# f2=$indir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_2_val_2.fq.gz
# 
# outdir="/fast_new/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6/"
# o1=$outdir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_1_val_1.fq.gz
# o2=$outdir/GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6_2_val_2.fq.gz
# 
# 
# reformat.sh reads=100000 sampleseed=1564 in=$f1 in2=$f2 out=$o1 out2=$o2 overwrite=true 2> $outdir/err.txt


rule subset_reads_pe:
    input:
      i1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
      i2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
    output:
      o1 = DIR_trimmed_subset+"{sample}/{sample}_1_val_1.fq.gz",
      o2 = DIR_trimmed_subset+"{sample}/{sample}_2_val_2.fq.gz",
    params:
        reads = " reads=100000 ",
        sampleseed=" sampleseed=1564 "
    log:
        DIR_trimmed_subset+'{sample}/reformat_{sample}.log'
    message: "Subset reads with reformat.sh"
    shell:
        "reformat.sh {params} in={input.i1} in2={input.i2} out={output.o1} out2={output.o2} overwrite=true 2> {log}.err"  



# ==========================================================================================
# Pre-mapping rules:

# Paired-end
include: "./Rules/Prealign_rules.py"

# # Single-end
# include: "./Rules/Prealign_rules_SE.py"






