#!/usr/bin/env python

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
# This pipeline is heavily based on the PiGx BSseq pipeline github.com/BIMSBbioinfo/pigx_bsseq


# TODO:
### add divide bwa-meth output into paired data and singletons and run multiqc on it
###

import glob, os, re
from snakemake.utils import R

# ==========================================================================================
# Input/output parameters
# ==========================================================================================

# Load parameters from a config file
inputdir = config["input"]
outputdir = config["output"]
genomedir = os.path.dirname(config["genome"])+"/"
genomefile =config["genome"]
chromsfile = config["chromsfile"]
chromcanonicalfile = config["chromcanonicalfile"]
tools = config['tools']
ARGS = config['args']

try:
    SAMPLES = config["samples"]
except KeyError:
    SAMPLES = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]

########################### TODO [START]
SAMPLES = [os.path.basename(x)[:-8] for x in glob.glob(inputdir+"*_1.fq.gz")]
SAMPLES = ["22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4"]
#SAMPLES = ["GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6"]
#SAMPLES =["QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7"]
##########################  TODO [END]

try:
  TREATMENT = config['treatment']
  TREATMENT_UNIQUE = set(TREATMENT)
  SAMPLE_TREAT_DICT = dict(zip(SAMPLES, TREATMENT))
except KeyError:
  print("No treatment supplied")
  TREATMENT, TREATMENT_UNIQUE,SAMPLE_TREAT_DICT = [],[],[]
  
CHROMS = [line.rstrip('\n').split("\t")[0] for line in open(chromsfile)] 
CHROMS_LENGTH = [line.rstrip('\n').split("\t")[1] for line in open(chromsfile)]
CHROMS_CANON = [line.rstrip('\n').split("\t")[0] for line in open(chromcanonicalfile)]

ASSEMBLY=config["assembly"]
MINCOV=ARGS["MINCOV"]
MINQUAL=ARGS["MINQUAL"]

SUBSET_READS = config['args']['subset_reads']==True
NOTRIMMING = config['args']['notrimming']==True

# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = '/fast/users/kwreczy_m/projects/makeWGBSnake/Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'07_bigwig_files/'
DIR_methcall    = outputdir+'06_methyl_calls/'
DIR_methcall_tabix    = outputdir+'06_methyl_calls/Tabix/'
DIR_deduped_picard     = outputdir+'05_deduplication/'
#DIR_mapped      = outputdir+'04_mapping/'     ############'04_mapping/' ##########################################
DIR_mapped      = outputdir+'04_mapping_bwameth/'     ############'04_mapping/' ##########################################
#DIR_mapped      = outputdir+'04_mapping_bwameth_notrimming/'     ############'04_mapping/' ##########################################
#DIR_mapped      = outputdir+'04_mapping_bismark/'     ############'04_mapping/' ##########################################
DIR_posttrim_QC = outputdir+'03_posttrimming_QC/'
#DIR_trimmed     = outputdir+'/02_trimming/' 
DIR_trimmed     = outputdir+'/02_trimming_nohardtrimming/' 
DIR_rawqc       = outputdir+'01_raw_QC/'
DIR_bam_per_chrom = DIR_mapped+'bam_per_chr/' 
DIR_seg = outputdir+'08_segmentation/'
DIR_diffmeth    = outputdir+'differential_methylation/'
DIR_ucsc_hub = outputdir+"09_ucsc_hub/"
DIR_multiqc = outputdir+"multiqc/"
DIR_ucschub = outputdir+"ucsc_hub/"

########################### TODO [START]
#if SUBSET_READS:
DIR_trimmed_subset=outputdir+'subset_reads/'
DIR_input_subset=inputdir+'subset_reads/'
#DIR_mapped_bwameth      = outputdir+'04_mapping_bwameth/' ##################
#DIR_mapped      = DIR_mapped_bwameth
########################### TODO [END]





# ==========================================================================================
# Output files


# Construct all the files we're eventually expecting to have.
FINAL_FILES = []

# # FASTQC
# FINAL_FILES.extend(
#   expand(DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )
# 

# # # trim
# FINAL_FILES.extend(
#    expand(DIR_trimmed+"{sample}/{sample}_{ext}_val_{ext}.fq.gz",sample=SAMPLES, ext=["1", "2"])
# )#


# # Fastqc afater trimming
# FINAL_FILES.extend(
#    expand(DIR_posttrim_QC+"{sample}/{sample}_{ext}_val_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )#


# # # Create genome bisulfite index
# FINAL_FILES.extend(
#    expand(genomedir+"Bisulfite_Genome/{din}_conversion/genome_mfa.{din}_conversion.fa",din=["CT","GA"])
# )

# # Create genome bisulfite index
# FINAL_FILES.extend(
#    expand(genomefile+".bwameth.{ext}",ext=["c2t.sa","c2t.amb","c2t.ann","c2t.pac","c2t.bwt","c2t"])
# )
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}.bwameth.bam",sample=SAMPLES)
# )
# 
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}.bwameth.flagstat.txt",sample=SAMPLES)
# )
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}.bwameth.stats.txt",sample=SAMPLES)
# )
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}.bwameth.idxstats.txt",sample=SAMPLES)
# )



# # # # Subset reads
# FINAL_FILES.extend(
#    expand(DIR_trimmed_subset+"{sample}/{sample}_{ext}_val_{ext}.fq.gz", sample=SAMPLES, ext=["1", "2"])
# )
# FINAL_FILES.extend(
#    expand(DIR_input_subset+"{sample}_1.fq.gz", sample=SAMPLES, ext=["1", "2"])
# )

# 
# #
# # Alignment
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}.bam",sample=SAMPLES)
# )

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

# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}_merged.flagstat.txt",sample=SAMPLES)
# )
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}_merged.stats.txt",sample=SAMPLES)
# )
# FINAL_FILES.extend(
#     expand(DIR_mapped+"{sample}/{sample}_merged.idxstats.txt",sample=SAMPLES)
# )


# # Deduplicate
FINAL_FILES.extend(
  expand(DIR_deduped_picard+"{sample}/{sample}.dedup.bam",sample=SAMPLES)
)

# # methylation calling
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz",sample=SAMPLES)
# )

# # unite_meth_calls 
FINAL_FILES.extend(
    [DIR_methcall+"methylBase/methylBase_cpg_dF.RDS"]
)


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
#include: "./Rules/Segmentation_merged_rules.py"
# 
# # ==========================================================================================
# # Differential methylation in a pairwise manner
# #
# # rules for PE and SE+PE
#include: "./Rules/DMC_pairwise.py"
# 
# 
# # ==========================================================================================
# # Export a bigwig file:
# #
# # rules for PE and SE+PE
#include: "./Rules/Export_BW.py"
# 
#        
# # ==========================================================================================
# # Methylation calling:
# # 
# # rules for PE and SE+PE
include: "./Rules/Meth_preprocessing_rules.py"
#include: "./Rules/Meth_preprocessing_merged_rules.py"
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

# ==========================================================================================
# Mapping + deduplication :  

# # treat unaligned reads as single-end, map the into a genome and merge them to a bam
# # file with aligned reads
# Process unaligned reads + deduplication
include: "./Rules/Unaligned_rules.py"
include: "./Rules/Align_bismark_rules.py"
include: "./Rules/Align_bwameth_rules.py"

           
# ==========================================================================================
# Subset reads:   

# reformat.sh is a part bbmap
if SUBSET_READS:

  rule subset_reads_trimmed_pe:
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

  rule subset_reads_pe:
      input:
        i1 = inputdir+"{sample}_1.fq.gz",
        i2 = inputdir+"{sample}_2.fq.gz",
      output:
        o1 = DIR_input_subset+"{sample}_1.fq.gz",
        o2 = DIR_input_subset+"{sample}_2.fq.gz",
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






