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

# Load parameters from a config file
inputdir = config["input"]
inputdirall = config["input_all"]
outputdir = config["output"]
genomedir = os.path.dirname(config["genome"])+"/"
genomefile =config["genome"]
chromsfile = config["chromsfile"]
chromcanonicalfile = config["chromcanonicalfile"]
tools = config['tools']
ARGS = config['args']

# Get sample names and which lane corresponds to which sample
# /path/samplename_A
# /path/samplename_A/A_lane1.fq.gz
# /path/samplename_A/A_lane2.fq.gz
import os
dirs_in_input_path = [x for x in os.listdir(inputdir) if os.path.isdir(inputdir+x)]
SAMPLES_LANES = {}
for adir in dirs_in_input_path:
    files_1 =  [x for x in glob.glob(inputdir + adir + "/*_1.fq.gz")]
    files_basenames = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in files_1]
    SAMPLES_LANES[adir] = files_basenames
SAMPLES = sum(SAMPLES_LANES.values(), [])

# Get info about chromosomes
CHROMS = [line.rstrip('\n').split("\t")[0] for line in open(chromsfile)] 
CHROMS_LENGTH = [line.rstrip('\n').split("\t")[1] for line in open(chromsfile)]
CHROMS_CANON = [line.rstrip('\n').split("\t")[0] for line in open(chromcanonicalfile)]

ASSEMBLY=config["assembly"]
MINCOV=ARGS["MINCOV"]
MINQUAL=ARGS["MINQUAL"]
MINCOV_FILTER=ARGS["MINCOV_FILTER"]

SUBSET_READS = ARGS['subset_reads']==True
NOTRIMMING = ARGS['notrimming']==True

# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = '/fast/users/kwreczy_m/work/projects/makeWGBSnake/Scripts/'
DIR_seg = outputdir+'08_segmentation/'
DIR_bigwig      = outputdir+'09_bigwig_files/'

if ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  DIR_methcall    = outputdir+'06_methyl_calls_bwameth/'
  DIR_deduped_picard     = outputdir+'05_deduplication_bwameth/'
  DIR_mapped      = outputdir+'04_mapping_bwameth/' 
  DIR_mapped_sample      = outputdir+'04_mapping_bwameth_persample/'
  DIR_mapped_persample_filtered      = outputdir+'04_mapping_bwameth_persample_filtered/'
  DIR_bigwig = outputdir+'07_bigwig_files_bwameth/'
  
if ARGS["run_bismark"] and not ARGS["run_bwameth"]:
  DIR_methcall    = outputdir+'06_methyl_calls_bismark/'
  DIR_deduped_picard     = outputdir+'05_deduplication_bismark/' 
  DIR_mapped      = outputdir+'04_mapping_bismark_nohardtrimming/'    
  DIR_mapped_sample      = outputdir+'04_mapping_bismark_sample/'

if ARGS["run_bsmap"] and not ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  DIR_methcall    = outputdir+'06_methyl_calls_bsmap/'
  DIR_deduped_picard     = outputdir+'05_deduplication_bsmap/' 
  DIR_mapped      = outputdir+'04_mapping_bsmap/'    
  DIR_mapped_sample      = outputdir+'04_mapping_bsmap_sample/'
    
DIR_posttrim_QC = outputdir+'03_posttrimming_QC/'
DIR_trimmed     = outputdir+'02_trimming/' 
DIR_rawqc       = outputdir+'01_raw_QC/'


if SUBSET_READS:
  DIR_trimmed_subset=outputdir+'subset_reads/'



# ==========================================================================================
# Output files
# ==========================================================================================

# Construct all the files we're eventually expecting to have.
FINAL_FILES = []


# ==========================================================================================
# ==========================================================================================
# Output: Pre-alignment


# # # Subset reads
# FINAL_FILES.extend(
#    expand(DIR_trimmed_subset+"{sample}/{sample}_{ext}_val_{ext}.fq.gz", sample=SAMPLES, ext=["1", "2"])
# )
# DIR_input_subset = DIR_trimmed_subset
# FINAL_FILES.extend(
#    expand(DIR_input_subset+"{sample}_1.fq.gz", sample=SAMPLES, ext=["1", "2"])
# )

# FINAL_FILES.extend(
#   expand(DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )

# # trim
# FINAL_FILES.extend(
#     expand(DIR_trimmed+"{sample}/{sample}_{ext}_val_{ext}.fq.gz", sample = SAMPLES, ext=["1", "2"])
# )#

# # Fastqc afater trimming
# FINAL_FILES.extend(
#    expand(DIR_posttrim_QC+"{sample}/{sample}_{ext}_val_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )#


# ==========================================================================================
# ==========================================================================================
# Output: Alignment


if ARGS["run_bismark"] and not ARGS["run_bwameth"]:

  # # # # Create genome bisulfite index
  # FINAL_FILES.extend(
  #     expand(genomedir+"Bisulfite_Genome/{din}_conversion/genome_mfa.{din}_conversion.fa",din=["CT","GA"])
  # )
  # 
  #
  # Alignment
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bam",sample=SAMPLES)
  )
  # Sorting
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
  )
  # # Align unmapped reads as sinle-end
  # FINAL_FILES.extend(
  #     expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}.bam",sample=SAMPLES, ext=["1", "2"])
  # )
  FINAL_FILES.extend(
     expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}_sorted.bam",sample=SAMPLES, ext=["1", "2"])
  )
  # Merge PE and SE reads
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}_sorted_merged.bam", sample=SAMPLES)
  )
  # Alignment stats
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}_merged.flagstat.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}_merged.stats.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}_merged.idxstats.txt",sample=SAMPLES)
  )
  pass

if ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  
  # # Create genome bisulfite index
  # FINAL_FILES.extend(
  #    expand(genomefile+".bwameth.{ext}",ext=["c2t.sa","c2t.amb","c2t.ann","c2t.pac","c2t.bwt","c2t"])
  # )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bwameth.bam",sample=SAMPLES)
  )
  FINAL_FILES.extend(
     expand(DIR_mapped+"{sample}/{sample}.bwameth_sorted.bam",sample=SAMPLES)
  )
  
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bwameth.flagstat.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bwameth.stats.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bwameth.idxstats.txt",sample=SAMPLES)
  )
    

if ARGS["run_bsmap"]:
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bsmap.bam",sample=SAMPLES)
  )   
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bsmap.flagstat.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bsmap.stats.txt",sample=SAMPLES)
  )
  FINAL_FILES.extend(
      expand(DIR_mapped+"{sample}/{sample}.bsmap.idxstats.txt",sample=SAMPLES)
  )



def read_treatment(file_path):
  """Read a two columns file. The first column indicates fastq_file names,
  the second sample names used to merge fastq files from the same sample but 
  different lanes into one file and the third refers to treatments for each sample.
  For example:
  WYKWK3-RUNID-FLOWCELL-LANE-8_2.fq.gz    WYKWK3 1
  WYKWK3-RUNID-FLOWCELL-LANE-8_1.fq.gz    WYKWK3 1
  WYKWK3-RUNID-FLOWCELL-LANE-7_2.fq.gz    WYKWK3 1
  WYKWK3-RUNID-FLOWCELL-LANE-7_1.fq.gz    WYKWK3 1
  WYKWK3-RUNID-FLOWCELL-LANE-8_1.fq.gz    WYKWK3 1
  WYKWK3-RUNID-FLOWCELL-LANE-8_2.fq.gz    WYKWK3 1
  It return a dict which keys indicate sample names and values treatment. For example:
  {'WYKWK3': 1]}
  """
  treatment_dict = dict()
  lanes_file=open(file_path, "r").readlines()
  for x in lanes_file:
    row=x.rstrip().split("\t")
    treatment_dict [ row[1] ] = row[2] # treatment value is a character 
  return treatment_dict
  from itertools import chain
SAMPLES_LANES_copy = SAMPLES_LANES
SAMPLES_LANES = {}
for key, value in SAMPLES_LANES_copy.items():
    if value[0] in SAMPLES:
        SAMPLES_LANES[key] = value

SAMPLES_TREATMENT = read_treatment(config['lanes_file'])



# ==========================================================================================
# ==========================================================================================
# Output: Post-alignment


if ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  # Deduplicate
  # FINAL_FILES.extend(
  #   expand(DIR_deduped_picard+"{sample}/{sample}.bwameth.dedup.sorted.bam",sample=SAMPLES_LANES.keys())
  # )
  #FINAL_FILES.extend(
  #   expand(DIR_methcall+"{sample}/tabix_{context}/{sample}.txt.bgz",
  #          sample=SAMPLES_LANES.keys(),
  #          context=[#'CpG', 
  #                   'CHG',
  #                   'CHH'
  #                   ])
  #)
  # FINAL_FILES.extend(
  #    expand(DIR_methcall+"{sample}/tabix_{context}/{sample}_{context}_filtered.txt.bgz",
  #           sample=SAMPLES_LANES.keys(),
  #           context=['CpG'#, 
  #                    #'CHG',
  #                    #'CHH'
  #                    ])
  # )
  # )

  # FINAL_FILES.extend(
  #    expand(DIR_methcall+"{sample}/tabix_CpG/{sample}_methyldacker_cpg_filtered.txt.bgz",sample=SAMPLES_LANES.keys())
  # )
  # # Mbias plot
  # FINAL_FILES.extend(
  #    expand(DIR_methcall+"{sample}/{sample}_mbias_methyldackel.txt",sample=SAMPLES_LANES.keys())
  # )
  pass



# # unite_meth_calls
# FINAL_FILES.extend(
#     [DIR_methcall+"methylBase/methylBase_CpG_dT.RDS"]
# )

# # BigWig files
# FINAL_FILES.extend(
#     expand(DIR_bigwig+"{sample}/{sample}.bw", sample=SAMPLES_LANES.keys())
# )




# This part won't work most likely:
# 
# # # Segmentation
# FINAL_FILES.extend(
#  expand(DIR_seg+"{sample}/{sample}.segments.bed", sample=SAMPLES_LANES.keys())
# )

# 
# # Differential methylation between subgroups, pairwise
# FINAL_FILES.extend(
#     expand(DIR_diffmeth+'diffmeth_{treat}.RDS', treat=TREATMENT_UNIQUE)
# )



rule target:
  input: FINAL_FILES


# ==========================================================================================
# Snakemake rules
# ==========================================================================================


  
# # ==========================================================================================
# # Segmentation:
# #
# # rules for PE and SE+PE
include: "./Rules/Segmentation_rules.py"
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
# rules for PE and SE+PE
include: "./Rules/Export_BW.py"
# 
#        
# # ==========================================================================================
# # Methylation calling:
# # 
# # rules for PE and SE+PE
if ARGS["run_bwameth"]:
  include: "./Rules/Meth_preprocessing_methyldackel_rules.py"
if ARGS["run_bismark"]:
  include: "./Rules/Meth_preprocessing_methylKit_rules.py"

# ==========================================================================================
# Deduplication :  
include: "./Rules/Deduplication_picard_rules.py"

# ==========================================================================================
# Mapping :  

# # treat unaligned reads as single-end, map the into a genome and merge them to a bam
# # file with aligned reads
# Process unaligned reads + deduplication
if ARGS["run_bismark"]:
  include: "./Rules/Unaligned_rules.py"
  include: "./Rules/Align_bismark_rules.py"
if ARGS["run_bwameth"]:
  include: "./Rules/Align_bwameth_rules.py"
if ARGS["run_bsmap"]:
  include: "./Rules/Align_bsmap_rules.py"

           
# ==========================================================================================
# Subset reads:   

# # reformat.sh is a part bbmap
# if SUBSET_READS:
# 
#   rule subset_reads_trimmed_pe:
#       input:
#         i1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
#         i2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
#       output:
#         o1 = DIR_trimmed_subset+"{sample}_1_val_1.fq.gz", #############
#         o2 = DIR_trimmed_subset+"{sample}_2_val_2.fq.gz", #############
#       params:
#           reads = " reads=100000 ",
#           sampleseed=" sampleseed=1564 "
#       log:
#           DIR_trimmed_subset+'{sample}/reformat_{sample}.log'
#       message: "Subset reads with reformat.sh"
#       shell:
#           "reformat.sh {params} in={input.i1} in2={input.i2} out={output.o1} out2={output.o2} overwrite=true 2> {log}.err"

  # rule subset_reads_pe:
  #     input:
  #       i1 = inputdir+"{sample}_1.fq.gz",
  #       i2 = inputdir+"{sample}_2.fq.gz",
  #     output:
  #       o1 = DIR_input_subset+"{sample}_1.fq.gz",
  #       o2 = DIR_input_subset+"{sample}_2.fq.gz",
  #     params:
  #         reads = " reads=100000 ",
  #         sampleseed=" sampleseed=1564 "
  #     log:
  #         DIR_trimmed_subset+'{sample}/reformat_{sample}.log'
  #     message: "Subset reads with reformat.sh"
  #     shell:
  #         "reformat.sh {params} in={input.i1} in2={input.i2} out={output.o1} out2={output.o2} overwrite=true 2> {log}.err"  


# ==========================================================================================
# Pre-mapping rules:

# Paired-end
include: "./Rules/Prealign_rules.py"

# # Single-end
# include: "./Rules/Prealign_rules_SE.py"



