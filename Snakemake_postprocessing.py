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
#SAMPLES = ["L3RXU9-RUNID-0157-FLOWCELL-AHFJKHCCXY-LANE-2"]
#SAMPLES = ["GW3LEP-RUNID-0143-FLOWCELL-BHFCTMCCXY-LANE-6"]
# WYKWK3
# SAMPLES =  ["WYKWK3-RUNID-0187-FLOWCELL-AHFCW5CCXY-LANE-8",
#            "WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-7",
#            "WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-8"]
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

SUBSET_READS = ARGS['subset_reads']==True
NOTRIMMING = ARGS['notrimming']==True

# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = '/fast/users/kwreczy_m/work/projects/makeWGBSnake/Scripts/'
#DIR_scripts   = '/home/kwreczy/projects/makeWGBSnake/Scripts/'
# DIR_plots = outputdir+'plots/'
# DIR_seg = outputdir+'08_segmentation/'
# DIR_diffmeth    = outputdir+'differential_methylation/'
# DIR_ucsc_hub = outputdir+"09_ucsc_hub/"
# DIR_multiqc = outputdir+"multiqc/"
# DIR_ucschub = outputdir+"ucsc_hub/"
# DIR_bigwig      = outputdir+'07_bigwig_files/'

if ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  DIR_methcall    = outputdir+'06_methyl_calls_bwameth/'
  DIR_deduped_picard     = outputdir+'05_deduplication_bwameth/'
  #DIR_mapped      = outputdir+'04_mapping_bwameth/'    ############
  DIR_mapped      = outputdir+'04_mapping_bwameth_nohardtrimming/'    ############
  DIR_mapped_sample      = outputdir+'04_mapping_bwameth_sample/'
  
if ARGS["run_bismark"] and not ARGS["run_bwameth"]:
  DIR_methcall    = outputdir+'06_methyl_calls_bismark/'
  DIR_deduped_picard     = outputdir+'05_deduplication_bismark/' 
  DIR_mapped      = outputdir+'04_mapping_bismark_nohardtrimming/'    
  DIR_mapped_sample      = outputdir+'04_mapping_bismark_sample/'
  

DIR_posttrim_QC = outputdir+'03_posttrimming_QC/'
#DIR_trimmed     = outputdir+'02_trimming/' 
DIR_trimmed     = outputdir+'02_trimming_nohardtrimming/' ############
DIR_rawqc       = outputdir+'01_raw_QC/'


if SUBSET_READS:
  DIR_trimmed_subset=outputdir+'subset_reads/'
#DIR_trimmed     = DIR_trimmed_subset#outputdir+'02_trimming_nohardtrimming/' ##################################




# ==========================================================================================
# Output files
# ==========================================================================================

# Construct all the files we're eventually expecting to have.
FINAL_FILES = []

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
# # FASTQC
# FINAL_FILES.extend(
#   expand(DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )

# # trim
# FINAL_FILES.extend(
#     expand(DIR_trimmed+"{sample}/{sample}_{ext}_val_{ext}.fq.gz",sample=SAMPLES, ext=["1", "2"])
# )#

# # Fastqc afater trimming
# FINAL_FILES.extend(
#    expand(DIR_posttrim_QC+"{sample}/{sample}_{ext}_val_{ext}_fastqc.html",sample=SAMPLES, ext=["1", "2"])
# )#

# ==========================================================================================
# Output: Alignment


if ARGS["run_bismark"] and not ARGS["run_bwameth"]:

# # # # Create genome bisulfite index
# FINAL_FILES.extend(
#     expand(genomedir+"Bisulfite_Genome/{din}_conversion/genome_mfa.{din}_conversion.fa",din=["CT","GA"])
# )
# 
# #
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


if ARGS["run_bwameth"] and not ARGS["run_bismark"]:
  
  pass

  # # Create genome bisulfite index
  # FINAL_FILES.extend(
  #    expand(genomefile+".bwameth.{ext}",ext=["c2t.sa","c2t.amb","c2t.ann","c2t.pac","c2t.bwt","c2t"])
  # )
  # FINAL_FILES.extend(
  #     expand(DIR_mapped+"{sample}/{sample}.bwameth.bam",sample=SAMPLES)
  # )
  # FINAL_FILES.extend(
  #    expand(DIR_mapped+"{sample}/{sample}.bwameth_sorted.bam",sample=SAMPLES)
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

### Merge lanes
def read_lanes_file(lanes_file_path):
  """Read a two columns file. The first column indicates fastq_file names and
  the second sample names used to merge fastq files from the same sample but 
  different lanes into one file. For example:
  WYKWK3-RUNID-FLOWCELL-LANE-8_2.fq.gz    WYKWK3
  WYKWK3-RUNID-FLOWCELL-LANE-8_1.fq.gz    WYKWK3
  WYKWK3-RUNID-FLOWCELL-LANE-7_2.fq.gz    WYKWK3
  WYKWK3-RUNID-FLOWCELL-LANE-7_1.fq.gz    WYKWK3
  WYKWK3-RUNID-FLOWCELL-LANE-8_1.fq.gz    WYKWK3
  WYKWK3-RUNID-FLOWCELL-LANE-8_2.fq.gz    WYKWK3
  It return a dict which keys indicate sample names and values 
  fastq files basenames. For example:
  {'WYKWK3': ['WYKWK3-RUNID-FLOWCELL-LANE-8', 'WYKWK3-RUNID-FLOWCELL-LANE-7']}
  """
  lanes_samples_dict = dict()
  lanes_file=open(lanes_file_path, "r").readlines()
  for x in lanes_file:
    row=x.rstrip().split("\t")
    # remove suffixes of fastq files
    row[0] = row[0].replace('_1.fq.gz', '').replace('_2.fq.gz', '')
    # write a dict = {"sample_name": [list_of_names_of_raw_files]}, eg:
    # {'WYKWK3': ['WYKWK3-RUNID-FLOWCELL-LANE-8', 'WYKWK3-RUNID-FLOWCELL-LANE-7']}
    try:
      if row[0] not in lanes_samples_dict [ row[1] ]:
        lanes_samples_dict [ row[1] ].append( row[0] ) 
    except KeyError:
      lanes_samples_dict [ row[1] ] = [row[0] ]
  return lanes_samples_dict

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
  
SAMPLES_LANES = read_lanes_file(config['lanes_file'])

from itertools import chain
SAMPLES_LANES_copy = SAMPLES_LANES
SAMPLES_LANES = {}
for key, value in SAMPLES_LANES_copy.items():
    if value[0] in SAMPLES:
        SAMPLES_LANES[key] = value


SAMPLES_TREATMENT = read_treatment(config['lanes_file'])

# FINAL_FILES.extend(
#    expand(DIR_mapped_sample+"{sample}/{sample}_sorted.bam",sample=SAMPLES_LANES.keys())
# )
# 
FINAL_FILES.extend(
    expand(DIR_mapped_sample+"{sample}/{sample}.flagstat.txt",sample=SAMPLES_LANES.keys())
)
FINAL_FILES.extend(
    expand(DIR_mapped_sample+"{sample}/{sample}.stats.txt",sample=SAMPLES_LANES.keys())
)
FINAL_FILES.extend(
    expand(DIR_mapped_sample+"{sample}/{sample}.idxstats.txt",sample=SAMPLES_LANES.keys())
)
# # 

# ==========================================================================================
# Output: Post-alignment


# Deduplicate
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/{sample}.bwameth.dedup.sorted.bam",sample=SAMPLES_LANES.keys())
# )
# # Methylation calling
## methylDacker
FINAL_FILES.extend(
   expand(DIR_methcall+"{sample}/{sample}_methyldacker_CpG.methylKit",sample=SAMPLES_LANES.keys()) ##########
)
FINAL_FILES.extend(
   expand(DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz",sample=SAMPLES_LANES.keys())
)
FINAL_FILES.extend(
   expand(DIR_methcall+"{sample}/tabix_CHG/{sample}.txt.bgz",sample=SAMPLES_LANES.keys())
)
FINAL_FILES.extend(
   expand(DIR_methcall+"{sample}/tabix_CHH/{sample}.txt.bgz",sample=SAMPLES_LANES.keys())
)
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/tabix_CpG/{sample}_methyldacker_cpg_filtered.txt.bgz",sample=SAMPLES_LANES.keys())
# )
# # # Mbias plot
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/{sample}_mbias_methyldackel.txt",sample=SAMPLES_LANES.keys())
# )


# Deduplicate
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam",sample=SAMPLES_LANES.keys())
# )
# # Methylation calling
## methylKit
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/tabix_{context}/{sample}_{context}_filtered.txt.bgz",
#           sample=SAMPLES_LANES.keys(),
#           context=['CpG', 'CHG','CHH'])
# )
#print(FINAL_FILES)

# # unite_meth_calls
# FINAL_FILES.extend(
#     [DIR_methcall+"methylBase/methylBase_cpg_dF.RDS"]
# )


# This part won't work most likely:
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
#include: "./Rules/Prealign_rules.py"

# # Single-end
# include: "./Rules/Prealign_rules_SE.py"



