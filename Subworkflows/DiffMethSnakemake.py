#!/usr/bin/env python

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
#

import glob, os, re

inputdir = config["input"]
outputdir = config["output"]
genomedir = config["genome"]
envs = config["env"]
tools = config['tools']
args = config['args']


try:
    SAMPLES = config["samples"]
except KeyError:
    SAMPLES = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]
print(config)

TREATMENT = config['treatment']
TREATMENT_UNIQUE = set(TREATMENT)
SAMPLE_TREAT_DICT = dict(zip(SAMPLES, TREATMENT))

#chrr='chroms_hg38_nohaplo.fa.txt' # hg38
chrr="chroms.txt" # hg19
CHROMS = [line.rstrip('\n') for line in open(genomedir+chrr)]
ASSEMBLY="hg19"
MINCOV=10
MINQUAL=20


WORKDIR = os.getcwd() + "/"
DIR_scripts   = './Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'07_bigwig_files/'
DIR_methcall    = outputdir+'06_methyl_calls/'
DIR_deduped     = outputdir+'05_deduplication/'
DIR_mapped      = outputdir+'04_mapping/'
DIR_posttrim_QC = outputdir+'02_posttrimming_QC/'
DIR_trimmed     = outputdir+'02_trimming/'
DIR_rawqc       = outputdir+'01_raw_QC/'
DIR_bam_per_chrom = DIR_mapped+'bam_per_chr/'
DIR_seg = outputdir+'08_segmentation/'
DIR_diffmeth    = outputdir+'differential_methylation/'
DIR_diffmeth_pairwise = outputdir+'differential_methylation_pairwise/'
DIR_ucsc_hub = outputdir+"09_ucsc_hub/"
DIR_multiqc = outputdir+"multiqc/"




# Construct all the files we're eventually expecting to have.
FINAL_FILES = [] 


FINAL_FILES.extend(
expand(DIR_methcall+'{sample}/diffmeth_{sample}_{treatment}.RDS', sample=SAMPLES, treatment = TREATMENT_UNIQUE)
)



rule diffmeth_pairwise:
   input:
     expand(DIR_methcall+'{sample}/{sample}_dedup_methylRaw.RDS', sample=SAMPLES)
   output:
     expand(DIR_methcall+'{sample}/diffmeth_{sample}_{treatment}.RDS', sample=SAMPLES, treatment= TREATMENT_UNIQUE)
   run:
     print(input[0])
     print(output[0])  



















