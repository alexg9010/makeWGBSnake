#!/usr/bin/env python

# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
#

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
envs = config["env"]
tools = config['tools']
ARGS = config['args']
try:
    SAMPLES = config["samples"]
except KeyError:
    SAMPLES = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]

TREATMENT = config['treatment']
TREATMENT_UNIQUE = set(TREATMENT)
SAMPLE_TREAT_DICT = dict(zip(SAMPLES, TREATMENT))

CHROMS = [line.rstrip('\n').split("\t")[0] for line in open(chromsfile)] 
CHROMS_LENGTH = [line.rstrip('\n').split("\t")[1] for line in open(chromsfile)]
CHROMS_CANON = [line.rstrip('\n').split("\t")[0] for line in open(chromcanonicalfile)]

ASSEMBLY=config["assembly"]
MINCOV=ARGS["MINCOV"]
MINQUAL=ARGS["MINQUAL"]


# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = './Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'07_bigwig_files/'
DIR_methcall    = outputdir+'06_methyl_calls/'
DIR_methcall_tabix    = outputdir+'06_methyl_calls/Tabix/'
DIR_deduped_picard     = outputdir+'05_deduplication/'
DIR_mapped      = outputdir+'04_mapping/'
DIR_posttrim_QC = outputdir+'02_posttrimming_QC/'
DIR_trimmed     = outputdir+'02_trimming/'
DIR_rawqc       = outputdir+'01_raw_QC/'
DIR_bam_per_chrom = DIR_mapped+'bam_per_chr/' 
DIR_seg = outputdir+'08_segmentation/'
DIR_diffmeth    = outputdir+'differential_methylation/'
DIR_ucsc_hub = outputdir+"09_ucsc_hub/"
DIR_multiqc = outputdir+"multiqc/"


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
# # Alignment
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )
# 
# # Split files per chromosome
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam", sample=SAMPLES, chrom=CHROMS_CANON)
# )
# 
# # Sorting
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )
# 
# 
# # Align unmapped reads as sinle-end
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}_sorted.bam",sample=SAMPLES, ext=["1", "2"])
# )

# Merge PE and SE reads
FINAL_FILES.extend(
  expand(DIR_mapped+"{sample}/{sample}_sorted_merged.bam", sample=SAMPLES)
)

# Sort merged PE and SE reads
FINAL_FILES.extend(
  expand(DIR_mapped+"{sample}/{sample}_{chrom}_merged_sorted.bam", sample=SAMPLES, chrom=CHROMS_CANON)
)

# # Deduplicate
# FINAL_FILES.extend(
#   expand(DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.dedup.sorted.bam",sample=SAMPLES, chrom=CHROMS_CANON)
# )
# 
# # Methylation call. files
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_cpg.txt.bgz",sample=SAMPLES, chrom=CHROMS_CANON)
# )
# 
# # Merge methylation calling
# FINAL_FILES.extend(
#    expand(DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz",sample=SAMPLES)
# )
# 
# # Segmentation
# ##FINAL_FILES.extend(
# ##  expand(DIR_seg+"{sample}/{sample}.deduped_meth_segments.bed", sample=SAMPLES)
# ##)
# FINAL_FILES.extend(
#   expand(DIR_seg+"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments.bed", sample=SAMPLES, chrom=CHROMS_CANON)
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

rule target:
  input: FINAL_FILES


# ==========================================================================================
# Snakemake rules
# ==========================================================================================

# ==========================================================================================
# Segmentation:

include: "./Rules/Segmentation_rules.py"

# ==========================================================================================
# Differential methylation in a pairwise manner

rule diffmeth_pairwise:
     input:
         destrandTfile = DIR_methcall+"methylBase_per_chrom/methylBase_cpg_dF.txt.bgz"
     output:
         outfile=DIR_diffmeth+'diffmeth_{treat}.RDS'
     params:
         treat="{treat}",
         cores = 4, # with more than 5 usually there is not enough memory
         treatments = TREATMENT,
         sampleids = SAMPLES,
         context = "CpG",
         assembly=ASSEMBLY,
         outputdir = DIR_diffmeth,
         suffix = "diffmeth_{treat}",
         save_db = True
     shell:
         """{tools}/Rscript {DIR_scripts}/MethDiff.R \
         {input} {output} \
         {params.treat} {params.cores} "{params.treatments}" "{params.sampleids}" \
         {params.context} {params.assembly} \
         {params.outputdir} {params.suffix} {params.save_db}"""
        

# ==========================================================================================
# Export a bigwig file:

rule export_bigwig:
    input:
        seqlengths = chromcanonicalfile,
        tabixfile    = DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz"
    params:
        sampleid = "{sample}"
    output:
        bw         = DIR_bigwig+"{sample}/{sample}.bw",
    message: "Exporting bigwig files."
    shell:
       """
         {tools}/Rscript {DIR_scripts}/export_bw.R \
                 {input.tabixfile} \
                 {input.seqlengths} \
                 {ASSEMBLY} \
                 {params.sampleid} \
                 {output}
       """
       
# ==========================================================================================
# Methylation calling:
# 
include: "./Rules/Meth_preprocessing_rules.py"


# ==========================================================================================
# Deduplication:

rule sort_index_dedup_perchr:
     input:
         DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.dedup.bam"
     output:
         DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.dedup.sorted.bam"
     params:
         sort_args = config['args']['sambamba_sort'],
         tmpdir=DIR_deduped_picard+"{sample}/per_chrom/"
     log:
         DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}_sort.log"
     shell:
         "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"

rule dedup_picard_perchr:
     input:
         DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam"
     output:
        outfile=DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.dedup.bam",
        metrics = DIR_deduped_picard+"{sample}/per_chrom/{sample}_{chrom}.deduplication.metrics.txt"
     params:
         picard_MarkDuplicates_args = config['args']['picard_MarkDuplicates_args']
     log:
         DIR_deduped_picard+"{sample}/per_chrom/{sample}_deduplication.{chrom}.log"
     message:
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          """{tools}/picard MarkDuplicates I={input} O={output.outfile} \
          M={output.metrics} \
          REMOVE_DUPLICATES=true AS=true {params.picard_MarkDuplicates_args} \
          > {log} \
          2> {log}.err"""

# ==========================================================================================
# Split bam file per chromosome:

rule split_bam_per_chr:
  input:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped+"{sample}/per_chrom/{sample}_{chrom}.bam"
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
# Process unaligned reads

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
         bam = DIR_mapped+"{sample}/{sample}.bam",
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
	       #'{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
         'ln -s '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
         'ln -s '+output.odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
         ]
         for c in commands:
            shell(c)
            

# rule split:
#      input:
#          DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz"
#      output:
#          
#      params:
#          lines = 1000000
#      message:
#          "Splitting file {input} into {params.parts} parts."
#      shell:
#        "zcat tmp.fq.gz | split --lines=10000000 - bigfile-split. --numeric-suffixes  --filter='gzip > $FILE.gz'"          
#        

         
# ==========================================================================================
# Pre-mapping rules:

include: "./Rules/Prealign_rules.py"


