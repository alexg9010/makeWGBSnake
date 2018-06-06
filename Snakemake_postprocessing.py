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

ASSEMBLY="hg19"
MINCOV=10 
MINQUAL=20  
NPARTS=5


# ==========================================================================================
# Output directories

WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = './Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'07_bigwig_files/'
DIR_methcall    = outputdir+'06_methyl_calls/'
DIR_methcall_tabix    = outputdir+'06_methyl_calls/Tabix/'
DIR_deduped     = outputdir+'05_deduplication/'
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

# Alignment
FINAL_FILES.extend(
   expand(DIR_mapped+"{sample}/{sample}.bam",sample=SAMPLES)
)
# 
# # # Sorting
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )


#Align unmapped reads as sinle-end
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}_sorted.bam",sample=SAMPLES, ext=["1", "2"])
# )
# 
# #Sorting
# FINAL_FILES.extend(
#    expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=SAMPLES)
# )
# 
# #Merge PE and SE reads
# FINAL_FILES.extend(
#   expand(DIR_mapped+"{sample}/{sample}_sorted_merged.bam", sample=SAMPLES)
# )


# Deduplicate
FINAL_FILES.extend(
  expand(DIR_deduped+"{sample}/{sample}_{chrom}.dedup.sorted.bam",sample=SAMPLES, chrom=CHROMS_CANON)
)

# Methylation calling
FINAL_FILES.extend(
   #expand(DIR_methcall+'{sample}/{sample}_cpg.txt.bgz', sample=SAMPLES)
   expand(DIR_methcall+"{sample}/{sample}_{chrom}_cpg.txt.bgz",sample=SAMPLES, chrom=CHROMS_CANON)
)

# # Filter meth.
FINAL_FILES.extend(
   expand(DIR_methcall+'{sample}/{sample}_{chrom}_cpg_filtered.txt.bgz', sample=SAMPLES, chrom=CHROMS_CANON)
)
 
# # Unite meth calls
FINAL_FILES.extend(
   expand(DIR_methcall+'methylBaseDB.obj_filtered_destrandF.{chrom}.RDS', chrom=CHROMS_CANON)
)
# 
# # # Segmentation
# FINAL_FILES.extend(
#   expand(DIR_seg+"{sample}/{sample}.deduped_meth_segments.bed", sample=SAMPLES)
# )
# 
# Differential methylation between subgroups, pairwise
# FINAL_FILES.extend(
#     expand(DIR_diffmeth+'diffmeth_{tret}.RDS', sample=SAMPLES, tret=TREATMENT_UNIQUE)
# )

# BigWig files
FINAL_FILES.extend(
    expand(DIR_bigwig+"{sample}/{sample}_{chrom}.bw", sample=SAMPLES, chrom=CHROMS_CANON)
)


#snakemake -s /fast/users/kwreczy_m/projects/makeWGBSnake/Snakemake_postprocessing.py -j 24  --configfile /fast/users/kwreczy_m/projects/makeWGBSnake/Config_files/test.json --printshellcmds  --forceall

print(FINAL_FILES)
rule target:
  input: FINAL_FILES


# ==========================================================================================
# Snakemake rules
# ==========================================================================================

# ==========================================================================================
# Segmentation:

# 
# rule meth_segments:
#      input:
#          inputfile     = DIR_methcall+"{sample}/{sample}_filtered.txt.bgz"
#      output:
#          grfile      = os.path.join(DIR_seg,"{sample}/{sample}.deduped_meth_segments_gr.RDS"),
#          bedfile     = os.path.join(DIR_seg,"{sample}/{sample}.deduped_meth_segments.bed")
#      params:
#          methSegPng = DIR_seg+"{sample}/{sample}.deduped_meth_segments.png",
#          assembly = ASSEMBLY,
#          sampleid = "{sample}"
#      log:
#          os.path.join(DIR_seg,"{sample}/{sample}.deduped_meth_segments.log")
#      message: "Segmenting methylation profile for {input.inputfile}."
#      shell:
#          """
#           {tools}/Rscript {DIR_scripts}/methSeg.R \
#                           {input.inputfile} \
#                           {output.grfile} \
#                           {output.bedfile} \
#                           {params.methSegPng} \
#                           {params.assembly} \
#                           {params.sampleid} \
#                           {log}
#           """


# ==========================================================================================
# Export a bigwig file:

rule export_bigwig:
    input:
        seqlengths = chromcanonicalfile,
        tabixfile    = expand(DIR_methcall+"{{sample}}/{{sample}}_{chrom}_cpg_filtered.txt.bgz", chrom=CHROMS_CANON)
    params:
        sampleid = "{sample}"
    output:
        bw         = expand(DIR_bigwig+"{{sample}}/{{sample}}_{chrom}.bw", chrom=CHROMS_CANON),
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
# Differential methylation in a pairwise manner

# rule diffmeth_pairwise:
#      input:
#          destrandTfile = DIR_methcall+"methylBaseDB.obj_filtered_destrandT.RDS"
#      output:
#          outfile=DIR_diffmeth+'diffmeth_{treat}.RDS'
#      params:
#          treatment="{treat}",
#          cors = 20
#      run:
#          R("""
#          library(methylKit)
#     
#          input = "{input.destrandTfile}"
#          output = "{output.outfile}"
#          t = as.numeric( "{params.treatment}" )
#          
#          
#          my.methylBase.obj = readRDS(input)
#          
#          # change treatment vector
#          my.methylBase.obj@treatment  = ifelse(my.methylBase.obj@treatment==t,  1, 0)
#          
#          myDiff<-calculateDiffMeth(my.methylBase.obj,
#                                   overdispersion="MN",
#                                   test="Chisq",
#                                   mc.cores=24)
#                                   
#          saveRDS(myDiff, output)    
#                  
#            """)
# 
# 

rule unite_meth_calls:
     input:
         [ expand(DIR_methcall+sample+"/"+sample+"_{chrom}_cpg_filtered.txt.bgz", chrom=CHROMS_CANON) for sample in SAMPLES]
     output:
         destrandTfile = expand(DIR_methcall+"methylBaseDB.obj_filtered_destrandT.{chrom}.RDS", chrom=CHROMS_CANON),
         destrandFfile = expand(DIR_methcall+"methylBaseDB.obj_filtered_destrandF.{chrom}.RDS", chrom=CHROMS_CANON)
     params:
         inputdir = DIR_methcall,
         samples = SAMPLES,
         treatments = TREATMENT,
         assembly=ASSEMBLY,
         cores=24
     log: expand(DIR_methcall+"meth_unite.{chrom}.log", chrom=CHROMS_CANON)
     shell:
       """
         {tools}/Rscript {DIR_scripts}/Unite_meth.R \
                 --inputfiles="{input}" \
                 --destrandTfile={output.destrandTfile} \
                 --destrandFfile={output.destrandFfile} \
                 --inputdir={params.inputdir} \
                 --samples="{params.samples}" \
                 --treatments="{params.treatments}" \
                 --assembly="{params.assembly}" \
                 --cores={params.cores} \
                 --logFile={log}
         """

rule filter_and_canonical_chroms:
     input:
         tabixfile     =  expand(DIR_methcall+"{{sample}}/{{sample}}_{chrom}_cpg.txt.bgz", chrom=CHROMS_CANON)
     output:
         outputfile    = expand(DIR_methcall+"{{sample}}/{{sample}}_{chrom}_cpg_filtered.txt.bgz", chrom=CHROMS_CANON)
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}/",
         sample_id = "{sample}",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=99,
         cores=10
     log:
         expand(DIR_methcall+"{{sample}}/{{sample}}_{chrom}.meth_calls_filter.log", chrom=CHROMS_CANON)
     message: ""
     shell:
       """
         {tools}/Rscript {DIR_scripts}/Filter_meth.R \
                 --tabixfile={input.tabixfile} \
                 --mincov={params.mincov} \
                 --hi_perc={params.hi_perc} \
                 --save_folder={params.save_folder} \
                 --sample_id={params.sample_id} \
                 --assembly={params.assembly} \
                 --cores={params.cores} \
                 --canon_chrs_file={params.canon_chrs_file} \
                 --logFile={log}
         """


rule methCall_CpG:
     input:
         bamfile = expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}.dedup.sorted.bam", chrom=CHROMS_CANON)
     output:
         callFile = expand(DIR_methcall+"{{sample}}/{{sample}}_{chrom}_cpg.txt.bgz", chrom=CHROMS_CANON)
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CpG",
         save_db      = True,
         save_folder = DIR_methcall+"{sample}/",
         sample_id = "{sample}"
     log:
         expand(DIR_methcall+"{{sample}}/{{sample}}.{chrom}.meth_calls.log", chrom=CHROMS_CANON)
     message: "Extract methylation calls from bam file."
     shell:
       """
          {tools}/Rscript {DIR_scripts}/methCall.R \
                 --inBam={input.bamfile} \
                 --assembly={params.assembly} \
                 --mincov={params.mincov} \
                 --minqual={params.minqual} \
                 --context={params.context} \
                 --save_db={params.save_db}  \
                 --save_folder={params.save_folder}  \
                 --logFile={log}
       """

rule sort_index_dedup:
     input:
         expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}.dedup.bam", chrom=CHROMS_CANON)
     output:
         expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}.dedup.sorted.bam", chrom=CHROMS_CANON)
     params:
         sort_args = config['args']['sambamba_sort'],
         tmpdir=DIR_deduped+"{sample}/"
     log:
         expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}sort.log", chrom=CHROMS_CANON)
     shell:
         "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"



rule deduplication:
     input:
         expand(DIR_mapped+"{{sample}}/{{sample}}_{chrom}.sorted.bam", chrom=CHROMS_CANON)
     output:
        outfile=expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}.dedup.bam", chrom=CHROMS_CANON),
        metrics = expand(DIR_deduped+"{{sample}}/{{sample}}_{chrom}.deduplication.metrics.txt", chrom=CHROMS_CANON)
     params:
         picard_MarkDuplicates_args = config['args']['picard_MarkDuplicates_args']
     log:
         expand(DIR_deduped+"{{sample}}/{{sample}}_deduplication.{chrom}", chrom=CHROMS_CANON)
     message:
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          """{tools}/picard MarkDuplicates I={input} O={output.outfile} \
          M={output.metrics} \
          REMOVE_DUPLICATES=true AS=true {params.picard_MarkDuplicates_args} \
          > {log}.log \
          2> {log}.err"""


rule split_bam_per_chr:
  input:  
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output: 
    expand(DIR_mapped+"{{sample}}/{{sample}}_{chrom}.sorted.bam", chrom=CHROMS_CANON)
  run:
    for chrom in CHROMS_CANON:
      mycmd = "%s/sambamba slice --output-filename %s %s %s" % (tools, output[0], input[0], chrom)
      shell(mycmd)
 
 
#      
# rule sort_index_bam_mapped:
#   input:
#     DIR_mapped+"{sample}/{sample}.bam"
#   output:
#     DIR_mapped+"{sample}/{sample}_sorted.bam"
#   params:
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+"{sample}/"
#   log:
#     DIR_mapped+"{sample}/{sample}_sort.log"
#   shell:
#     "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"

# 
# rule align_pe:
#      input:
#          refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
#          refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
#          fin1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
#          fin2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
#          #qc   = [ DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.html",
#          #        DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.html"]
#      output:
#          bam = DIR_mapped+"{sample}/{sample}.bam",
#          report = DIR_mapped+"{sample}/{sample}_report.txt",
#          un1 = DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz",
#          un2 = DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz",
#          odir = DIR_mapped+"{sample}/"
#      params:
#         # Bismark parameters
#          bismark_args = config['args']['bismark'],
#          genomeFolder = "--genome_folder " + genomedir,
#          outdir = "--output_dir  "+DIR_mapped+"{sample}/",
#          #nucCov = "--nucleotide_coverage",
#          pathToBowtie = "--path_to_bowtie " + config['tools'],
#          useBowtie2  = "--bowtie2 ",
#          samtools    = "--samtools_path "+ config['tools']+'samtools',
#          tempdir     = "--temp_dir "+DIR_mapped+"/{sample}"
#      log:
#          DIR_mapped+"{sample}/{sample}_bismark_pe_mapping.log"
#      message: "Mapping paired-end reads to genome."
#      run:
#          commands = [
# 	   #'{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err',
#          'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
#          'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
#          'ln -s '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
#          'ln -s '+output.odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
#          ]
#          for c in commands:
#             shell(c)
#           
# 
# 
