# WGBS pipeline
#
# Copyright © 2018 Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
#

import glob, os, re

inputdir = config["input"]
outputdir = config["output"]
genomedir = config["genome"]
envs = config["env"]

config["samples"] = [re.sub('\\_1.fq.gz$', '', os.path.basename(x)) for x in glob.glob(inputdir+"*_1.fq.gz")]

print(config)

WORKDIR = os.getcwd() + "/"                         

DIR_scripts   = './Scripts/'

#DIR_diffmeth    = output+'differential_methylation/'
#DIR_seg         = output+'segmentation/'
DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'bigwig_files/'
DIR_methcall    = outputdir+'methyl_calls/'
DIR_deduped     = outputdir+'deduplication/'
DIR_mapped      = outputdir+'bismark/'
DIR_posttrim_QC = outputdir+'posttrimming_QC/'
DIR_trimmed     = outputdir+'trimming/'
DIR_rawqc       = outputdir+'raw_QC/'


rule all:
    input:
        #expand(DIR_mapped+"{sample}_1_val_1_bismark_bt2_pe.bam",dataset=config["samples"])
        #expand(DIR_rawqc+"{sample}_fastqc.html",sample=config["samples"])
        # QC after trimming
        expand(DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",sample=config["samples"]),
	# create genome index
        genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
        # map
        expand(DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe.bam",sample=config["samples"]),
        expand(DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam",sample=config["samples"]),
        expand(DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam",sample=config["samples"]),
	# sort
	expand(DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2_sorted.bam",sample=config["samples"]),
	expand(DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2_sorted.bam",sample=config["samples"]),
	expand(DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe_sorted.bam",sample=config["samples"]),
	expand(DIR_mapped+"{sample}/{sample}_merged.bam",sample=config["samples"])

print(config)

# rule methseg:
#     ## paths inside input and output should be relative
#     input:
#         rdsfile     = os.path.join(DIR_methcall,"{prefix}.deduped_methylRaw.RDS")
#     output:
#         grfile      = os.path.join(DIR_seg,"{prefix}.deduped_meth_segments_gr.RDS"),
#         bedfile     = os.path.join(DIR_seg,"{prefix}.deduped_meth_segments.bed")
#     params:
#         methCallRDS = os.path.join(WORKDIR,DIR_methcall,"{prefix}.deduped_methylRaw.RDS"),
#         methSegGR       = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments_gr.RDS"),
#         methSegBed      = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments.bed"),
#         methSegPng      = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments.png")
#     log:
#         os.path.join(DIR_seg,"{prefix}.deduped_meth_segments.log")
#     message: fmt("Segmenting methylation profile for {input.rdsfile}.")
#     shell:
#         nice('Rscript', ["{DIR_scripts}/methSeg.R",
#                          "--rds={params.methCallRDS}",
#                          "--grds={params.methSegGR}",
#                          "--outBed={params.methSegBed}",
#                          "--png={params.methSegPng}",
#                          "--logFile={log}"])
#                          
# rule export_bigwig_pe:
#     input:
#         seqlengths = os.path.join(DIR_mapped,   "Refgen_"+ASSEMBLY+"_chromlengths.csv"),
#         rdsfile    = os.path.join(DIR_methcall, "{prefix}_1_val_1_bt2.sorted.deduped_methylRaw.RDS")
#     output:
#         bw         = os.path.join(DIR_bigwig,   "{prefix}_pe.bw")
#     message: fmt("exporting bigwig files from paired-end stream.")
#     shell:
#         nice('Rscript', ["{DIR_scripts}/export_bw.R",
#                          "{input.rdsfile}",
#                          "{input.seqlengths}",
#                          ASSEMBLY,
#                          "{output}"])
# 
# rule bam_methCall:
#     input:
#         bamfile     = os.path.join(DIR_deduped,"{prefix}.deduped.bam")
#     output:
#         rdsfile     = os.path.join(DIR_methcall,"{prefix}.deduped_methylRaw.RDS"),
#         callFile    = os.path.join(DIR_methcall,"{prefix}.deduped_CpG.txt")
#     params:
#         ## absolute path to bamfiles
#         inBam       = os.path.join(WORKDIR,DIR_deduped,"{prefix}.deduped.bam"),
#         assembly    = ASSEMBLY,
#         mincov      = int(config['general']['methylation-calling']['minimum-coverage']),
#         minqual     = int(config['general']['methylation-calling']['minimum-quality']),
#         ## absolute path to output folder in working dir
#         rds         = os.path.join(WORKDIR,DIR_methcall,"{prefix}.deduped_methylRaw.RDS")
#     log:
#         os.path.join(DIR_methcall,"{prefix}.deduped_meth_calls.log")
#     message: fmt("Extract methylation calls from bam file.")
#     shell:
#         nice('Rscript', ["{DIR_scripts}/methCall.R",
#                          "--inBam={params.inBam}",
#                          "--assembly={params.assembly}",
#                          "--mincov={params.mincov}",
#                          "--minqual={params.minqual}",
#                          "--rds={params.rds}",
#                          "--logFile={log}"])
# 
# rule split_bam_per_chr:
# 
# 
rule deduplication:
     input:
         DIR_mapped+"{sample}/{sample}_merged_sorted.bam"
     output:
         DIR_deduped+"{sample}/{sample}_merged_sorted.deduped.bam"
     params:
         metrics=DIR_deduped+"{sample}/{sample}_dup_metrics.txt"
     log:
         DIR_deduped+"{sample}/{sample}_deduplication.log"
     message: 
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          config['tools']+"picard MarkDuplicates I={input} O={output} M={params.metrics} REMOVE_DUPLICATES=true AS=true > {log} 2> {log}.err"
 
 

#rule split_per_chr:
#     input:
#        DIR_mapped+"{sample}/{sample}_merged_sorted.bam"
#     output:
#	[chrom for chrom in ]


rule merge:
    input:
        DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2_sorted.bam",
        DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2_sorted.bam",
        DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe_sorted.bam"
    output:
        DIR_mapped+"{sample}/{sample}_merged.bam"
    shell:
        config['tools']+"sambamba merge {output} {input}"


rule sort_pe_as_se2:
  input:
    DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam"
  output: 
    DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/",
    #t=1
  log:
    DIR_mapped+"{sample}/{sample}_sort.log"
  shell:
    config['tools']+"sambamba sort {input} {params} -o {output} > {log} 2> {log}.err"

 
rule sort_pe_as_se1:
  input:
    DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam"
  output: 
    DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/",
    #t=1
  log:
    DIR_mapped+"{sample}/{sample}_sort.log"
  shell:
    config['tools']+"sambamba sort {input} {params} -o {output} > {log} 2> {log}.err"
 
rule sort:
  input: 
    DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe.bam"
  output: 
    DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/",
    #t=1
  log:
    DIR_mapped+"{sample}/{sample}_sort.log"
  shell:
    config['tools']+"sambamba sort {input} {params} -o {output}  > {log} 2> {log}.err"
  

rule bismark_align_unmappedpe_as_se1:
    input:
        fin1 = DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1.fq.gz"
    output:
        out1 = DIR_mapped+"{sample}/{sample}_1_val_1.fq.gz_unmapped_reads_1_bismark_bt2.bam"
    params:
        bismark_args = config['args']['bismark'],
        genomeFolder = "--genome_folder " + genomedir,
        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
        pathToBowtie = "--path_to_bowtie " + config['tools'],
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ config['tools']+"samtools",
        tempdir     = "--temp_dir "+DIR_mapped+"{sample}/",
    log:
        DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_1.log"
    message: "Mapping paired-end reads as single-end to genome."
    shell:
        config['tools']+"bismark {params} {input.fin1} > {log} 2> {log}.err"


rule bismark_align_unmappedpe_as_se2:
    input:
        fin2 = DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2.fq.gz"
    output:
        out2 = DIR_mapped+"{sample}/{sample}_2_val_2.fq.gz_unmapped_reads_2_bismark_bt2.bam"
    params:
        bismark_args = config['args']['bismark'],
        genomeFolder = "--genome_folder " + genomedir,
        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
        pathToBowtie = "--path_to_bowtie " + config['tools'],
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ config['tools']+"samtools",
        tempdir     = "--temp_dir "+DIR_mapped+"{sample}/",
    log:
        DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2.log"
    message: "Mapping paired-end reads as single-end to genome."
    shell:
        config['tools']+"bismark {params} {input.fin2} > {log} 2> {log}.err"



rule bismark:
     input:
         refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
 	 refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
         fin1 = DIR_trimmed+"{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed+"{sample}_2_val_2.fq.gz",
         qc   = [ DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
                  DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"]
     output:
         DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_pe.bam",
         DIR_mapped+"{sample}/{sample}_1_val_1_bismark_bt2_PE_report.txt"
     params:
         bismark_args = config['args']['bismark'],
         genomeFolder = "--genome_folder " + genomedir,
         outdir = "--output_dir  "+DIR_mapped+"/{sample}",
         #nucCov = "--nucleotide_coverage",
         pathToBowtie = "--path_to_bowtie " + config['tools'],
         useBowtie2  = "--bowtie2 ",
         samtools    = "--samtools_path "+ config['tools']+'samtools',
         tempdir     = "--temp_dir "+DIR_mapped+"/{sample}"
     log:
         DIR_mapped+"{sample}/{sample}_bismark_pe_mapping.log"
     message: "Mapping paired-end reads to genome."
     shell:
         config['tools']+"bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err"
 
 
rule bismark_genome_preparation:
     input:
         ancient(genomedir)
     output:
         genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
     params:
         bismark_genome_preparation_args = config['args']['bismark_genome_preparation'],
         pathToBowtie = "--path_to_bowtie "+ config['tools'],
         useBowtie2 = "--bowtie2 ",
         verbose = "--verbose "
     log:
         outputdir+'bismark_genome_preparation.log'
     #message: 
     #	 "Converting Genome {input} into Bisulfite analogue"
     shell:
         "bismark_genome_preparation {params} {input} > {log} 2> {log}.err"
 

rule fastqc_after_trimming_pe:
     input:
         DIR_trimmed+"{sample}_1_val_1.fq.gz",
         DIR_trimmed+"{sample}_2_val_2.fq.gz"
     output:
     	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
     	DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
     	DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip",
         DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"
     params:
         fastqc_args = config['args']['fastqc'],
         outdir = "--outdir "+DIR_posttrim_QC
     log:
    	    DIR_posttrim_QC+"{sample}_trimmed_fastqc.log"
     message:
       "Quality checking trimmmed paired-end data from {input}"
     shell:
         config['tools']+"fastqc {params} {input} > {log} 2> {log}.err"

rule trim_reads_pe:
     input:
         #qc    = [ DIR_rawqc+"{sample}_1_fastqc.html",
         #          DIR_rawqc+"{sample}_2_fastqc.html"],
         files = [ inputdir+"{sample}_1.fq.gz",
                   inputdir+"{sample}_2.fq.gz"]
     output:
         DIR_trimmed+"{sample}_1_val_1.fq.gz", 
         DIR_trimmed+"{sample}_2_val_2.fq.gz",
     params:
         extra          = config['args']['trim_galore'],
         outdir         = "--output_dir "+DIR_trimmed,
         phred          = "--phred33",
         gz             = "--gzip",
         cutadapt       = "--path_to_cutadapt " + config["tools"]+"cutadapt",
         paired         = "--paired"
     log:
         DIR_trimmed+"{sample}.trimgalore.log"
     message:
         "Trimming raw paired-end read data from {input}"
     shell:
       config["tools"]+"trim_galore {params} {input.files} > {log} 2> {log}.err"

#rule fastqc_raw:
#    input:
#        inputdir+"{sample}.fq.gz"
#    output:
#        DIR_rawqc+"{sample}_fastqc.html",
#        DIR_rawqc+"{sample}_fastqc.zip"
#    params:
#        fastqc_args = config['args']['fastqc'],
#        outdir = "--outdir "+ DIR_rawqc
#    log:
#        DIR_rawqc+"{sample}_fastqc.log"
#    shell:
#        config["tools"]+"fastqc {params} {input} > {log} 2> {log}.err"
        

