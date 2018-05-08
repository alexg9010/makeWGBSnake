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
config['samples'] = SAMPLES

CHROMS = [line.rstrip('\n') for line in open(genomedir+"chroms.txt")]
ASSEMBLY="hg19"
MINCOV=10
MINQUAL=20



WORKDIR = os.getcwd() + "/"                         
DIR_scripts   = './Scripts/'

DIR_plots = outputdir+'plots/'
DIR_bigwig      = outputdir+'bigwig_files/'
DIR_methcall    = outputdir+'methyl_calls/'
DIR_deduped     = outputdir+'deduplication/'
DIR_mapped      = outputdir+'bismark/'
DIR_posttrim_QC = outputdir+'posttrimming_QC/'
DIR_trimmed     = outputdir+'trimming/'
DIR_rawqc       = outputdir+'raw_QC/'
DIR_bam_per_chrom = outputdir+'methyl_calls/bam_per_chr/' 
DIR_seg = outputdir+'segmentation/'
#DIR_diffmeth    = output+'differential_methylation/'
DIR_ucsc_hub = outputdir+"ucsc_hub/"
DIR_multiqc = outputdir+"multiqc/"



# Construct all the files we're eventually expecting to have.
FINAL_FILES = []

# FASTQC
#FINAL_FILES.extend(
#   expand(DIR_rawqc+"{sample}_{ext}_fastqc.html",sample=config["samples"], ext=["1", "2"])
#)

# Alignment
FINAL_FILES.extend(
   expand(DIR_mapped+"{sample}/{sample}.bam",sample=config["samples"])
)

FINAL_FILES.extend(
    expand(DIR_mapped+"{sample}/{sample}_unmapped_{ext}_sorted.bam",sample=config["samples"], ext=["1", "2"])
)


# Sorting
FINAL_FILES.extend(
   expand(DIR_mapped+"{sample}/{sample}_sorted.bam",sample=config["samples"])
)

# Merge PE and SE reads
FINAL_FILES.extend(
   expand(DIR_mapped+"{sample}/{sample}_sorted_merged.bam", sample=config["samples"])
)

# Multiqc
#FINAL_FILES.extend(
#   expand(DIR_multiqc+"multiqc.html")
#)


# Deduplicate
#FINAL_FILES.extend(
#   expand(DIR_deduped+"{sample}/{sample}_sorted_dedup.bam",sample=config["samples"])
#)

# Split files
#FINAL_FILES.extend(
#   expand(DIR_bam_per_chrom+'{sample}/{sample}_sorted_dedup_{chrom}.bam',sample=config["samples"], chrom=['chr1'])
#)

#FINAL_FILES.extend(
#   expand(DIR_bam_per_chrom+'{sample}/{sample}_sorted_dedup_{chrom}.bam', sample=config["samples"],chrom=CHROMS)
#)

# Methylation calling
#FINAL_FILES.extend(
#   expand(DIR_bam_per_chrom+'{sample}/{sample}_sorted_dedup_{chrom}.bam', sample=config["samples"],chrom=CHROMS)
#)

# Merge methyl. calling files and create a tabix file
#FINAL_FILES.extend(
# expand(DIR_methcall+'{sample}/Tabix/{sample}_methyl.txt.bgz', sample=config["samples"])
#)


# Create BigWig files
#FINAL_FILES.extend(
# expand(DIR_bigwig+'{sample}/{sample}_{chrom}.bw', sample=config["samples"], chrom=CHROMS)
#)

print(FINAL_FILES)


rule target:
  input: FINAL_FILES



#ru leall:
#    input:
#        # QC
#        #expand(DIR_rawqc+"{sample}_fastqc.html",sample=config["samples"]),
#        # QC after trimming
#        #expand(DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",sample=config["samples"]),
#	# create genome index
#        #genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",



rule meth_segments:
     input:
         rdsfile     = os.path.join(DIR_methcall,"{prefix}.deduped_methylRaw.RDS") ####### TODO
     output:
         grfile      = os.path.join(DIR_seg,"{prefix}.deduped_meth_segments_gr.RDS"),
         bedfile     = os.path.join(DIR_seg,"{prefix}.deduped_meth_segments.bed")
     params:
         methCallRDS = os.path.join(WORKDIR,DIR_methcall,"{prefix}.deduped_methylRaw.RDS"),
         methSegGR       = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments_gr.RDS"),
         methSegBed      = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments.bed"),
         methSegPng = os.path.join(WORKDIR,DIR_seg,"{prefix}.deduped_meth_segments.png")
     log:
         os.path.join(DIR_seg,"{sample}/{sample}.deduped_meth_segments.log")
     message: "Segmenting methylation profile for {input.rdsfile}."
     shell:
         """
          Rscript {DIR_scripts}/methSeg.R
                          --rds={params.methCallRDS}
                          --grds={params.methSegGR}
                          --outBed={params.methSegBed}
                          --png={params.methSegPng}
                          --logFile={log}
          """                


rule export_bigwig:
     input:
         seqlengths = os.path.join(DIR_mapped,   "Refgen_"+ASSEMBLY+"_chromlengths.csv"),
         rdsfile    = [os.path.join(DIR_methcall, "{sample}/{sample}_sorted_dedup_"+chrom+"_methylRaw.RDS") for chrom in CHROMS]
     output:
         bw         = [os.path.join(DIR_bigwig, "{sample}_"+chrom+".bw") for chrom in CHROMS]
     message: "Exporting bigwig files"
     shell:
         """
          {tools}/Rscript {DIR_scripts}/export_bw.R \
                          {input.rdsfile} \
                          {input.seqlengths} \
                          ASSEMBLY \
                          {output} \
          """ 


# I dont use it currently
rule filterSNPs:
     input: 
           bam = expand(DIR_methcall+"{{sample}}/{{sample}}_sorted_dedup_{chrom}_methylRaw.RDS", chrom=CHROMS),
           bed = "SNPFile.bed",
           tabixfile = expand(DIR_methcall+"{{sample}}/Tabix/{{sample}}_{chrom}.txt.bgz", chrom=CHROMS)
     output:
            expand(DIR_methcall+"{{sample}}/{{sample}}.deduped_{chrom}_filteredSNP.RDS", chrom=CHROMS)
     message: "Filtering SNPs"
     shell:
            "{tools}/Rscript {DIR_scripts}/FilterSNP.R {input.bam} {input.bed} {output}"


# I dont use it currently
rule makeTabix:
     input:
         rdsfiles = [DIR_methcall+"{sample}/{sample}_sorted_dedup_"+chrom+"_methylRaw.RDS" for chrom in CHROMS],
         txtfiles = [DIR_methcall+"{sample}/{sample}_sorted_dedup_"+chrom+"_CpG.txt" for chrom in CHROMS]
     output:
         tabixdir = DIR_methcall+"{sample}/Tabix/",
         output = DIR_methcall+"{sample}/Tabix/{sample}_methyl.txt.bgz"
     #conda:
	# "Envs/env_min.yaml" 
     shell:
         """
         {tools}/Rscript {DIR_scripts}/makeTabix.R {output} {input}
         """  



rule bam_methCall:
     input:
         bamfile     = expand(DIR_bam_per_chrom+'{{sample}}/{{sample}}_sorted_dedup_{chrom}.bam', chrom=CHROMS)
     output:
         rdsfile     = expand(DIR_methcall+"{{sample}}/{{sample}}_sorted_dedup_{chrom}_methylRaw.RDS", chrom=CHROMS),
         callFile    = expand(DIR_methcall+"{{sample}}/{{sample}}_sorted_dedup_{chrom}_CpG.txt", chrom=CHROMS)
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CpG",
         #savedb      = False #TODO,
         savefolder = DIR_methcall+"{sample}/"
     log:
         os.path.join(DIR_methcall,"{sample}/{sample}.deduped_meth_calls.log")
     message: "Extract methylation calls from bam file."
     shell:
         """
         {tools}/Rscript {DIR_scripts}/methCall.R \
                 --inBam={input.bamfile} \
                 --assembly={params.assembly} \
                 --mincov={params.mincov} \
                 --minqual={params.minqual} \
                 --rds={output.rdsfile} \
                 --logFile={log}
         """


rule bam_sort_index_per_chr:
     input:
       expand(DIR_bam_per_chrom+'{{sample}}/{{sample}}_dedup_{chrom}.bam', chrom=CHROMS)
     output:
       expand(DIR_bam_per_chrom+'{{sample}}/{{sample}}_sorted_dedup_{chrom}.bam', chrom=CHROMS)
     params:
       sort_args = config['args']['sambamba_sort'],
       tmpdir=DIR_bam_per_chrom+"{sample}/"
     shell:
       "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}"



# https://github.com/daler/enhancer-snakemake-demo/blob/master/Snakefile
rule split_bam_per_chr:
     input:
        DIR_deduped+"{sample}/{sample}_sorted_dedup.bam"
     output:
        expand(DIR_bam_per_chrom+'{{sample}}/{{sample}}_dedup_{chrom}.bam', chrom=CHROMS)
     run:
        sampleid=list(wildcards)[0]
        for chrom in CHROMS:
          cmd=tools+"/sambamba slice "+input[0]+" " +chrom+" -o "+DIR_bam_per_chrom+sampleid+"/"+sampleid+"_dedup_"+chrom+".bam"
          #print(cmd)
          shell(cmd)



rule sort_index_dedup:
     input:
         DIR_deduped+"{sample}/{sample}_dedup.bam"
     output:
         DIR_deduped+"{sample}/{sample}_sorted_dedup.bam"
     params:
         sort_args = config['args']['sambamba_sort'],
         tmpdir=DIR_deduped+"{sample}/"
     log:
         DIR_deduped+"{sample}/{sample}_sort.log"
     shell:
         "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"



rule deduplication:
     input:
         DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
     output:
         DIR_deduped+"{sample}/{sample}_dedup.bam"
     params:
         metrics=DIR_deduped+"{sample}/{sample}_dup_metrics.txt"
     log:
         DIR_deduped+"{sample}/{sample}_deduplication.log"
     message: 
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          "{tools}/picard MarkDuplicates I={input} O={output} M={params.metrics} REMOVE_DUPLICATES=true AS=true > {log} 2> {log}.err"


#rule multiqc:
#    input:
#        DIR_deduped+"{sample}/{sample}_dup_metrics.txt",
#        DIR_trimmed+"{sample}_1.fastq.gz_trimming_report.txt",
#        DIR_trimmed+"{sample}_2.fastq.gz_trimming_report.txt",
#        DIR_posttrim_QC+"{sample}_1_val_1_fastqc.zip",
#        DIR_posttrim_QC+"{sample}_2_val_2_fastqc.zip"
#    output:
#        "{DIR_multiqc}multiqc.html"
#    params:
#        ""  # Optional: extra parameters for multiqc.
#    log:
#        "{DIR_multiqc}multiqc.log"
#    wrapper:
#        "0.23.1/bio/multiqc"



rule merge_bam_se_and_pe:
  input:
    DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam",
    DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam",
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
  shell:
    "{tools}/sambamba merge {output} {input}"



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
  


rule align_sort_index_unmapped2_se:
    input:
        DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
    output:
        outfile=DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam",
        outdir=DIR_mapped+"{sample}/"
    params:
        bismark_args = config['args']['bismark_unmapped'],
        genomeFolder = "--genome_folder " + genomedir,
        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
        pathToBowtie = "--path_to_bowtie " + config['tools'],
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ config['tools']+"samtools",
        tmpdir     = "--temp_dir "+DIR_mapped+"{sample}/",
        sort_args = config['args']['sambamba_sort']
    log:
        align=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2.log",
        sort=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2_sort.log" 
    message: "Mapping paired-end reads as single-end to genome."
    shell:
        """
        {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
        {tools}/sambamba sort {output.outdir}{wildcards.sample}_unmapped_2_bismark_bt2.bam --tmpdir={params.tmpdir} -o {output.outfile} {params.sort_args} > {log.sort} 2> {log.sort}.err
        """


rule align_sort_index_unmapped1_se:
    input:
        DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz"
    output:
        outfile=DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam",
        outdir=DIR_mapped+"{sample}/"
    params:
        bismark_args = config['args']['bismark_unmapped'],
        genomeFolder = "--genome_folder " + genomedir,
        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
        pathToBowtie = "--path_to_bowtie " + config['tools'],
        useBowtie2  = "--bowtie2 ",
        samtools    = "--samtools_path "+ config['tools']+"samtools",
        tmpdir     = "--temp_dir "+DIR_mapped+"{sample}/",
        sort_args = config['args']['sambamba_sort']
    log:
        align=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_1.log",
        sort=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_1_sort.log"
    message: "Mapping paired-end reads as single-end to genome."
    shell:
        """
        {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
        {tools}/sambamba sort {output.outdir}{wildcards.sample}_unmapped_1_bismark_bt2.bam --tmpdir={params.tmpdir} -o {output.outfile} {params.sort_args} > {log.sort} 2> {log.sort}.err
        """


#def input_bismark_align_unmapped_pe_as_se(wildcards):
#    if (wildcards.ext=="1"):
#       input=DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz"
#    elif (wildcards.ext=="2"):
#       input=DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#    return(input)

#rule bismark_align_unmapped_pe_as_se:
#    input:
#        DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#    output:
#        out = DIR_mapped+"{sample}/{sample}_unmapped_{ext}.bam",
#        odir = DIR_mapped+"{sample}/"
#    params:
#        # intermediate files
#        inter_bam = DIR_mapped+"{sample}_unmapped_2_bismark_bt2.bam",
#        # bismark parameters
#        bismark_args = config['args']['bismark_unmapped'],
#        genomeFolder = "--genome_folder " + genomedir,
#        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
#        pathToBowtie = "--path_to_bowtie " + config['tools'],
#        useBowtie2  = "--bowtie2 ",
#        samtools    = "--samtools_path "+ config['tools']+"samtools",
#        tempdir     = "--temp_dir "+DIR_mapped+"{sample}/",
#    log:
#        DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2.log"
#    message: "Mapping paired-end reads as single-end to genome."
#    run:
#        cmds=[  
#        '{tools}/bismark {params} {input.fin2} > {log} 2> {log}.err'
#        'ln -s ' + output.odir + os.path.basename(input.fin2)+'_bismark_bt2_pe.bam'+ + ' {output.out2}'
#        ]
#        for c in cmds:
#           shell(c)

        


rule align_pe:
     input:
         refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
 	 refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
         fin1 = DIR_trimmed+"{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed+"{sample}_2_val_2.fq.gz",
         qc   = [ DIR_posttrim_QC+"{sample}_1_val_1_fastqc.html",
                  DIR_posttrim_QC+"{sample}_2_val_2_fastqc.html"]
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
         outdir = "--output_dir  "+DIR_mapped+"/{sample}",
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
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
         'ln -s '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
         'ln -s '+output.odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
         ]
         for c in commands:
            shell(c)
         
 
 
#rule bismark_genome_preparation:
#     input:
#         ancient(genomedir)
#     output:
#         genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
#         genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
#     params:
#         bismark_genome_preparation_args = config['args']['bismark_genome_preparation'],
#         pathToBowtie = "--path_to_bowtie "+ config['tools'],
#         useBowtie2 = "--bowtie2 ",
#         verbose = "--verbose "
#     log:
#         outputdir+'bismark_genome_preparation.log'
     #message: 
     #	 "Converting Genome {input} into Bisulfite analogue"
#     shell:
#         "{tools}/bismark_genome_preparation {params} {input} > {log} 2> {log}.err"
 

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
         "{tools}/fastqc {params} {input} > {log} 2> {log}.err"

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
         cutadapt       = "--path_to_cutadapt " + tools +"cutadapt",
         paired         = "--paired"
     log:
         DIR_trimmed+"{sample}.trimgalore.log"
     message:
         "Trimming raw paired-end read data from {input}"
     shell:
       "{tools}/trim_galore {params} {input.files} > {log} 2> {log}.err"
   

rule fastqc_raw:
    input:
       inputdir+"{sample}.fq.gz"
    output:
        DIR_rawqc+"{sample}_fastqc.html",
        DIR_rawqc+"{sample}_fastqc.zip"
    params:
       fastqc_args = config['args']['fastqc'],
       outdir = "--outdir "+ DIR_rawqc
    log:
        DIR_rawqc+"{sample}_fastqc.log"
    shell:
        config["tools"]+"fastqc {params} {input} > {log} 2> {log}.err"
        

#rule fastqc_raw:
#    input:
#        inputdir+"{sample}.fq.gz"
#    output:
#        html="{DIR_rawqc}{sample}.html",
#       zip="{DIR_rawqc}{sample}.zip"
#    params: config['args']['fastqc']
#    wrapper:
#        "0.23.1/bio/fastqc"


