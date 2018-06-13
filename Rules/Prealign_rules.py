# ==========================================================================================
# Pre-mapping rules:

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
    message:
        "Converting Genome {input} into Bisulfite analogue"
    shell:
        "{tools}/bismark_genome_preparation {params} {input} > {log} 2> {log}.err"


rule fastqc_after_trimming_pe:
     input:
         DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
         DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz"
     output:
     	DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.html",
     	DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.zip",
     	DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.zip",
         DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.html"
     params:
         fastqc_args = config['args']['fastqc'],
         outdir = "--outdir "+DIR_posttrim_QC + "{sample}/"
     log:
    	    DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.log"
     message:
       "Quality checking trimmmed paired-end data from {input}"
     shell:
         "{tools}/fastqc {params} {input} > {log} 2> {log}.err"

rule trim_reads_pe:
     input:
         files = [ inputdir+"{sample}_1.fq.gz",
                   inputdir+"{sample}_2.fq.gz"]
     output:
         DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz", 
         DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
     params:
         extra          = config['args']['trim_galore'],
         outdir         = "--output_dir "+DIR_trimmed+"{sample}/",
         phred          = "--phred33",
         gz             = "--gzip",
         cutadapt       = "--path_to_cutadapt " + tools +"cutadapt",
         paired         = "--paired"
     log:
         DIR_trimmed+"{sample}/{sample}.log"
     message:
         "Trimming raw paired-end read data from {input}"
     shell:
       "{tools}/trim_galore {params} {input.files} > {log} 2> {log}.err"
   

rule fastqc_raw:
    input:
       inputdir+"{sample}_{ext}.fq.gz"
    output:
        DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",
        DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.zip"
    params:
       fastqc_args = config['args']['fastqc'],
       outdir = "--outdir "+ DIR_rawqc +"{sample}/"
    log:
        DIR_rawqc+"{sample}/{sample}_{ext}.log"
    shell:
        "{tools}/fastqc {params} {input} > {log} 2> {log}.err"
        
