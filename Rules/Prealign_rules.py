# ==========================================================================================
# Pre-mapping rules for paired-end reads:

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
     benchmark:
        outputdir+"benchmarks/{sample}.fastqc_after_trimming_pe.benchmark.txt"
     shell:
         "{tools}/fastqc {params} {input} > {log} 2> {log}.err"

rule trim_reads_pe:
     input:
         files = [ inputdirall+"{sample}_1.fq.gz", ###########
                   inputdirall+"{sample}_2.fq.gz"]############
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
     benchmark:
        outputdir+"benchmarks/{sample}.trim_reads_pe.benchmark.txt"
     shell:
       "{tools}/trim_galore {params} {input.files} > {log} 2> {log}.err"


rule fastqc_raw_pe:
    input:
       inputdirall+"{sample}_{ext}.fq.gz"
    output:
        DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.html",
        DIR_rawqc+"{sample}/{sample}_{ext}_fastqc.zip"
    params:
       fastqc_args = config['args']['fastqc'],
       outdir = "--outdir "+ DIR_rawqc +"{sample}/"
    log:
        DIR_rawqc+"{sample}/{sample}_{ext}.log"
    #benchmark:
    #    outputdir+"benchmarks/{sample}.fastqc_raw_pe.benchmark.txt"
    shell:
        "{tools}/fastqc {params} {input} > {log} 2> {log}.err"
       
       

