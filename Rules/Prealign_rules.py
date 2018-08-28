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


rule fastqc_raw_pe:
    input:
       inputdir+"{sample}_{ext}.fastq.gz"
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
       
       
# rule split_reads_2_parts:
#      input:
#          files = [ inputdir+"{sample}_1.fq.gz",
#                    inputdir+"{sample}_2.fq.gz"]
#      output:
#          DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz", 
#          DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz",
#      params:
#          extra          = config['args']['trim_galore'],
#          outdir         = "--output_dir "+DIR_trimmed+"{sample}/",
#          phred          = "--phred33",
#          gz             = "--gzip",
#          cutadapt       = "--path_to_cutadapt " + tools +"cutadapt",
#          paired         = "--paired"
#      log:
#          DIR_trimmed+"{sample}/{sample}.log"
#      message:
#          "Splitting paired-end reads {input} into parts"
#      shell:
#        "{tools}/Rscript {DIR}Split_fastqfile.R {params} {input.files} > {log} 2> {log}.err"
#           
        
