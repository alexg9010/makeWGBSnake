# ==========================================================================================
# Pre-mapping rules for paired-end reads:

rule fastqc_after_trimming_pe:
     input:
         DIR_trimmed+"{sample}/{sample}_1_part{n}_val_1.fq.gz",
         DIR_trimmed+"{sample}/{sample}_2_part{n}_val_2.fq.gz"
     output:
         DIR_posttrim_QC+"{sample}/{sample}_1_part{n}_val_1_fastqc.html",
         DIR_posttrim_QC+"{sample}/{sample}_1_part{n}_val_1_fastqc.zip",
         DIR_posttrim_QC+"{sample}/{sample}_2_part{n}_val_2_fastqc.zip",
         DIR_posttrim_QC+"{sample}/{sample}_2_part{n}_val_2_fastqc.html"
     params:
         fastqc_args = config['args']['fastqc'],
         outdir = "--outdir "+DIR_posttrim_QC + "{sample}/"
     log:
    	   DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.log"
     message:
         "Quality checking trimmmed paired-end data from {input}"
     shell:
         "{tools}/fastqc {params} {input} > {log} 2> {log}.err"


rule trim_reads_pe_part:
     input:
         files = [ DIR_inputdir_parts+"{sample}_1_part{n}.fq.gz",
                   DIR_inputdir_parts+"{sample}_2_part{n}.fq.gz"]
     output:
         DIR_trimmed+"{sample}/{sample}_1_part{n}_val_1.fq.gz", 
         DIR_trimmed+"{sample}/{sample}_2_part{n}_val_2.fq.gz"
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
       
rule split_PE_reads_into_parts:
     input:
         in1=  inputdir+"{sample}_1.fq.gz",
         in2 = inputdir+"{sample}_2.fq.gz"
     output:
         out1 = [DIR_inputdir_parts+"{sample}_1_part"+str(n)+".fq.gz" for n in range(1,SPLIT_INTO_PARTS+1)],
         out2 = [DIR_inputdir_parts+"{sample}_2_part"+str(n)+".fq.gz" for n in range(1,SPLIT_INTO_PARTS+1)]
     params:
         parts = 5,
         cores = 2,
         outdir = DIR_inputdir_parts
     message:
         "Splitting paired-end reads {input} into parts"
     shell:
       "{tools}/Rscript {DIR_scripts}/Split_fastqfiles.R  {input.in1} {input.in2} {params.outdir} {params.parts}"
          
        
rule fastqc_raw_pe:
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
       
      
