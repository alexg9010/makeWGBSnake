# ==========================================================================================
# Pre-mapping rules for single-end reads:


rule bismark_align_and_map_se:
    input:
        refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
	refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
        fin1 = DIR_trimmed+"{sample}/{sample}_trimmed.fq.gz"
        #qc     = DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.html"
    output:
         bam = DIR_mapped+"{sample}/{sample}.bam",
         report = DIR_mapped+"{sample}/{sample}_report.txt",
         un1 = DIR_mapped+"{sample}/{sample}_unmapped.fq.gz",
         odir = DIR_mapped+"{sample}/"
    params:
         # Bismark parameters
         bismark_args = config['args']['bismark'],
         genomeFolder = "--genome_folder " + genomedir,
         outdir = "--output_dir  "+DIR_mapped+"{sample}/",
         nucCov = "--nucleotide_coverage",
         pathToBowtie = "--path_to_bowtie " + config['tools'],
         useBowtie2  = "--bowtie2 ",
         samtools    = "--samtools_path "+ config['tools']+'samtools',
         tempdir     = "--temp_dir "+DIR_mapped+"/{sample}"      
    log:
        DIR_mapped+"{sample}/{sample}_bismark_pe_mapping.log"
    message: "Mapping single-end reads to genome {ASSEMBLY}"
    run:
         print(config['tools']+'bismark {params} {input.fin1} > {log} 2> {log}.err')
         commands = [
	       #config['tools']+'bismark {params} {input.fin1} > {log} 2> {log}.err',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2.bam {output.bam}',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_SE_report.txt {output.report}',
         'ln -s '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads.fq.gz {output.un1}'
         ]
         for c in commands:
            shell(c)


rule fastqc_after_trimming_se:
    input:
        DIR_trimmed+"{sample}/{sample}_trimmed.fq.gz",
    output:
    	DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.html",
    	DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.zip"
    params:
        fastqc_args = config['args']['fastqc'],
        outdir = "--outdir "+DIR_posttrim_QC
    log:
   	    DIR_posttrim_QC+"{sample}/{sample}_trimmed_fastqc.log"
    message: "Quality checking trimmmed single-end data from {input}"
    shell:
        config['tools']+"fastqc {params} {input} {log}"
        
        
rule trim_reads_se:
    input:
       qc   = DIR_rawqc+"{sample}/{sample}_fastqc.html",
       file = inputdir+"{sample}.fq.gz"
    output:
       DIR_trimmed+"{sample}/{sample}_trimmed.fq.gz" 
    params:
       extra          = config['args']['trim_galore'],
       outdir = "--output_dir "+ DIR_trimmed + "/{sample}",
       phred = "--phred33",
       gz = "--gzip",
       cutadapt = "--path_to_cutadapt " + config['tools']+'cutadapt'
    log:
       DIR_trimmed+"{sample}/{sample}.trimgalore.log"
    message: "Trimming raw single-end read data from {input}"
    shell:
       config['tools']+"trim_galore {params} {input.file} {log}"
       
      
