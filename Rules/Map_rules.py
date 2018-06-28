# ==========================================================================================
# 1A Mapping:

if SPLIT_INTO_PARTS == 1:
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
	       '{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
         'ln -s '+output.odir+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
         'ln -s '+output.odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
         'ln -s '+output.odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
         ]
         for c in commands:
            shell(c)
            
# ==========================================================================================
# 1B Mapping split input files:        

if SPLIT_INTO_PARTS > 1: 
  rule merge_bam_parts:
    input:
      [DIR_mapped+"{sample}/{sample}_part"+str(n)+"_sorted.bam" for n in range(1,SPLIT_INTO_PARTS+1)]
    output:
      DIR_mapped+"{sample}/{sample}_sorted.bam"
    params:
      sort_args = config['args']['sambamba_sort'],
      tmpdir=DIR_mapped+"{sample}/",
      unsorted = DIR_mapped+"{sample}/{sample}.bam"
    shell:
      """
      {tools}/sambamba merge {params.unsorted} {input}
      {tools}/sambamba sort {params.unsorted} --tmpdir={params.tmpdir} -o {output} {params.sort_args}
      """

rule sort_index_bam_mapped_parts:
  input:
    DIR_mapped+"{sample}/{sample}_part{n}.bam"
  output:
    DIR_mapped+"{sample}/{sample}_part{n}_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_part{n}_sort.log"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"


rule align_pe_parts:
     input:
         refconvert_CT = genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
         refconvert_GA = genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa",
         fin1 = DIR_trimmed+"{sample}/{sample}_1_part{n}_val_1.fq.gz",
         fin2 = DIR_trimmed+"{sample}/{sample}_2_part{n}_val_2.fq.gz",
         qc   = [ DIR_posttrim_QC+"{sample}/{sample}_1_part{n}_val_1_fastqc.html",
                 DIR_posttrim_QC+"{sample}/{sample}_2_part{n}_val_2_fastqc.html"]
     output:
         bam = DIR_mapped+"{sample}/{sample}_part{n}.bam",
         report = DIR_mapped+"{sample}/{sample}_part{n}_report.txt",
         un1 = DIR_mapped+"{sample}/{sample}_unmapped_1_part{n}.fq.gz",
         un2 = DIR_mapped+"{sample}/{sample}_unmapped_2_part{n}.fq.gz"
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
         DIR_mapped+"{sample}/{sample}_part{n}_bismark_pe_mapping.log"
     message: "Mapping paired-end reads to genome."
     run:
         commands = [
	       '{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err',
         'ln -s '+DIR_mapped+"{wildcards.sample}/"+os.path.basename(input.fin1[:-6])+'_bismark_bt2_pe.bam {output.bam}',
         'ln -s '+DIR_mapped+"{wildcards.sample}/"+os.path.basename(input.fin1[:-6])+'_bismark_bt2_PE_report.txt {output.report}',
         'ln -s '+DIR_mapped+"{wildcards.sample}/"+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1}',
         'ln -s '+DIR_mapped+"{wildcards.sample}/"+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2}'
         ]
         for c in commands:
            shell(c)


            
# ==========================================================================================
# Genome preparation:

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
            
