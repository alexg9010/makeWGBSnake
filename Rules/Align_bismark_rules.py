

# # ==========================================================================================
# # Mapping:
#


rule idxstats_bismark:
  input:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.idxstats.txt"
  shell:
    "samtools idxstats {input} > {output}"


rule stat_bismark:
  input:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.stats.txt"
  shell:
    "samtools stats {input} > {output}"


rule flagstat_bismark:
  input:
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.flagstat.txt"
  shell:
    "samtools flagstat {input} > {output}"



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
  benchmark:
    outputdir+"benchmarks/{sample}.sort_index_bam_mapped.benchmark.txt"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"


rule bismark_align_pe:
     input:
         refconvert_CT = ancient(genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa"),
         refconvert_GA = ancient(genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"),
         fin1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz" if not SUBSET_READS else DIR_trimmed_subset+"{sample}_1.fq.gz",############
         fin2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz" if not SUBSET_READS else DIR_trimmed_subset+"{sample}_2.fq.gz",############
         #qc   = [ DIR_posttrim_QC+"{sample}/{sample}_1_val_1_fastqc.html",
         #        DIR_posttrim_QC+"{sample}/{sample}_2_val_2_fastqc.html"]
     output:
         bam = DIR_mapped+"{sample}/{sample}.bam",
         report = DIR_mapped+"{sample}/{sample}_bismark_bt2_PE_report.txt",
         un1 = DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz",
         un2 = DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
         #odir = DIR_mapped+"{sample}/"
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
     benchmark:
        outputdir+"benchmarks/{sample}.align_pe.benchmark.txt"
     run:
         sample_name = os.path.basename(input.fin1[:-8]) # remove suffix '_1.fq.gz'
         output_odir = DIR_mapped+sample_name+"/"
         commands = [
	       '{tools}/bismark {params} -1 {input.fin1} -2 {input.fin2} > {log} 2> {log}.err;',
         'mv '+output_odir+sample_name+'_1_val_1_bismark_bt2_pe.bam {output.bam};' if not SUBSET_READS else 'mv '+output_odir+sample_name+'_1_bismark_bt2_pe.bam {output.bam};',
         'mv '+output_odir+sample_name+'_1_val_1_bismark_bt2_PE_report.txt {output.report};' if not SUBSET_READS else 'mv '+output_odir+sample_name+'_1_bismark_bt2_PE_report.txt {output.report};',
         'mv '+output_odir+os.path.basename(input.fin1)+'_unmapped_reads_1.fq.gz {output.un1};',
         'mv '+output_odir+os.path.basename(input.fin2)+'_unmapped_reads_2.fq.gz {output.un2};'
         ]
         for c in commands:
           shell(c)


      
# 
# # ==========================================================================================
# # Generate methyl-converted version of the reference genome:
#            
# 
# rule bismark_genome_preparation:
#     input:
#         ancient(genomedir)
#     output:
#         genomedir+"Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa",
#         genomedir+"Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa"
#     params:
#         pathToBowtie = "--path_to_bowtie "+ tools,
#         useBowtie2   = "--bowtie2 ",
#         verbose      = "--verbose "
#     log:
#         genomedir+'bismark_genome_preparation_'+ASSEMBLY+'.log'
#     message: "Converting {ASSEMBLY} Genome into Bisulfite analogue"
#     shell:
#         "bismark_genome_preparation {params} {input} > {log} 2> {log}.err"

