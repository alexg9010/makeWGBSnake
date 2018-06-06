#  
# ################################ Process UNALIGNED START
# 
# rule merge_bam_se_and_pe:
#   input:
#     DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam",
#     DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam",
#     DIR_mapped+"{sample}/{sample}_sorted.bam"
#   output:
#     DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
#   shell:
#     "{tools}/sambamba merge {output} {input}"
# 
# 
# 
# rule sort_index_bam_unmapped2_se:
#   input:
#     DIR_mapped+"{sample}/{sample}_unmapped_2.bam"
#   output:
#     DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam"
#   params:
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+"{sample}/"
#   log:
#     DIR_mapped+"{sample}/{sample}_sort2.log"
#   shell:
#     "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"
# 
# 
# rule align_unmapped2_se:
#     input:
#         DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#     output:
#         outfile=DIR_mapped+"{sample}/{sample}_unmapped_2.bam",
#         outdir=DIR_mapped+"{sample}/"
#     params:
#         bismark_args = config['args']['bismark_unmapped'],
#         genomeFolder = "--genome_folder " + genomedir,
#         outdir = "--output_dir  "+DIR_mapped+"{sample}/",
#         pathToBowtie = "--path_to_bowtie " + config['tools'],
#         useBowtie2  = "--bowtie2 ",
#         samtools    = "--samtools_path "+ config['tools']+"samtools",
#         tmpdir     = "--temp_dir "+DIR_mapped+"{sample}/",
#         sort_args = config['args']['sambamba_sort']
#     log:
#         align=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2.log",
#         sort=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2_sort.log"
#     message: "Mapping paired-end reads as single-end to genome."
#     shell:
#         """
#         {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
#         ln -s {output.outdir}{wildcards.sample}_unmapped_2_bismark_bt2.bam {output.outfile}
#         """
# 
# 
# rule sort_index_bam_unmapped1_se:
#   input:
#     DIR_mapped+"{sample}/{sample}_unmapped_1.bam"
#   output:
#     DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam"
#   params:
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+"{sample}/"
#   log:
#     DIR_mapped+"{sample}/{sample}_sort1.log"
#   shell:
#     "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"
# 
# 
# rule align_unmapped1_se:
#     input:
#         DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz"
#     output:
#         outfile=DIR_mapped+"{sample}/{sample}_unmapped_1.bam",
#         outdir=DIR_mapped+"{sample}/"
#     params:
#         bismark_args = config['args']['bismark_unmapped'],
#         genomeFolder = "--genome_folder " + genomedir,
#         outdir = "--output_dir  "+DIR_mapped+"{sample}/",
#         pathToBowtie = "--path_to_bowtie " + config['tools'],
#         useBowtie2  = "--bowtie2 ",
#         samtools    = "--samtools_path "+ config['tools']+"samtools",
#         tmpdir     = "--temp_dir "+DIR_mapped+"{sample}/",
#         sort_args = config['args']['sambamba_sort']
#     log:
#         align=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_1.log",
#         sort=DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_1_sort.log"
#     message: "Mapping paired-end reads as single-end to genome."
#     shell:
#         """
#         {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
#         ln -s {output.outdir}{wildcards.sample}_unmapped_1_bismark_bt2.bam {output.outfile}
#         """
# 
# ################################ Process UNALIGNED END
