

# ==========================================================================================
# Deduplication:
# 


rule sort_index_dedup_pe_se:
     input:
         DIR_deduped_picard+"{sample}/{sample}_merged.dedup.bam"
     output:
         DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam"
     params:
         sort_args = config['args']['sambamba_sort'],
         tmpdir=DIR_deduped_picard+"{sample}/"
     benchmark:
          outputdir+"benchmarks/{sample}.sort_index_dedup_pe_se.benchmark.txt"
     shell:
         "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}"

rule dedup_picard_pe_se:
     input:
         DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
     output:
        outfile=temp(DIR_deduped_picard+"{sample}/{sample}_merged.dedup.bam"),
        metrics = DIR_deduped_picard+"{sample}/{sample}_merged.deduplication.metrics.txt"
     params:
         picard_MarkDuplicates_args = config['args']['picard_MarkDuplicates_args']
     log:
         DIR_deduped_picard+"{sample}/{sample}_merged_deduplication.log"
     message:
          "Deduplicating paired-end aligned reads from {input}"
     benchmark:
          outputdir+"benchmarks/{sample}.dedup_picard_pe_se.benchmark.txt"
     shell:
          """{tools}/picard MarkDuplicates I={input} O={output.outfile} \
          M={output.metrics} \
          REMOVE_DUPLICATES=true AS=true {params.picard_MarkDuplicates_args} \
          > {log} \
          2> {log}.err"""


# # ==========================================================================================
# # Merge bam files from the same sample but different lanes
#

rule merge_lanes_idxstats:
  input:
    DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped_sample+"{sample}/{sample}.idxstats.txt"
  shell:
    "samtools idxstats {input} > {output}"


rule merge_lanes_stat:
  input:
    DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped_sample+"{sample}/{sample}.stats.txt"
  shell:
    "samtools stats {input} > {output}"

rule merge_lanes_flagstat:
  input:
    DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
  output:
    DIR_mapped_sample+"{sample}/{sample}.flagstat.txt"
  shell:
    "samtools flagstat {input} > {output}"


rule merge_lanes_bismark:
     input:
         config['lanes_file'],
         infiles=lambda sample: [DIR_mapped+x+"/"+x+"_sorted_merged.bam" for x in SAMPLES_LANES[sample[0]] ]
     output:
         DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
     shell:
         "{tools}/sambamba merge -t 5 {output} {input.infiles}"



# # ==========================================================================================
# # Merge SE ad PE:
# 

rule idxstats_bismark_merged:
  input:
    DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
  output:
    DIR_mapped+"{sample}/{sample}_merged.idxstats.txt"
  shell:
    "samtools idxstats {input} > {output}"


rule stat_bismark_merged:
  input:
    DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
  output:
    DIR_mapped+"{sample}/{sample}_merged.stats.txt"
  shell:
    "samtools stats {input} > {output}"


rule flagstat_bismark_merged:
  input:
    DIR_mapped+"{sample}/{sample}_sorted_merged.bam"
  output:
    DIR_mapped+"{sample}/{sample}_merged.flagstat.txt"
  shell:
    "samtools flagstat {input} > {output}"

rule merge_bam_pe_and_se:
  input:
    DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam",
    DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam",
    DIR_mapped+"{sample}/{sample}_sorted.bam"
  output:
    temp(DIR_mapped+"{sample}/{sample}_sorted_merged.bam")
  params:
    merged=DIR_mapped+"{sample}/{sample}_merged.bam",
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}"
  benchmark:
    outputdir+"benchmarks/{sample}.merge_bam_pe_and_se.benchmark.txt"
  shell:
    """
    {tools}/sambamba merge {params.merged} {input}
    {tools}/sambamba sort {params.merged} --tmpdir={params.tmpdir} -o {output} {params.sort_args}
    """


# ==========================================================================================
# Map unaliagned reads:

## PBAT for R2
# ### Paired-end
# rule sort_index_bam_unmapped2_se_pbat:
#   input:
#     DIR_mapped+"{sample}/{sample}_unmapped_2_pbat.bam"
#   output:
#     DIR_mapped+"{sample}/{sample}_unmapped_2_pbat_sorted.bam"
#   params:
#     sort_args = config['args']['sambamba_sort'],
#     tmpdir=DIR_mapped+"{sample}/"
#   shell:
#     "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}"
# 
# rule align_unmapped2_se_pbat:
#     input:
#         DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#     output:
#         outfile=temp(DIR_mapped+"{sample}/{sample}_unmapped_2_pbat.bam"),
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
#     message: "Mapping unmapped reads as single-end to genome."
#     shell:
#         """
#         ln -s {input} {DIR_mapped}{wildcards.sample}/{wildcards.sample}_unmapped_2_pbat.fq.gz
#         {tools}/bismark {params.bismark_args} --pbat {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {DIR_mapped}{wildcards.sample}/{wildcards.sample}_unmapped_2_pbat.fq.gz
#         mv {output.outdir}{wildcards.sample}_unmapped_2_pbat_bismark_bt2.bam {output.outfile}
#         """

### Paired-end
rule sort_index_bam_unmapped2_se:
  input:
    DIR_mapped+"{sample}/{sample}_unmapped_2.bam"
  output:
    DIR_mapped+"{sample}/{sample}_unmapped_2_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort']#,
    #tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_sort2.log"
  benchmark:
    outputdir+"benchmarks/{sample}.sort_index_bam_unmapped2_se.benchmark.txt"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={DIR_mapped}/{wildcards.sample}/ -o {output} {params.sort_args}  > {log} 2> {log}.err"


rule align_unmapped2_se:
    input:
        DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
    output:
        outfile=temp(DIR_mapped+"{sample}/{sample}_unmapped_2.bam")
        #outdir=DIR_mapped+"{sample}/"
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
    message: "Mapping unmapped reads as single-end to genome."
    benchmark:
        outputdir+"benchmarks/{sample}.align_unmapped2_se.benchmark.txt"
    shell:
        """
        {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
        mv {DIR_mapped}/{wildcards.sample}/{wildcards.sample}_unmapped_2_bismark_bt2.bam {output.outfile}
        """

rule sort_index_bam_unmapped1_se:
  input:
    DIR_mapped+"{sample}/{sample}_unmapped_1.bam"
  output:
    DIR_mapped+"{sample}/{sample}_unmapped_1_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_sort1.log"
  benchmark:
    outputdir+"benchmarks/{sample}.sort_index_bam_unmapped1_se.benchmark.txt"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"


rule align_unmapped1_se:
    input:
        DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz"
    output:
        outfile=temp(DIR_mapped+"{sample}/{sample}_unmapped_1.bam")
        #outdir=DIR_mapped+"{sample}/"
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
    benchmark:
        outputdir+"benchmarks/{sample}.align_unmapped1_se.benchmark.txt"
    message: "Mapping unmapped reads as single-end to genome."
    shell:
        """
        {tools}/bismark {params.bismark_args} {params.genomeFolder} {params.outdir} {params.pathToBowtie} {params.samtools} {params.tmpdir} {input} > {log.align} 2> {log.align}.err
        mv {DIR_mapped}/{wildcards.sample}/{wildcards.sample}_unmapped_1_bismark_bt2.bam {output.outfile}
        """





