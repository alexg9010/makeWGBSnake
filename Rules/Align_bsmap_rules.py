# 
# # # ==========================================================================================
# # # Mapping:
# # 


# # ==========================================================================================
# # Merge bam files from the same sample but different lanes
#  

# rule merge_lanes_idxstats:
#   input:
#     DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
#   output:
#     DIR_mapped_sample+"{sample}/{sample}.idxstats.txt"
#   shell:
#     "samtools idxstats {input} > {output}"


# rule merge_lanes_stat:
#   input:
#     DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
#   output:
#     DIR_mapped_sample+"{sample}/{sample}.stats.txt"
#   shell:
#     "samtools stats {input} > {output}"

# rule merge_lanes_flagstat:
#   input:
#     DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
#   output:
#     DIR_mapped_sample+"{sample}/{sample}.flagstat.txt"
#   shell:
#     "samtools flagstat {input} > {output}"
    
    
# rule merge_lanes_bwameth:
#      input:
#          config['lanes_file'],
#          infiles=lambda sample: [DIR_mapped+x+"/"+x+".bwameth_sorted.bam" for x in SAMPLES_LANES[sample[0]] ]
#      output:
#          DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
#      shell:
#          "{tools}/sambamba merge -t 5 {output} {input.infiles}"


# # ==========================================================================================
# # Align. stats
#  

rule idxstats_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bsmap_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.bsmap.idxstats.txt"
  shell:
    "samtools idxstats {input} > {output}"


rule stat_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bsmap_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.bsmap.stats.txt"
  shell:
    "samtools stats {input} > {output}"


rule flagstat_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bsmap_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.bsmap.flagstat.txt"
  shell:
    "samtools flagstat {input} > {output}"


rule sort_index_bam_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bsmap.bam"
  output:
    DIR_mapped+"{sample}/{sample}.bsmap_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_bwameth_sort.log"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"


# # ==========================================================================================
# # Align whole data sets
#    

rule bsmap_align_pe:
    input:
        genomefile=genomefile,
        fin1 = DIR_trimmed+"{sample}/{sample}_1_val_1.fq.gz",
        fin2 = DIR_trimmed+"{sample}/{sample}_2_val_2.fq.gz"
    output:
        bam = DIR_mapped+"{sample}/{sample}.bsmap.bam"
    #params:
        # bwa-meth parameters
        #bwameth_args = config['args']['bsmap']
    log:
        DIR_mapped+"{sample}/{sample}_bsmap_pe_mapping.log"
    message: "Mapping paired-end reads to genome using BSMAP."
    shell:
      """
       set -o pipefail
       {tools}/bsmap \\
       -d  {input.genomefile} \\
       -a {input.fin1} -b {input.fin2} \\
       -o {output.bam} \\
       -2 out_upair.bsp -p 8 -v 5 -l 8
      """
        
   



        
        
        
        
        
        
