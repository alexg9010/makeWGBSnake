     
# ==========================================================================================
# Deduplication:
#
rule sort_index_dedup:
     input:
         DIR_deduped_picard+"{sample}/{sample}.dedup.bam"
     output:
         DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     params:
         sort_args = config['args']['sambamba_sort'],
         tmpdir=DIR_deduped_picard+"{sample}/"
     log:
         DIR_deduped_picard+"{sample}/{sample}_sort.log"
     shell:
         "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"

rule dedup_picard_bwameth:
     input:
         DIR_mapped_sample+"{sample}/{sample}_sorted.bam"
     output:
        outfile=DIR_deduped_picard+"{sample}/{sample}.dedup.bam",
        metrics = DIR_deduped_picard+"{sample}/{sample}.deduplication.metrics.txt"
     params:
         picard_MarkDuplicates_args = config['args']['picard_MarkDuplicates_args'],
         avail_mem = '10g'
     log:
         DIR_deduped_picard+"{sample}/{sample}.deduplication.log"
     message:
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          #"""{tools}/picard MarkDuplicates I={input} O={output.outfile} \ ######################TODO
          """{tools}/java -Xmx{params.avail_mem} -jar {tools}/../share/picard-2.19.0-0/picard.jar MarkDuplicates I={input} O={output.outfile} \
          M={output.metrics} \
          REMOVE_DUPLICATES=true AS=true {params.picard_MarkDuplicates_args} \
          > {log} \
          2> {log}.err"""


