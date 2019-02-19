# 
# # # ==========================================================================================
# # # Mapping:
# # 


rule flagstat_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bwameth_sorted.bam"
  output:
    DIR_mapped+"{sample}/{sample}.flagstat.txt"
  shell:
    "samtools flagstat {input} > {output}"


rule sort_index_bam_bwameth:
  input:
    DIR_mapped+"{sample}/{sample}.bwameth.bam"
  output:
    DIR_mapped+"{sample}/{sample}.bwameth_sorted.bam"
  params:
    sort_args = config['args']['sambamba_sort'],
    tmpdir=DIR_mapped+"{sample}/"
  log:
    DIR_mapped+"{sample}/{sample}_bwameth_sort.log"
  shell:
    "{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}  > {log} 2> {log}.err"


rule bwameth2bam_sort:
     input:
         DIR_mapped+"{sample}/{sample}.bwameth.sam"
     output:
         DIR_mapped+"{sample}/{sample}.bwameth.bam"
     shell:
         "samtools view -Sb {input} > {output}; rm {input}"


rule bwameth_align_pe:
     input:
         genomefile+".bwameth.c2t.sa",
         genomefile+".bwameth.c2t.amb",
         genomefile+".bwameth.c2t.ann",
         genomefile+".bwameth.c2t.pac",
         genomefile+".bwameth.c2t.bwt",
         genomefile+".bwameth.c2t",
         fin1 = DIR_trimmed_subset+"{sample}/{sample}_1_val_1.fq.gz",
         fin2 = DIR_trimmed_subset+"{sample}/{sample}_2_val_2.fq.gz",
     output:
         bam = DIR_mapped+"{sample}/{sample}.bwameth.sam"
     params:
        # Bismark parameters
         bwameth_args = config['args']['bwameth']
     log:
         DIR_mapped+"{sample}/{sample}_bwameth_pe_mapping.log"
     message: "Mapping paired-end reads to genome using bwa-meth."
     shell:
         "bwameth.py {params.bwameth_args} --reference {genomefile} {input.fin1} {input.fin2} > {output}"   
   
   
# reference=/fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Base/Genomes/hg38/hg38_canonical/hg38.sel.fa
# bwameth.py index $reference
# bwameth.py -t 12 --reference $reference /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7_1_val_1.fq.gz /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7_2_val_2.fq.gz   
#  
# bwameth.py -t 12 --reference /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Base/Genomes/hg38/hg38_canonical/hg38.sel.fa /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7_1_val_1.fq.gz /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane/subset_reads/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7_2_val_2.fq.gz > /fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane/04_mapping_bwameth/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7/QMQHSB-RUNID-0195-FLOWCELL-BHFMKYCCXY-LANE-7.bwameth.sam
      
# 
# # ==========================================================================================
# # Generate methyl-converted version of the reference genome:
#       

# rule bwameth_genome_preparation:
#     input:
#         ancient(genomefile)
#     output:
#         genomefile+".bwameth.c2t.sa",
#         genomefile+".bwameth.c2t.amb",
#         genomefile+".bwameth.c2t.ann",
#         genomefile+".bwameth.c2t.pac",
#         genomefile+".bwameth.c2t.bwt",
#         genomefile+".bwameth.c2t"
#     log:
#         genomedir+'bismark_genome_preparation_'+ASSEMBLY+'.log'
#     message: "Converting {ASSEMBLY} Genome into Bisulfite analogue with bwa-meth"
#     shell:
#         "bwameth.py index {input}"
#        
# bwameth.py index $REF
# bwameth.py --reference $REF some_R1.fastq.gz some_R2.fastq.gz > some.output.sam        
#         




        
        
        
        
        
        
