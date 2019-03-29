
# ==========================================================================================
# Methylation calling:
# 


rule unite_meth_calls:
     input:
         [DIR_methcall+sample+"/"+sample+"_cpg_filtered.txt.bgz" for sample in SAMPLES]
     output:
         destrandTfileT = DIR_methcall+"methylBase/methylBase_cpg_dT.RDS",
         destrandFfileF = DIR_methcall+"methylBase/methylBase_cpg_dF.RDS",
         destrandTfile_tbxT = DIR_methcall+"methylBase/methylBase_cpg_dT.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
         destrandFfile_tbxF = DIR_methcall+"methylBase/methylBase_cpg_dF.txt.bgz"
     params:
         inputdir = DIR_methcall,
         samples = SAMPLES,
         treatments = TREATMENT,
         assembly=ASSEMBLY,
         cores=24,
         savedb=True,
         dbdir = DIR_methcall+"methylBase/",
         suffixT = "cpg_dT",
         suffixF = "cpg_dF",
     log: DIR_methcall+"methylBase/meth_unite_cpg.log"
     shell:
       """
         {tools}/Rscript {DIR_scripts}/Unite_meth.R \
                 --inputfiles="{input}" \
                 --destrandTfile={output.destrandTfileT} \
                 --destrandFfile={output.destrandFfileF} \
                 --inputdir={params.inputdir} \
                 --samples="{params.samples}" \
                 --treatments="{params.treatments}" \
                 --assembly="{params.assembly}" \
                 --cores={params.cores} \
                 --savedb={params.savedb} \
                 --logFile={log} \
                 --dbdir={params.dbdir} \
                 --suffixT={params.suffixT} \
                 --suffixF={params.suffixF}
         """


rule filter_and_canon:
     input:
         tabixfile     =  DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}",
         sample_id = "{sample}_cpg",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=99.9,
         cores=10
     log:
         DIR_methcall+"{sample}/{sample}.meth_calls_filter.log"
     message: ""
     shell:
       """
         {tools}/Rscript {DIR_scripts}/Filter_meth.R \
                 --tabixfile={input.tabixfile} \
                 --mincov={params.mincov} \
                 --hi_perc={params.hi_perc} \
                 --save_folder={params.save_folder} \
                 --sample_id={params.sample_id} \
                 --assembly={params.assembly} \
                 --cores={params.cores} \
                 --canon_chrs_file={params.canon_chrs_file} \
                 --logFile={log}
         """



rule methCall_CpG:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CpG",
         save_db      = True,
         save_folder = DIR_methcall+"{sample}/",
         sample_id = "{sample}"
     log:
         DIR_methcall+"{sample}/{sample}.meth_calls.log"
     message: "Extract methylation calls from bam file."
     shell:
       """
          {tools}/Rscript {DIR_scripts}/methCall.R \
                 --inBam={input.bamfile} \
                 --assembly={params.assembly} \
                 --mincov={params.mincov} \
                 --minqual={params.minqual} \
                 --context={params.context} \
                 --save_db={params.save_db}  \
                 --save_folder={params.save_folder}  \
                 --sample_id={params.sample_id} \
                 --logFile={log}
       """

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

rule dedup_picard:
     input:
         DIR_mapped+"{sample}/{sample}.bwameth_sorted.bam" #####################################TODO ###### BWA-meth
     output:
        outfile=DIR_deduped_picard+"{sample}/{sample}.dedup.bam",
        metrics = DIR_deduped_picard+"{sample}/{sample}.deduplication.metrics.txt"
     params:
         picard_MarkDuplicates_args = config['args']['picard_MarkDuplicates_args']
     log:
         DIR_deduped_picard+"{sample}/{sample}_deduplication.log"
     message:
          "Deduplicating paired-end aligned reads from {input}"
     shell:
          """{tools}/picard MarkDuplicates I={input} O={output.outfile} \
          M={output.metrics} \
          REMOVE_DUPLICATES=true AS=true {params.picard_MarkDuplicates_args} \
          > {log} \
          2> {log}.err"""





