
# ==========================================================================================
# Methylation calling:
# 

####################################################################### METHYLKIT [START] ################################ 
# 
rule unite_meth_calls_cpg:
     input:
         [DIR_methcall+sample+"/tabix_CpG/"+sample+"_cpg_filtered.txt.bgz" for sample in SAMPLES_LANES.keys()]
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


def filter_and_canon_in_func(w):
     if w.context=="CpG":
         return(DIR_methcall+w.sample+"/tabix_CpG/"+w.sample+"_cpg.txt.bgz")
     if w.context=="CHG":
         return(DIR_methcall+w.sample+"/tabix_CHG/"+w.sample+"_chg.txt.bgz")
     if w.context=="CHH":
         return(DIR_methcall+w.sample+"/tabix_CHH/"+w.sample+"_chh.txt.bgz")
         
rule filter_and_canon_chh:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix_CHH/{sample}_chh.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix_CHH/{sample}_CHH_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}/tabix_CHH/",
         sample_id = "{sample}_CHH",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=99.9,
         cores=10
     log:
         DIR_methcall+"{sample}/{sample}.CHH.meth_calls_filter.log"
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
         
rule filter_and_canon_chg:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix_CHG/{sample}_chg.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix_CHG/{sample}_CHG_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}/tabix_CHG/",
         sample_id = "{sample}_CHG",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=99.9,
         cores=10
     log:
         DIR_methcall+"{sample}/{sample}.CHG.meth_calls_filter.log"
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
rule filter_and_canon_cpg:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix_CpG/{sample}_cpg.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix_CpG/{sample}_CpG_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}/tabix_CpG/",
         sample_id = "{sample}_CpG",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=99.9,
         cores=10
     log:
         DIR_methcall+"{sample}/{sample}.CpG.meth_calls_filter.log"
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

rule methCall_CHH:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/tabix_CHH/{sample}_chh.txt.bgz"
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CHH",
         save_db      = True,
         save_folder = DIR_methcall+"{sample}/tabix_CHH/",
         sample_id = "{sample}"
     log:
         DIR_methcall+"{sample}/tabix_CHH/{sample}.CHH.meth_calls.log"
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

rule methCall_CHG:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/tabix_CHG/{sample}_chg.txt.bgz"
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CHG",
         save_db      = True,
         save_folder = DIR_methcall+"{sample}/tabix_CHG/",
         sample_id = "{sample}"
     log:
         DIR_methcall+"{sample}/tabix_CHG/{sample}.CHG.meth_calls.log"
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

rule methCall_CpG:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}_merged.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/tabix_CpG/{sample}_cpg.txt.bgz"
     params:
         assembly    = ASSEMBLY,
         mincov      = MINCOV,
         minqual     = MINQUAL,
         context     = "CpG",
         save_db      = True,
         save_folder = DIR_methcall+"{sample}/tabix_CpG/",
         sample_id = "{sample}"
     log:
         DIR_methcall+"{sample}/tabix_CpG/{sample}.CpG.meth_calls.log"
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
####################################################################### METHYLKIT [END] ################################ 








