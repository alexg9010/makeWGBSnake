# ==========================================================================================
# Differential methylation in a pairwise manner



rule diffmeth_pairwise:
     input:
         destrandTfile = DIR_methcall+"methylBase_per_chrom/methylBase_cpg_dT.txt.bgz"
     output:
         outfile=DIR_diffmeth+'diffmeth_{treat}.RDS'
     params:
         treat="{treat}",
         cores = 4, # with more than 5 usually there is not enough memory
         treatments = TREATMENT,
         sampleids = SAMPLES,
         context = "CpG",
         assembly=ASSEMBLY,
         outputdir = DIR_diffmeth,
         suffix = "diffmeth_{treat}",
         save_db = True
     shell:
         """{tools}/Rscript {DIR_scripts}/MethDiff.R \
         {input} {output} \
         {params.treat} {params.cores} "{params.treatments}" "{params.sampleids}" \
         {params.context} {params.assembly} \
         {params.outputdir} {params.suffix} {params.save_db}"""
        
        
rule diffmeth_pairwise_merged:
     input:
         destrandTfile = DIR_methcall+"methylBase_per_chrom/methylBase_merged_cpg_dT.txt.bgz"
     output:
         outfile=DIR_diffmeth+'diffmeth_merged_{treat}.RDS'
     params:
         treat="{treat}",
         cores = 4, # with more than 5 usually there is not enough memory
         treatments = TREATMENT,
         sampleids = SAMPLES,
         context = "CpG",
         assembly=ASSEMBLY,
         outputdir = DIR_diffmeth,
         suffix = "diffmeth_{treat}",
         save_db = True
     shell:
         """{tools}/Rscript {DIR_scripts}/MethDiff.R \
         {input} {output} \
         {params.treat} {params.cores} "{params.treatments}" "{params.sampleids}" \
         {params.context} {params.assembly} \
         {params.outputdir} {params.suffix} {params.save_db}"""
        
