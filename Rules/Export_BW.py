
# ==========================================================================================
# Export a bigwig file:


rule export_bigwig:
    input:
        seqlengths = chromcanonicalfile,
        tabixfile    = DIR_methcall+"{sample}/tabix_CpG/{sample}_CpG_filtered.txt.bgz"
    params:
        sampleid = "{sample}"
    output:
        bw         = DIR_bigwig+"{sample}/{sample}.bw",
    message: "Exporting bigwig files."
    shell:
       """
         {tools}/Rscript {DIR_scripts}/export_bw.R \
                 {input.tabixfile} \
                 {input.seqlengths} \
                 {ASSEMBLY} \
                 {params.sampleid} \
                 {output}
       """
       
