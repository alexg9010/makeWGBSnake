
# ==========================================================================================
# Segmentation:

rule meth_segments:
     input:
         inputfile     = DIR_methcall+"{sample}/tabix_CpG/{sample}_CpG_filtered.txt.bgz"
     output:
         grfile      = os.path.join(DIR_seg,"{sample}/{sample}_segments_gr.RDS"),
         bedfile     = os.path.join(DIR_seg,"{sample}/{sample}.segments.bed")
     params:
         methSegPng = DIR_seg+"{sample}/{sample}.segments.png",
         assembly = ASSEMBLY,
         sampleid = "{sample}"
     log:
         os.path.join(DIR_seg,"{sample}/{sample}_segments.log")
     message: "Segmenting methylation profile for {input.inputfile}."
     shell:
         """
          {tools}/Rscript {DIR_scripts}/methSeg.R \
                          {input.inputfile} \
                          {output.grfile} \
                          {output.bedfile} \
                          {params.methSegPng} \
                          {params.assembly} \
                          {params.sampleid} \
                          {log} 2> {log}.err
          """





