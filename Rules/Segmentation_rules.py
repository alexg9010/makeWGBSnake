
# ==========================================================================================
# Segmentation:

rule meth_segments_perchr:
     input:
         inputfile     = DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_cpg_filtered.txt.bgz"
     output:
         grfile      = os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments_gr.RDS"),
         bedfile     = os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments.bed")
     params:
         methSegPng = DIR_seg+"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments.png",
         assembly = ASSEMBLY,
         sampleid = "{sample}"
     log:
         os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}.deduped_meth_segments.log")
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
                          {log}
          """





