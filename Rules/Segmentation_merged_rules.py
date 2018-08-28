
# ==========================================================================================
# Segmentation:



# rule merge_segments_pe_se:
#      input:
#        [DIR_seg + "{sample}/per_chrom/{sample}_"+chrom+"_merged.deduped_meth_segments_gr.RDS" for chrom in CHROMS_CANON]
#      output:
#        grfile      = os.path.join(DIR_seg,"{sample}/{sample}.pe.se.deduped_meth_segments_gr.RDS"),
#        bedfile     = os.path.join(DIR_seg,"{sample}/{sample}.pe.se.deduped_meth_segments.bed")
#      params:
#        assembly = ASSEMBLY,
#        sampleid = "{sample}",
#        cores = 20
#      run:
#        R("""
#        library(GenomicRanges)
#        library(parallel)
# 
#        inputfiles = "{input}"
#        outgrfile="{output.grfile}"
#        outbedfile = "{output.bedfile}"
#        sample = "{params.sampleid}"
#        assembly = "{params.assembly}"
#        cores = as.numeric("{params.cores}")
# 
#        inputs <- strsplit(inputfiles, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
# 
#        list.gr = mclapply(inputs, readRDS, mc.cores=cores)
#        print(list.gr)
# 
# 
#        """)

rule segment_pe_se:
     input:
         inputfile     = DIR_methcall+"{sample}/{sample}_merged_cpg_filtered.txt.bgz"
     output:
         grfile      = os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments_gr.RDS"),
         bedfile     = os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments.bed")
     params:
         methSegPng = DIR_seg+"{sample}/{sample}_merged.deduped_meth_segments.png",
         assembly = ASSEMBLY,
         sampleid = "{sample}",
         joinneighbours = True,
         initializeonsubset=0.4,
         maxInt=100
     log:
         os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments.log")
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
                          {params.joinneighbours} \
                          {params.initializeonsubset} \
                          {params.maxInt} \
                          {log}
          """


rule segment_perchr_pe_se:
     input:
         inputfile     = DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_merged_cpg_filtered.txt.bgz"
     output:
         grfile      = os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}_merged.deduped_meth_segments_gr.RDS"),
         bedfile     = os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}_merged.deduped_meth_segments.bed")
     params:
         methSegPng = DIR_seg+"{sample}/per_chrom/{sample}_{chrom}_merged.deduped_meth_segments.png",
         assembly = ASSEMBLY,
         sampleid = "{sample}"
     log:
         os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_{chrom}_merged.deduped_meth_segments.log")
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


# rule segments_perchr_pe_se:
#      input:
#          inputfile     = DIR_methcall+"{sample}/{sample}_merged_cpg_filtered.txt.bgz"
#      output:
#          grfile      = os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments_gr.RDS"),
#          bedfile     = os.path.join(DIR_seg,"{sample}/{sample}_merged.deduped_meth_segments.bed")
#      params:
#          methSegPng = DIR_seg+"{sample}/{sample}_merged.deduped_meth_segments.png",
#          assembly = ASSEMBLY,
#          sampleid = "{sample}",
#          cores=24
#      log:
#          os.path.join(DIR_seg,"{sample}/per_chrom/{sample}_merged.deduped_meth_segments.log")
#      message: "Segmenting methylation profile for {input.inputfile}."
#      shell:
#          """
#           {tools}/Rscript {DIR_scripts}/methSeg.R \
#                           {input.inputfile} \
#                           {output.grfile} \
#                           {output.bedfile} \
#                           {params.methSegPng} \
#                           {params.assembly} \
#                           {params.sampleid} \
#                           {params.cores} \ ############################################# ADDDDDDDDDDDDDDDDD
#                           {log}
#           """




