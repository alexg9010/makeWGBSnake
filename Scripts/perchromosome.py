


# rule merge_united_methcalls:
#      input:
#          destrandTfile_perchrom = [DIR_methcall+"methylBase_per_chrom/"+chrom+"/methylBaseDB.obj_filtered_destrandF."+chrom+".RDS" for chrom in CHROMS_CANON]
#          #destrandFfile_perchrom = expand(DIR_methcall+"methylBase_per_chrom/methylBaseDB.obj_filtered_destrandF.{chrom}.RDS",chrom=CHROMS_CANON)
#      output:
#          destrandTfile_perchrom = DIR_methcall+"methylBase_per_chrom/methylBaseDB.obj_filtered_destrandF.RDS"
#          #destrandFfile_perchrom = DIR_methcall+"methylBase_per_chrom/methylBaseDB.obj_filtered_destrandF.RDS"
#      run:
#        R("""
# 
#        inputfiles = "{input.destrandTfile_perchrom}"
#        outputfile = "{output.destrandTfile_perchrom}"
#        print(inputfiles)
#        
#        library(methylKit)
#        
#        inputs <- strsplit(inputfiles, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
#        a = lapply(inputs, readRDS)
#        #print(a)
#        """)
# 
# rule unite_meth_calls_perchr:
#      input:
#          [DIR_methcall+sample+"/per_chrom/"+sample+"_{chrom}_cpg_filtered.txt.bgz" for sample in SAMPLES]
#      output:
#          destrandTfile_perchrom = DIR_methcall+"methylBase_per_chrom/{chrom}/methylBaseDB.obj_filtered_destrandT.{chrom}.RDS",
#          destrandFfile_perchrom = DIR_methcall+"methylBase_per_chrom/{chrom}/methylBaseDB.obj_filtered_destrandF.{chrom}.RDS",
#          destrandTfile_perchrom_tbx = DIR_methcall+"methylBase_per_chrom/{chrom}/methylBase_filtered_destrandT.{chrom}.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
#          destrandFfile_perchrom_tbx = DIR_methcall+"methylBase_per_chrom/{chrom}/methylBase_filtered_destrandF.{chrom}.txt.bgz"
#      params:
#          inputdir = DIR_methcall,
#          samples = SAMPLES,
#          treatments = TREATMENT,
#          assembly=ASSEMBLY,
#          cores=24,
#          savedb=True,
#          adbdir = DIR_methcall+"methylBase_per_chrom/{chrom}/",
#          suffixT = "filtered_destrandT.{chrom}",
#          suffixF = "filtered_destrandF.{chrom}",
#      log: DIR_methcall+"methylBase_per_chrom/meth_unite.log"
#      shell:
#        """
#          {tools}/Rscript {DIR_scripts}/Unite_meth.R \
#                  --inputfiles="{input}" \
#                  --destrandTfile={output.destrandTfile_perchrom} \
#                  --destrandFfile={output.destrandFfile_perchrom} \
#                  --inputdir={params.inputdir} \
#                  --samples="{params.samples}" \
#                  --treatments="{params.treatments}" \
#                  --assembly="{params.assembly}" \
#                  --cores={params.cores} \
#                  --savedb={params.savedb} \
#                  --logFile={log} \
#                  --adbdir={params.adbdir} \
#                  --suffixT={params.suffixT} \
#                  --suffixF={params.suffixF}
# 
#          """
# 
# 
# rule filter_and_canon_chroms_perchr:
#      input:
#          tabixfile     =  DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_cpg.txt.bgz"
#      output:
#          outputfile    = DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}_cpg_filtered.txt.bgz"
#      params:
#          mincov      = MINCOV,
#          save_folder = DIR_methcall+"{sample}/per_chrom/",
#          sample_id = "{sample}",
#          canon_chrs_file = chromcanonicalfile,
#          assembly    = ASSEMBLY,
#          hi_perc=99.9,
#          cores=10
#      log:
#          DIR_methcall+"{sample}/per_chrom/{sample}_{chrom}.meth_calls_filter.log"
#      message: ""
#      shell:
#        """
#          {tools}/Rscript {DIR_scripts}/Filter_meth.R \
#                  --tabixfile={input.tabixfile} \
#                  --mincov={params.mincov} \
#                  --hi_perc={params.hi_perc} \
#                  --save_folder={params.save_folder} \
#                  --sample_id={params.sample_id} \
#                  --assembly={params.assembly} \
#                  --cores={params.cores} \
#                  --canon_chrs_file={params.canon_chrs_file} \
#                  --logFile={log};
#          """


#############

# rule unite_meth_calls:
#      input:
#          [DIR_methcall+sample+"/"+sample+"_cpg_filtered.txt.bgz" for sample in SAMPLES]
#      output:
#          destrandTfile_perchrom = DIR_methcall+"methylBase/methylBaseDB.obj_filtered_destrandT.RDS",
#          destrandFfile_perchrom = DIR_methcall+"methylBase/methylBaseDB.obj_filtered_destrandF.RDS",
#          destrandTfile_perchrom_tbx = DIR_methcall+"methylBase/methylBase_filtered_destrandT.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
#          destrandFfile_perchrom_tbx = DIR_methcall+"methylBase/methylBase_filtered_destrandF.txt.bgz"
#      params:
#          inputdir = DIR_methcall,
#          samples = SAMPLES,
#          treatments = TREATMENT,
#          assembly=ASSEMBLY,
#          cores=24,
#          savedb=True,
#          adbdir = DIR_methcall+"methylBase/",
#          suffixT = "filtered_destrandT",
#          suffixF = "filtered_destrandF",
#      log: DIR_methcall+"methylBase/meth_unite.log"
#      shell:
#        """
#          {tools}/Rscript {DIR_scripts}/Unite_meth.R \
#                  --inputfiles="{input}" \
#                  --destrandTfile={output.destrandTfile_perchrom} \
#                  --destrandFfile={output.destrandFfile_perchrom} \
#                  --inputdir={params.inputdir} \
#                  --samples="{params.samples}" \
#                  --treatments="{params.treatments}" \
#                  --assembly="{params.assembly}" \
#                  --cores={params.cores} \
#                  --savedb={params.savedb} \
#                  --logFile={log} \
#                  --adbdir={params.adbdir} \
#                  --suffixT={params.suffixT} \
#                  --suffixF={params.suffixF}
# 
#          """


# rule filter_and_canon:
#      input:
#          tabixfile     =  DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
#      output:
#          outputfile    = DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz"
#      params:
#          mincov      = MINCOV,
#          save_folder = DIR_methcall+"{sample}/",
#          sample_id = "{sample}",
#          canon_chrs_file = chromcanonicalfile,
#          assembly    = ASSEMBLY,
#          hi_perc=99,
#          cores=10
#      log:
#          DIR_methcall+"{sample}/{sample}.meth_calls_filter.log"
#      message: ""
#      shell:
#        """
#          {tools}/Rscript {DIR_scripts}/Filter_meth.R \
#                  --tabixfile={input.tabixfile} \
#                  --mincov={params.mincov} \
#                  --hi_perc={params.hi_perc} \
#                  --save_folder={params.save_folder} \
#                  --sample_id={params.sample_id} \
#                  --assembly={params.assembly} \
#                  --cores={params.cores} \
#                  --canon_chrs_file={params.canon_chrs_file} \
#                  --logFile={log};
#          """

# rule merge_methCall_CpG:
#      input:
#          callFiles = [DIR_methcall+"{sample}/per_chrom/{sample}_"+chrom+"_cpg.txt.bgz" for chrom in CHROMS_CANON]
#      output:
#          callFiles_merged = DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
#      params:
#          cores=20,
#          outfilename = "{sample}_cpg.txt",
#          outdir = DIR_methcall+"{sample}/"
#      run:
#         R("""
#         inputfiles = "{input.callFiles}"
#         cores = as.numeric("{params.cores}")
#         outfilename = "{params.outfilename}"
#         outdir = "{params.outdir}"
#         
#         inputs <- strsplit(inputfiles, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
#         
#         
#         library(methylKit)
#         
#        # this works only on metylBase object I guess  
#        # methylKit:::mergeTabix(tabixList=inputs,
#        #                                  dir=outdir,
#        #                                  filename=outfilename,
#        #                                  mc.cores=cores,
#        #                                  all=FALSE)
#        
#        
#        methylKit:::mergeTabix(tabixList=c("/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg19/snakemake_test/06_methyl_calls/methylBase_per_chrom/chr1/methylBase_filtered_destrandF.chr1.txt.bgz","/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg19/snakemake_test/06_methyl_calls/methylBase_per_chrom/chr2/methylBase_filtered_destrandF.chr2.txt.bgz"
#        ),
#                                         dir="/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg19/snakemake_test/06_methyl_calls/",
#                                         filename="tmp",
#                                         mc.cores=20,
#                                         all=FALSE)
#         """)















############################## PER CHROMOSOME
rule
export_bigwig:
input:
seqlengths = os.path.join(DIR_mapped, "Refgen_" + ASSEMBLY + "_chromlengths.csv"),
rdsfile = [os.path.join(DIR_methcall, "{sample}/{sample}_sorted_dedup_" + chrom + "_methylRaw.RDS") for chrom in CHROMS]
output:
bw = [os.path.join(DIR_bigwig, "{sample}_" + chrom + ".bw") for chrom in CHROMS]
message: "Exporting bigwig files"
shell:
"""
 {tools}/Rscript {DIR_scripts}/export_bw.R \
                 {input.rdsfile} \
                 {input.seqlengths} \
                 ASSEMBLY \
                 {output} \
 """

# I dont use it currently
rule
filterSNPs:
input:
bam = expand(DIR_methcall + "{{sample}}/{{sample}}_sorted_dedup_{chrom}_methylRaw.RDS", chrom=CHROMS),
bed = "SNPFile.bed",
tabixfile = expand(DIR_methcall + "{{sample}}/Tabix/{{sample}}_{chrom}.txt.bgz", chrom=CHROMS)
output:
expand(DIR_methcall + "{{sample}}/{{sample}}.deduped_{chrom}_filteredSNP.RDS", chrom=CHROMS)
message: "Filtering SNPs"
shell:
"{tools}/Rscript {DIR_scripts}/FilterSNP.R {input.bam} {input.bed} {output}"

# I dont use it currently
rule
makeTabix:
input:
rdsfiles = [DIR_methcall + "{sample}/{sample}_sorted_dedup_" + chrom + "_methylRaw.RDS" for chrom in CHROMS],
txtfiles = [DIR_methcall + "{sample}/{sample}_sorted_dedup_" + chrom + "_CpG.txt" for chrom in CHROMS]
output:
tabixdir = DIR_methcall + "{sample}/Tabix/",
output = DIR_methcall + "{sample}/Tabix/{sample}_methyl.txt.bgz"
# conda:
# "Envs/env_min.yaml"
shell:
"""
{tools}/Rscript {DIR_scripts}/makeTabix.R {output} {input}
"""


def getTreatment(samplename):
    SAMPLE_TREAT_DICT[samplename]


# rule merge_methCall:
#     input:
#        expand(DIR_methcall+"{{sample}}/{{sample}}_sorted_dedup_{chrom}_methylRaw.RDS", chrom=CHROMS)
#     output:
#        DIR_methcall+"{sample}/{sample}_sorted_dedup_methylRaw.RDS"
#     params:
#        treatment = lambda wildcards: SAMPLE_TREAT_DICT[wildcards.sample],
#        mincov = MINCOV,
#     shell:
#        "{tools}/Rscript {DIR_scripts}/methCall_merge.R {output} {params.treatment} {params.mincov} {input}"


rule
bam_methCall_per_chr:
input:
bamfile = expand(DIR_bam_per_chrom + '{{sample}}/{{sample}}_sorted_dedup_{chrom}.bam', chrom=CHROMS)
output:
rdsfile = expand(DIR_methcall + "{{sample}}/{{sample}}_sorted_dedup_{chrom}_methylRaw.RDS", chrom=CHROMS),
callFile = expand(DIR_methcall + "{{sample}}/{{sample}}_sorted_dedup_{chrom}_CpG.txt", chrom=CHROMS)
params:
assembly = ASSEMBLY,
mincov = MINCOV,
minqual = MINQUAL,
context = "CpG",
# savedb      = False #TODO,
savefolder = DIR_methcall + "{sample}/"
log:
os.path.join(DIR_methcall, "{sample}/{sample}.deduped_meth_calls.log")
message: "Extract methylation calls from bam file."
shell:
"""
{tools}/Rscript {DIR_scripts}/methCall.R \
        --inBam={input.bamfile} \
        --assembly={params.assembly} \
        --mincov={params.mincov} \
        --minqual={params.minqual} \
        --rds={output.rdsfile} \
        --logFile={log}
"""

rule
bam_sort_index_per_chr:
input:
expand(DIR_bam_per_chrom + '{{sample}}/{{sample}}_dedup_{chrom}.bam', chrom=CHROMS)
output:
expand(DIR_bam_per_chrom + '{{sample}}/{{sample}}_sorted_dedup_{chrom}.bam', chrom=CHROMS)
params:
sort_args = config['args']['sambamba_sort'],
tmpdir = DIR_bam_per_chrom + "{sample}/"
shell:
"{tools}/sambamba sort {input} --tmpdir={params.tmpdir} -o {output} {params.sort_args}"

# https://github.com/daler/enhancer-snakemake-demo/blob/master/Snakefile
rule
split_bam_per_chr:
input:
DIR_deduped + "{sample}/{sample}_sorted_dedup.bam"
output:
expand(DIR_bam_per_chrom + '{{sample}}/{{sample}}_dedup_{chrom}.bam', chrom=CHROMS)
run:
sampleid = list(wildcards)[0]
for chrom in CHROMS:
    cmd = tools + "/sambamba slice " + input[
        0] + " " + chrom + " -o " + DIR_bam_per_chrom + sampleid + "/" + sampleid + "_dedup_" + chrom + ".bam"
    # print(cmd)
    shell(cmd)

##########################################################


#def input_bismark_align_unmapped_pe_as_se(wildcards):
#    if (wildcards.ext=="1"):
#       input=DIR_mapped+"{sample}/{sample}_unmapped_1.fq.gz"
#    elif (wildcards.ext=="2"):
#       input=DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#    return(input)

#rule bismark_align_unmapped_pe_as_se:
#    input:
#        DIR_mapped+"{sample}/{sample}_unmapped_2.fq.gz"
#    output:
#        out = DIR_mapped+"{sample}/{sample}_unmapped_{ext}.bam",
#        odir = DIR_mapped+"{sample}/"
#    params:
#        # intermediate files
#        inter_bam = DIR_mapped+"{sample}_unmapped_2_bismark_bt2.bam",
#        # bismark parameters
#        bismark_args = config['args']['bismark_unmapped'],
#        genomeFolder = "--genome_folder " + genomedir,
#        outdir = "--output_dir  "+DIR_mapped+"{sample}/",
#        pathToBowtie = "--path_to_bowtie " + config['tools'],
#        useBowtie2  = "--bowtie2 ",
#        samtools    = "--samtools_path "+ config['tools']+"samtools",
#        tempdir     = "--temp_dir "+DIR_mapped+"{sample}/",
#    log:
#        DIR_mapped+"{sample}/{sample}_bismark_pe_mapping_unmapped_reads_2.log"
#    message: "Mapping paired-end reads as single-end to genome."
#    run:
#        cmds=[
#        '{tools}/bismark {params} {input.fin2} > {log} 2> {log}.err'
#        'ln -s ' + output.odir + os.path.basename(input.fin2)+'_bismark_bt2_pe.bam'+ + ' {output.out2}'
#        ]
#        for c in cmds:
#           shell(c)
