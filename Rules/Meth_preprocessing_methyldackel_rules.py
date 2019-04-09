
rule unite_meth_calls:
     input:
         [DIR_methcall+sample+"/tabix_CpG/"+sample+"_methyldacker_cpg_filtered.txt.bgz" for sample in SAMPLES_LANES.keys()]
     output:
         destrandTfileT = DIR_methcall+"methylBase/methylBase_cpg_dT.RDS",
         destrandFfileF = DIR_methcall+"methylBase/methylBase_cpg_dF.RDS",
         destrandTfile_tbxT = DIR_methcall+"methylBase/methylBase_cpg_dT.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
         destrandFfile_tbxF = DIR_methcall+"methylBase/methylBase_cpg_dF.txt.bgz"
     params:
         inputdir = DIR_methcall,
         samples = " ".join(SAMPLES_LANES.keys()),
         treatments = [SAMPLES_TREATMENT[sample] for sample in SAMPLES_LANES.keys()],
         assembly=ASSEMBLY,
         cores=24,
         savedb=True,
         dbdir = DIR_methcall+"methylBase/",
         suffixT = "_cpg_dT",
         suffixF = "_cpg_dF",
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
# "H9H1KG"
# > methylRawListDB.obj_filtered[[18]]
# methylRawDB object with 1 rows
# --------------
#   chr    start      end strand coverage numCs numTs
# 1 chr3 93470645 93470645      +       12     1    11

# in1="/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/H9H1KG/tabix/H9H1KG_methyldacker_cpg_filtered.txt.bgz"
# in1="/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/H9H1KG/tabix/H9H1KG.txt.bgz"
# 
# tmp = methRead(in1, 
#          "a" , 
#          "a", 
#          dbtype='tabix')
# 
# filterByCoverage(tmp,
#                  lo.count=10,
#                  lo.perc=NULL,
#                  hi.count=NULL,
#                  hi.perc=NULL,
#                  #dbdir=save_folder,
#                  dbtype="tabix"#,
#                  #save.db = TRUE
#                  )



         

rule filter_and_chromcanon_meth_calls:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix_CpG/{sample}_methyldacker_cpg_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}"+"/tabix_CpG/",
         sample_id = "{sample}_methyldacker_cpg",
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

rule tabix_methylDackerfile_CpG:
     input:
         DIR_methcall+"{sample}/{sample}_methyldacker_CpG.methylKit"
     output:
         DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz",
         DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz.tbi"
     params:
         sampleid="{sample}",
         assembly=ASSEMBLY,
         treatment=lambda sampleid: SAMPLES_TREATMENT[sampleid[0]],
         subdir="/tabix_CpG/"
     log:
         DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CpG',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
         """
         
rule tabix_methylDackerfile_CHG:
     input:
         DIR_methcall+"{sample}/{sample}_methyldacker_CHG.methylKit"
     output:
         DIR_methcall+"{sample}/tabix_CHG/{sample}.txt.bgz",
         DIR_methcall+"{sample}/tabix_CHG/{sample}.txt.bgz.tbi"
     params:
         sampleid="{sample}",
         assembly=ASSEMBLY,
         treatment=lambda sampleid: SAMPLES_TREATMENT[sampleid[0]],
         subdir="/tabix_CHG/"
     log:
         DIR_methcall+"{sample}/tabix_CHG/{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CHG',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
         """
         
rule tabix_methylDackerfile_CHH:
     input:
         DIR_methcall+"{sample}/{sample}_methyldacker_CHH.methylKit"
     output:
         DIR_methcall+"{sample}/tabix_CHH/{sample}.txt.bgz",
         DIR_methcall+"{sample}/tabix_CHH/{sample}.txt.bgz.tbi"
     params:
         sampleid="{sample}",
         assembly=ASSEMBLY,
         treatment=lambda sampleid: SAMPLES_TREATMENT[sampleid[0]],
         subdir="/tabix_CHH/"
     log:
         DIR_methcall+"{sample}/tabix_CHH//{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CHH',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
         """


rule methylDacker_mbias:
     input:
         DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     output:
         txt=DIR_methcall+"{sample}/{sample}_mbias_methyldackel.txt",
         p1=DIR_methcall+"{sample}/{sample}_mbias_OB.svg",
         p2=DIR_methcall+"{sample}/{sample}_mbias_OT.svg"
     params:
         methylDacker_args = config['args']['methylDacker_mbias'],
         basename=DIR_methcall+"{sample}/{sample}_mbias"
     log:
         DIR_methcall+"{sample}/{sample}.methyldacker_mbias.log"
     message: "Calculate methylation bias using MethylDackel."
     shell:
       """
       {tools}/MethylDackel mbias {genomefile} {input} {params.basename} --txt {params.methylDacker_args} > {output.txt} 2> {log}.err
       """  


rule methylDacker_CHH:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/{sample}_methyldacker_CHH.methylKit"
     params:
         methylDacker_args = config['args']['methylDacker_methcalling'],
         out=DIR_methcall+"{sample}/{sample}_methyldacker"
     log:
         DIR_methcall+"{sample}/{sample}.methyldacker_calls.log"
     message: "Extract methylation calls in CHH context from bam file using MethylDackel."
     shell:
       """
       {tools}/MethylDackel extract --CHH {genomefile} {input.bamfile} -o {params.out} {params.methylDacker_args}
       """ 


rule methylDacker_CHG:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/{sample}_methyldacker_CHG.methylKit"
     params:
         methylDacker_args = config['args']['methylDacker_methcalling'],
         out=DIR_methcall+"{sample}/{sample}_methyldacker"
     log:
         DIR_methcall+"{sample}/{sample}.methyldacker_calls.log"
     message: "Extract methylation calls in CHG context from bam file using MethylDackel."
     shell:
       """
       {tools}/MethylDackel extract --CHG {genomefile} {input.bamfile} -o {params.out} {params.methylDacker_args}
       """ 

rule methylDacker_CpG:
     input:
         bamfile = DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
     output:
         callFile = DIR_methcall+"{sample}/{sample}_methyldacker_CpG.methylKit"
     params:
         methylDacker_args = config['args']['methylDacker_methcalling'],
         out=DIR_methcall+"{sample}/{sample}_methyldacker"
     log:
         DIR_methcall+"{sample}/{sample}.methyldacker_calls.log"
     message: "Extract methylation calls from bam file using MethylDackel."
     shell:
       """
       {tools}/MethylDackel extract {genomefile} {input.bamfile} -o {params.out} {params.methylDacker_args}
       """ 
