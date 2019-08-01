
# print(SAMPLES_TREATMENT.keys())
# #print(SAMPLES_LANES.keys())

# rule unite_meth_calls:
#      input:
#          [DIR_methcall+sample+"/tabix_CpG/"+sample+"_CpG_filtered.txt.bgz" for sample in [x for x in SAMPLES_TREATMENT.keys()]]
#      output:
#          destrandTfileT = DIR_methcall+"methylBase/methylBase_CpG_dT.RDS",
#          destrandFfileF = DIR_methcall+"methylBase/methylBase_CpG_dF.RDS",
#          destrandTfile_tbxT = DIR_methcall+"methylBase/methylBase_CpG_dT.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
#          destrandFfile_tbxF = DIR_methcall+"methylBase/methylBase_CpG_dF.txt.bgz"
#      params:
#          inputdir = DIR_methcall,
#          samples = " ".join(SAMPLES_LANES.keys()),
#          treatments = [SAMPLES_TREATMENT[sample] for sample in SAMPLES_LANES.keys()],
#          assembly=ASSEMBLY,
#          cores=24,
#          savedb=True,
#          dbdir = DIR_methcall+"methylBase/",
#          suffixT = "CpG_dT",
#          suffixF = "CpG_dF",
#      log: DIR_methcall+"methylBase/meth_unite_CpG.log"
#      shell:
#        """
#          {tools}/Rscript {DIR_scripts}/Unite_meth.R \
#                  --inputfiles="{input}" \
#                  --destrandTfile={output.destrandTfileT} \
#                  --destrandFfile={output.destrandFfileF} \
#                  --inputdir={params.inputdir} \
#                  --samples="{params.samples}" \
#                  --treatments="{params.treatments}" \
#                  --assembly="{params.assembly}" \
#                  --cores={params.cores} \
#                  --savedb={params.savedb} \
#                  --logFile={log} \
#                  --dbdir={params.dbdir} \
#                  --suffixT={params.suffixT} \
#                  --suffixF={params.suffixF}
#          """

rule filter_and_chromcanon_meth_calls:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix_CpG/{sample}_CpG_filtered.txt.bgz"
     params:
         mincov      = MINCOV_FILTER,
         save_folder = DIR_methcall+"{sample}"+"/tabix_CpG/",
         sample_id = "{sample}_CpG",
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
         subdir="/tabix_CpG/",
         mincov=config['args']['MINCOV']
     log:
         DIR_methcall+"{sample}/tabix_CpG/{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',mincov=as.numeric('{params.mincov}'), sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CpG',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
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
         subdir="/tabix_CHG/",
         mincov=config['args']['MINCOV']

     log:
         DIR_methcall+"{sample}/tabix_CHG/{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',mincov=as.numeric('{params.mincov}'),sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CHG',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
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
         subdir="/tabix_CHH/",
         mincov=config['args']['MINCOV']
     log:
         DIR_methcall+"{sample}/tabix_CHH/{sample}.txt.bgz.log"
     shell: # code below is an act of frustration of snakemake and R.
         """
         echo "require(methylKit); myobjDB=methRead('{input}',mincov=as.numeric('{params.mincov}'),sample.id='{params.sampleid}',assembly='{params.assembly}',treatment='{params.treatment}',context='CHH',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla > {log} 2> {log}.err
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
