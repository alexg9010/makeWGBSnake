
# ==========================================================================================
# Methylation calling:
# 



####################################################################### METHYLKIT [START] ################################ 
# 
# rule unite_meth_calls:
#      input:
#          [DIR_methcall+sample+"/"+sample+"_cpg_filtered.txt.bgz" for sample in SAMPLES]
#      output:
#          destrandTfileT = DIR_methcall+"methylBase/methylBase_cpg_dT.RDS",
#          destrandFfileF = DIR_methcall+"methylBase/methylBase_cpg_dF.RDS",
#          destrandTfile_tbxT = DIR_methcall+"methylBase/methylBase_cpg_dT.txt.bgz", # snakemake pretends that this file doesnt exist and removes it
#          destrandFfile_tbxF = DIR_methcall+"methylBase/methylBase_cpg_dF.txt.bgz"
#      params:
#          inputdir = DIR_methcall,
#          samples = SAMPLES,
#          treatments = TREATMENT,
#          assembly=ASSEMBLY,
#          cores=24,
#          savedb=True,
#          dbdir = DIR_methcall+"methylBase/",
#          suffixT = "cpg_dT",
#          suffixF = "cpg_dF",
#      log: DIR_methcall+"methylBase/meth_unite_cpg.log"
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
# 
# 
# rule filter_and_canon:
#      input:
#          tabixfile     =  DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
#      output:
#          outputfile    = DIR_methcall+"{sample}/{sample}_cpg_filtered.txt.bgz"
#      params:
#          mincov      = MINCOV,
#          save_folder = DIR_methcall+"{sample}",
#          sample_id = "{sample}_cpg",
#          canon_chrs_file = chromcanonicalfile,
#          assembly    = ASSEMBLY,
#          hi_perc=99.9,
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
#                  --logFile={log}
#          """
#          
# rule methCall_CpG:
#      input:
#          bamfile = DIR_deduped_picard+"{sample}/{sample}.dedup.sorted.bam"
#      output:
#          callFile = DIR_methcall+"{sample}/{sample}_cpg.txt.bgz"
#      params:
#          assembly    = ASSEMBLY,
#          mincov      = MINCOV,
#          minqual     = MINQUAL,
#          context     = "CpG",
#          save_db      = True,
#          save_folder = DIR_methcall+"{sample}/",
#          sample_id = "{sample}"
#      log:
#          DIR_methcall+"{sample}/{sample}.meth_calls.log"
#      message: "Extract methylation calls from bam file."
#      shell:
#        """
#           {tools}/Rscript {DIR_scripts}/methCall.R \
#                  --inBam={input.bamfile} \
#                  --assembly={params.assembly} \
#                  --mincov={params.mincov} \
#                  --minqual={params.minqual} \
#                  --context={params.context} \
#                  --save_db={params.save_db}  \
#                  --save_folder={params.save_folder}  \
#                  --sample_id={params.sample_id} \
#                  --logFile={log}
#        """
####################################################################### METHYLKIT [END] ################################ 


####################################################################### MethylDackel [START] ################################ 


#TMP NOTES:
#MethylDackel mbias reference_genome.fa alignments.sorted.bam output_prefix
#MethylDackel mbias $allcontexts $ignoreFlags $fasta $bam ${bam.baseName} --txt > ${bam.baseName}_methyldackel.txt 
# 
#genome=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Base/Genomes/hg38/nohaplo/hg38_nohaplo.fa
# input=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/05_deduplication/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4.dedup.sorted.bam
# output=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4_methyldacker.bedGraph
# MethylDackel extract $genome $input -o $output
# 
#input=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/05_deduplication/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4.dedup.sorted.bam
#basename=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4_methyldacker_mbias
#output=/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4_methyldacker.bedGraph
#MethylDackel mbias $genome $input $basename --txt > $basename'_methyldackel.txt'

# # MethylDackel
# inbam=/fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/05_deduplication/WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-7/WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-7.dedup.sorted.bam
# outfile=/fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/05_deduplication/WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-7/WYKWK3-RUNID-0188-FLOWCELL-BHF5H3CCXY-LANE-7.dedup.sorted.methyldacker.bedGraph
# genome=/fast/AG_Akalin/kwreczy/Projects/BIH_Neuroblastoma/Base/Genomes/hg38/hg38_canonical/hg38.sel.fa
# 
# MethylDackel extract $genome $inbam -o $outfile --methylKit --keepSingleton  --keepDiscordant --keepDiscordant -@ 20 --chunkSize 1000000
# #-d 10 -p 20
# #--minDepth 10
# MethylDackel mbias $genome $inbam mbias --keepSingleton  --keepDiscordant --keepDiscordant -@ 20


rule unite_meth_calls:
     input:
         [DIR_methcall+sample+"/tabix/"+sample+"_methyldacker_cpg_filtered.txt.bgz" for sample in SAMPLES] 
         #[DIR_methcall+sample+"/tabix/"+sample+".txt.bgz" for sample in SAMPLES]  #####################################
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
         
         

rule filter_and_chromcanon_meth_calls:
     input:
         tabixfile     =  DIR_methcall+"{sample}/tabix/{sample}.txt.bgz"
     output:
         outputfile    = DIR_methcall+"{sample}/tabix/{sample}_methyldacker_cpg_filtered.txt.bgz"
     params:
         mincov      = MINCOV,
         save_folder = DIR_methcall+"{sample}"+"/tabix/",
         sample_id = "{sample}_methyldacker_cpg",
         canon_chrs_file = chromcanonicalfile,
         assembly    = ASSEMBLY,
         hi_perc=100,#99.9, ############################################################################################
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


rule tabix_methylDackerfile:
     input:
         DIR_methcall+"{sample}/{sample}_methyldacker_CpG.methylKit"
     output:
         DIR_methcall+"{sample}/tabix/{sample}.txt.bgz",
         DIR_methcall+"{sample}/tabix/{sample}.txt.bgz.tbi"
     params:
         sampleid="{sample}",
         assembly=ASSEMBLY,
         subdir="/tabix/"
     shell: # code below is an act of frustration of snakemake and R.
       """
       echo "require(methylKit); myobjDB=methRead('{input}',sample.id='{params.sampleid}',assembly='{params.assembly}',treatment=0,context='CpG',dbtype = 'tabix',dbdir = paste0('{DIR_methcall}','{params.sampleid}','{params.subdir}'))" | {tools}/R --vanilla 
       """

       
########################## 
# myobjDB=methRead('{input}',
#            sample.id='{sample}',
#            assembly='{ASSEMBLY}',
#            treatment=0,
#            context="CpG",
#            dbtype = "tabix",
#            dbdir = DIR_methcall+"{sample}/tabix/"
#            )
           
#       # """
#       # echo "require(methylKit); methylKit:::makeMethTabix('{input}')" | {tools}/R --vanilla 
#       # """
       
# myobjDB=methRead('/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/236L96-RUNID-0163-FLOWCELL-AHFMF5CCXY-LANE-3/tabix/236L96-RUNID-0163-FLOWCELL-AHFMF5CCXY-LANE-3.txt.bgz',
#            sample.id='sample',
#            assembly='ASSEMBLY',
#            treatment=0,
#            context="CpG",
#            dbtype = "tabix"
#            #dbdir = DIR_methcall+"{sample}/tabix/"
#            )
       
# # # read the files to a methylRawListDB object: myobjDB 
# # # and save in databases in folder methylDB
# myobjDB=methRead(input,
#            sample.id="22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4",
#            assembly="hg18",
#            #treatment=1,
#            context="CpG",
#            dbtype = "tabix",
#            dbdir = "/fast/work/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg38/per_run_flowcell_lane_notrimming/06_methyl_calls/22X3H1-RUNID-0144-FLOWCELL-AHF3YNCCXY-LANE-4/"
#            )
# #        
# > methylRawDB.obj
# methylRawDB object with 154796 rows
# --------------
#   chr start    end strand coverage numCs numTs
# 1 chr1.266377  chr1 266377      R        1   100     0
# 2 chr1.280144  chr1 280144      R        1     0   100
# 3 chr1.280179  chr1 280179      R        1   100     0
# 4 chr1.280199  chr1 280199      R        1   100     0
# 5 chr1.870672  chr1 870672      R        1   100     0
# --------------
  


#################      
         
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

####################################################################### MethylDackel [END] ################################ 

       
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

rule dedup_picard_bwameth:
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




