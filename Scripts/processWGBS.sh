#!/bin/bash

# A desperate intitiative to run things fast.
# This script will be rewritten to snakemake.

id=$1
outputdir=$2
librarytype=$3
maxins=$4
score=$5
clip_R1=$6
three_prime_clip_R1=$7
clip_R2=$8
three_prime_clip_R2=$9
memorylimit="200G"

echo $outputdir
echo $score

# For all reads
inputdir=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_merged_lanes/
#outputdir=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/hg19/bsseq_pigx/

# For the random 10M read per lane
#inputdir=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Data/Raw_merged_random10Mlanespe/
#outputdir="/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/random10Mperlane_hg19_se/"
#outputdir="/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/random10Mperlane_hg19_se_gap2k_L-0.6-0.6/"

mkdir -p $outputdir 

assembly="hg19"
genomeprepoutput=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/hg19/bsseq_pigx/genome_prep.log 
genomedir=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Base/Genomes/hg19/

chromsizes="/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Base/Annotation/hg19/hg19.chrom.sizes"

ptspath=/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Scripts/Rscripts/
myprogs="/fast/users/kwreczy_m/programs/bin/miniconda3/envs/mybase/bin/"
rscriptspath="/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Scripts/Rscripts/"

#-----------------------------------------------------------------

fastqc=false
trimgalore=false
fastqcaftertrimgalore=false
bismark=true
dedup=true
splitperchr=true
methcall=true
bigwig=false
multiqc=false

#-----------------------------------------------------------------

#-------------------------------
# Quality of reads

outputqc=$outputdir/01_raw_QC/$id/
mkdir -p $outputqc
if [ "$fastq" == true ] ; then
 # Fastqc
 echo "FASTQC" 
 
 if [ "$librarytype" == "paired" ] ; then

  $myprogs/fastqc $inputdir/$id'_1'.fastq.gz --outdir $outputqc 2> $outputqc/$id'_1'.log
  $myprogs/fastqc $inputdir/$id'_2'.fastq.gz --outdir $outputqc 2> $outputqc/$id'_2'.log

 elif [ "$librarytype" == "single" ] ; then
 
  $myprogs/fastqc $inputdir/$id.fastq.gz --outdir $outputqc 2> $outputqc/$id.log

 fi
fi

#-------------------------------
# Trim reads

trimgaloredir=$outputdir/02_trimming/$id/

if [ "$trimgalore" == true ] ; then
 # Trimgalore
 echo "TRIMGALORE"

 mkdir -p $trimgaloredir
 
 echo $librarytype

 if [ "$librarytype" == "paired" ] ; then
 
  $myprogs/trim_galore --clip_R1 $clip_R1 --three_prime_clip_R1 $three_prime_clip_R1 --clip_R2 $clip_R2 --three_prime_clip_R2 $three_prime_clip_R2 --paired $inputdir/$id'_1'.fastq.gz $inputdir/$id'_2'.fastq.gz  --output_dir $trimgaloredir/ --phred33 --gzip --path_to_cutadapt $myprogs/cutadapt 2> $trimgaloredir/$id.log

 elif [ "$librarytype" == "single" ] ; then
   echo "SINGLE-END"
   $myprogs/trim_galore --clip_R1 $clip_R1 --three_prime_clip_R1 $three_prime_clip_R1 $inputdir/$id'.fastq.gz'  --output_dir $trimgaloredir/ --phred33 --gzip --path_to_cutadapt $myprogs/cutadapt 2> $trimgaloredir/$id.log

 fi
fi

#-------------------------------
fastqafter=$outputdir/03_posttrimming_QC/$id/
mkdir -p $fastqafter

if [ "$fastqcaftertrimgalore" == true ]; then
 # Fastqc after trim galore
 echo "FASTQC after TRIMGALORE"
 
  if [ "$librarytype" == "paired" ]; then
  
   $myprogs/fastqc $trimgaloredir/$id'_1_val_1.fq.gz'  --outdir $fastqafter 2> $fastqafter/$id'_1'.log
   $myprogs/fastqc $trimgaloredir/$id'_2_val_2.fq.gz' --outdir $fastqafter 2> $fastqafter/$id'_2'.log

  elif [ "$librarytype" == "single" ]; then
  
   $myprogs/fastqc $trimgaloredir/$id'_trimmed.fq.gz' --outdir $fastqafter 2> $fastqafter/$id.log

 fi
fi


#-------------------------------
## Prepare genome index for bismark
# myprogs/bismark_genome_preparation --path_to_bowtie $myprogs  --bowtie2  --verbose $genomedir  > $genomeprepoutput 2> $genomeprepoutput.err

#-------------------------------
# Mapping step

mapoutputdir=$outputdir/04_mapping/
if [ "$bismark" == true ]; then
 ## Map into a given genome
 echo "BISMARK"
 mkdir -p $mapoutputdir
 maplog=$mapoutputdir/$id.log
 maperr=$mapoutputdir/$id.err

 if [ "$librarytype" == "paired" ] ; then
  file1=$trimgaloredir/$id'_1_val_1.fq.gz'
  file2=$trimgaloredir/$id'_2_val_2.fq.gz'
  
  $myprogs/bismark $directional --maxins $maxins --score_min $score --ambig_bam  --unmapped --ambiguous  --genome_folder $genomedir --output_dir $mapoutputdir/$id/  --bowtie2 --path_to_bowtie $myprogs  --samtools_path  $myprogs/samtools --temp_dir $mapoutputdir/$id/ --multicore 6 -N 0 -L 15 -1 $file1 -2 $file2 > $maplog 2> $maperr  

  # Sort bam files after mapping
  $myprogs/sambamba sort $mapoutputdir/$id/$id'_1_val_1_bismark_bt2_pe'.bam --tmpdir=$mapoutputdir/$id/ -o $mapoutputdir/$id/$id'_sorted'.bam -t 1 --memory-limit=$memorylimit
  # this is too slow #myprogs/samtools sort -@ 1 $mapoutputdir/$id/$id'_1_val_1_bismark_bt2_pe'.bam > $mapoutputdir/$id/$id'_1_val_1_bismark_bt2_pe_sorted'.bam

 elif [ "$librarytype" == "single" ] ; then
  file1=$trimgaloredir/$id'_trimmed.fq.gz'
  
  $myprogs/bismark --non_directional --score_min $score --ambig_bam  --unmapped --ambiguous  --genome_folder $genomedir --output_dir $mapoutputdir/$id/  --bowtie2 --path_to_bowtie $myprogs  --samtools_path  $myprogs/samtools --temp_dir $mapoutputdir/$id/ --multicore 6 -N 0 -L 15 $file1 > $maplog 2> $maperr  

  # Sort bam files after mapping
  $myprogs/sambamba sort $mapoutputdir/$id/$id'_trimmed_bismark_bt2'.bam --tmpdir=$mapoutputdir/$id/ -o $mapoutputdir/$id/$id'_sorted'.bam -t 1 --memory-limit=$memorylimit

 fi

fi


#------------------------ 
# Map unligned reads as single-end


#-----------------------
# Include mapped unliagned reads as single end into a original bam file


#-------------------------
multiqcdir=$outputdir/MultiQC/$id/
if [ "$multiqc" == true ]; then
  echo "multiqc"
  $myprogs/multiqc $trimgaloredir $fastqafter $mapoutputdir -o $multiqcdir 
fi

#-------------------------------

dedupdir=$outputdir/05_picard_MarkDuplicates/$id/
if [ "$dedup" == true ]; then
 # Samtools deduplication
 echo "DEDUPLICATION"

 #dedupdir=$outputdir/05_samtools_dedup/$id/
 #mkdir -p $dedupdir
 # http://www.htslib.org/doc/samtools.html
 # The first sort can be omitted if the file is already name ordered
 #myprogs/sambamba sort --tmpdir $dedupdir/  --memory-limit=30G -n -o $dedupdir/$id'_namesort'.bam $mapoutputdir/$id/$id'_1_val_1_bismark_bt2_pe'.bam
 # Add ms and MC tags for markdup to use later
 #myprogs/samtools fixmate -m $dedupdir/$id'_namesort'.bam $dedupdir/$id'_fixmate'.bam
 # Markdup needs position order
 #myprogs/sambamba sort -o $dedupdir/$id'_positionsort'.bam $dedupdir/$id'_fixmate'.bam --memory-limit=30G --tmpdir $dedupdir/
 # Finally mark duplicates
#myprogs/samtools markdup -r -T $id -s $dedupdir/$id'_positionsort'.bam $dedupdir/$id'_markdup'.bam > $dedupdir/$id'_markdup'.log 2> $dedupdir/$id'_markdup'.err

 # Picard's deduplication
 mkdir -p $dedupdir
 $myprogs/picard MarkDuplicates I=$mapoutputdir/$id/$id'_sorted'.bam O=$dedupdir/$id'_dedup'.bam M=$dedupdir/$id'_dup_metrics'.txt REMOVE_DUPLICATES=true AS=true > $dedupdir/$id'_dedup'.log 2> $dedupdir/$id'_dedup'.err
 $myprogs/sambamba sort $dedupdir/$id'_dedup'.bam > $dedupdir/$id'_dedup.sorted'.bam --memory-limit=30G --tmpdir=$dedupdir/ -t 1
 $myprogs/samtools index $dedupdir/$id'_dedup.sorted'.bam

fi

if [ "$splitperchr" == true ]; then
echo "Split a bam file by chromosome"

# this is too slow 
# $myprogs/bamtools split -in $dedupdir/$id'_dedup.sorted'.bam -reference

# sambamba is faster than samtools
for chromosome in $(cat /fast/projects/peifer_wgs/work/2017-12-19_WGBS/Base/Genomes/hg19/chroms.txt | cut -c 2-); do 
   $myprogs/sambamba slice $dedupdir/$id'_dedup.sorted'.bam $chromosome -o $dedupdir/$chromosome'_dedup.sorted'.bam

done  

fi


#-------------------------------
if [ "$methcall" == true ]; then
 # Methylation calling
 echo "METH. CALLING"

 methoutput=$outputdir/07_methyl_calls/
 mkdir -p $methoutput/$id/
 for chromosome in $(cat /fast/projects/peifer_wgs/work/2017-12-19_WGBS/Base/Genomes/hg19/chroms.txt | cut -c 2-); do 
   $myprogs/Rscript --vanilla $rscriptspath/methCall.R --inBam=$dedupdir/$chromosome'_dedup.sorted'.bam --assembly=hg19 --mincov=10  --minqual=20 --rds=$methoutput/$id/$chromosome.RDS --logFile=$methoutput/$id/$chromosome.log
 done
fi 

#-------------------------------
# TODO: put output of meth.calling all together, stich chromosome together



#-------------------------------
if [ "$bigwig" == true ]; then
 # Methylation calling
# echo "Save as BigWig"

# bwoutput=$outputdir/08_bigwig_files/$id/
# bwoutfile=$bwoutput/$id.bw
# mkdir -p $bwoutput/$id
# $myprogs/Rscript --vanilla $rscriptspath/export_bw.R $methoutput/$id/$id.RDS $chromsizes $assembly $bwoutfile > $bwoutput/$id.log 2> $bwoutput/$id.err

fi 


