# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>,
# Katarzyna Wreczycka katarzyna.wreczycka@mdc-berlin.de
#
# This file is part of the PiGx BSseq Pipeline.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 1) {
  args <- c("--help")
}

## Help section
if("--help" %in% args) {
  cat("
      Render to report
      
      Arguments:
      --inBam location of input bam file
      --assembly assembly used to map the reads
      --mincov minimum coverage (default: 10)
      --minqual minimum base quality (default: 20)
      --rds name of the RDS output file
      --logFile file to print the logs to
      --help              - print this text
      
      Example:
      ./test.R --arg1=1 --arg2='output.txt' --arg3=TRUE \n\n")
  
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")

argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))

names(argsL) <- argsDF$V1
#saveRDS(argsL, "~/argsL.RDS") 

## catch output and messages into log file
# out <- file(argsL$logFile, open = "wt")
# sink(out,type = "output")
# sink(out, type = "message")

input     <- argsL$tabixfile
sample_id  <- argsL$sample_id
mincov    <- as.numeric(argsL$mincov)
hi_perc   <- as.numeric(argsL$hi_perc)
save_folder <- argsL$save_folder
assembly  <- argsL$assembly
rdsfile <- argsL$rdsfile
cores<- argsL$cores
chromcanonicalfile<-argsL$canon_chrs_file

# Extract sample id
#sample_id = strsplit(basename(input), "[.]")[[1]][1]

# Run Functions -----------------------------------------------------------


### Methylation Calling

## load methylKit
library("methylKit")

## Read Methylation Calls

methylRawDB.obj = methRead(input,
             sample.id=sample_id,
             assembly=assembly,
             dbtype='tabix')

## Get only canonical chromosomes

# read a file with info about chromosomes
df.chroms = read.table(chromcanonicalfile, sep="\t", stringsAsFactors = FALSE)
# From data.frame to GRanges:
df.chroms <- data.frame(chrom=df.chroms[,1], start=1, end=df.chroms[,2])
gr.chroms <- as(df.chroms, "GRanges")

# There can be situation when there will be no canonical chromosomses in the methylRawDB.obj
# and then selectByOverlap throws an error
# such as 'Error: scanTabix: 'chr9' not present in tabix index'
get.canon.chroms.methylRawDB.obj = function(methylRawDB.obj, gr.chroms, canonical_chromosomes, cores){
  
  n.ranges.per.chr = unlist(mclapply(canonical_chromosomes, function(chr){
      tryCatch( # there might be canonical chromosomes that are not in methylRawDB.obj object
               nrow(methylKit:::getTabixByChr(methylRawDB.obj@dbpath, #Rsamtools::countTabix
                                   chr=chr,
                                   return.type="data.table")), 
               error=function(e) 
                 return(NA)
                 )
  }, mc.cores=cores))
  can.chrs.methylRawDB = canonical_chromosomes[!(n.ranges.per.chr==0 | is.na(n.ranges.per.chr) | is.null(n.ranges.per.chr) )]
  if(length(can.chrs.methylRawDB)==0){
    warning("There are no ranges from canonical chromosomes in the methylation calling file.")
    can.chrs.methylRawDB=c() # its for testing, in real life shouldnt happen
  }
  return(can.chrs.methylRawDB)
}

can.chrs.methylRawDB = get.canon.chroms.methylRawDB.obj(methylRawDB.obj, 
                                                        gr.chroms, 
                                                        as.character(df.chroms$chrom), 
                                                        cores)
gr.chroms.canon.in.methylRawDB = gr.chroms[as.character(seqnames(gr.chroms)) %in% can.chrs.methylRawDB,]


methRaw_canon = selectByOverlap(methylRawDB.obj, gr.chroms.canon.in.methylRawDB )

## Filer an methRaw object
methRaw_canon_filtered_DB=filterByCoverage(methRaw_canon,
                                           lo.count=mincov,
                                           lo.perc=NULL,
                                           hi.count=NULL,
                                           hi.perc=hi_perc,
                                           dbdir=save_folder,
                                           dbtype="tabix",
                                           save.db = TRUE)



