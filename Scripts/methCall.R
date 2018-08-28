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

## catch output and messages into log file
#out <- file(argsL$logFile, open = "wt")
#sink(out,type = "output")
#sink(out, type = "message")


# Run Functions -----------------------------------------------------------


### Methylation Calling

## load methylKit
library("methylKit")

input     <- argsL$inBam
sample_id  <- argsL$sample_id
assembly  <- argsL$assembly
mincov    <- as.numeric(argsL$mincov)
minqual   <- as.numeric(argsL$minqual)
save_folder <- argsL$save_folder
save_db <- argsL$save_db

### Extract Methylation Calls

## read bam file into methylKit object
methRawDB = processBismarkAln(location = input,
                            sample.id = sample_id,
                            assembly = assembly,
                            mincov = mincov,
                            minqual = minqual,
                            save.context = "CpG",
                            save.folder = save_folder,
                            save.db = TRUE
)
# 
# #!/bin/bash
# START=$(date +%s)
# # do something
# # start your script work here
# /home/kwreczy/programs/R-3.5.0/bin/Rscript /fast/AG_Akalin/kwreczy/NB_BIH/methCall.R --inBam=/fast/AG_Akalin/kwreczy/NB_BIH/22X3H1_chr5_merged.dedup.sorted.bam --sample_id=22X3H1 --assembly=hg19 --mincov=10 --minqual=20 --save_folder=/fast/AG_Akalin/kwreczy/NB_BIH/ --save_db=True
# 
# # your logic ends here
# END=$(date +%s)
# DIFF=$(( $END - $START ))
# echo "It took $DIFF seconds"


#qsub -V -cwd -b y  -l h_vmem=10g -pe smp 16 -N methCall1 "/home/kwreczy/programs/R-3.5.0/bin/Rscript /fast/AG_Akalin/kwreczy/NB_BIH/methCall.R --inBam=/fast/AG_Akalin/kwreczy/NB_BIH/22X3H1_chr5_merged.dedup.sorted.bam --sample_id=22X3H1 --assembly=hg19 --mincov=10 --minqual=20 --save_folder=/fast/AG_Akalin/kwreczy/NB_BIH/ --save_db=True"

# cat todownload_XXX.txt | xargs -I{} -n 1-P 20 scp kwreczy_m@med-login1.bihealth.org:{}. 








