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
out <- file(argsL$logFile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")


library(methylKit)

## Load variables
inputs    <- strsplit(argsL$inputfiles, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
destrandTfile <- argsL$destrandTfile
destrandFfile <- argsL$destrandFfile
inputdir <- argsL$inputdir
samples <- strsplit(argsL$samples, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]
treatments <- as.numeric( strsplit(argsL$treatments, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]] )
cores <- as.numeric(argsL$cores)
assembly = argsL$assembly

## Read data
methylRawDB.list.obj_filtered = mclapply(1:length(inputs), function(i)
                 methRead(inputs[i], 
                          samples[i] , 
                          assembly, 
                          dbtype='tabix'), 
                 mc.cores=cores)
methylRawListDB.obj_filtered <- as(methylRawDB.list.obj_filtered, "methylRawListDB")

methylRawListDB.obj_filtered@treatment = treatments

## Unite
if( length(methylRawListDB.obj_filtered)>1 ){
  # destranded=TRUE
  meth.deT=unite(methylRawListDB.obj_filtered, destrand=TRUE, save.db = FALSE)
  # its faster to save methylBaseDB object into RDS and then to read it
  # than save it again as a list of tabix files.
  saveRDS(meth.deT, destrandTfile)
  
  # destranded=FALSE
  meth.deF=unite(methylRawListDB.obj_filtered, destrand=FALSE, save.db = FALSE)
  saveRDS(meth.deF, destrandFfile)
  
}else{ # is there is only 1 sample

  methylRawDB.2.methylBase = function(object, 
                                      sample.ids="sampleid",
                                      treatments=0, 
                                      destranded=TRUE,
                                      resolution="base"){
    
    new("methylBase",getData(object)[,-1],
        sample.ids=sample.ids,
        assembly=object@assembly,
        context=object@context,
        treatment=treatments,
        coverage.index=5,
        numCs.index=6,
        numTs.index=7,
        destranded=destranded,
        resolution=resolution
    )
  }
  # destranded=TRUE
  object = methylRawListDB.obj_filtered[[1]]
  meth.deT = methylRawDB.2.methylBase(object, 
                                      sample.ids=samples, 
                                      treatments=treatments, 
                                      destranded=TRUE)
  saveRDS(meth.deT, destrandTfile)
  
  # destranded=FALSE
  meth.deF = methylRawDB.2.methylBase(object, 
                                      sample.ids=samples, 
                                      treatments=treatments, 
                                      destranded=FALSE)
  saveRDS(meth.deF, destrandFfile)
}









