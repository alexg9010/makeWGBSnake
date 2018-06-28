# PiGx BSseq Pipeline.
#
# Copyright © 2018 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
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


# Run Functions -----------------------------------------------------------

## Segmentation

## load methylKit
library("methylKit")

input     <- args[1]#DIR_methcall+"{sample}/{sample}_filtered.txt.bgz"
output    <- args[3]#argsL$outBed
grFile    <- args[2]#argsL$grds
pngFile   <- args[4]#argsL$png
assembly <- args[5]#argsL$png
sampleid <- args[6]#argsL$png
logfile<- args[7]#argsL$png

# a=list(input,
# output,
# grFile ,
# pngFile,
# assembly,
# sampleid,
# logfile)
# print(a)
# saveRDS(a, "~/params.RDS")


## catch output and messages into log file
# out <- file(logfile, open = "wt")
# sink(out,type = "output")
# sink(out, type = "message")

## read input file
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

if( substrRight(input, 8)==".txt.bgz" ){
  # if input is a methylRawDB object
  
  ## read input methylRaw
  methRawDB = methRead(input,
                       sampleid , 
                       assembly, 
                       dbtype='tabix')
  
  ## convert to GRanges
  methRaw.gr= as(methRawDB,"GRanges")
  ## calculate methylation score 
  mcols(methRaw.gr)$meth=100*methRaw.gr$numCs/methRaw.gr$coverage
  ##destrand
  strand(methRaw.gr) <- "*"
  ##sort 
  methRaw.gr <- sort(methRaw.gr[,"meth"]) 
}else{
  # if input is a Granges object with a column `meth` that indicates % methylation
  methRaw.gr <- readRDS(input)
}

### Segmentation of methylation profile

# Remove chromosomes that contain <= 2 ranges
# It's a requirement for methSeg(). 
methRaw.gr.per.chr = split(methRaw.gr , seqnames(methRaw.gr))
methRaw.gr.per.chr.len = sapply(methRaw.gr.per.chr, length)
methRaw.gr.per.chr.len.2.remove = names(methRaw.gr.per.chr.len[methRaw.gr.per.chr.len<=5])

if( length(methRaw.gr.per.chr.len.2.remove)>=1 ){
  methRaw.gr <- dropSeqlevels(methRaw.gr, 
                              methRaw.gr.per.chr.len.2.remove, 
                              pruning.mode="coarse")
}

png(filename = pngFile,units = "in",width = 8,height = 4.5,res=300)
res.gr = methSeg(methRaw.gr, diagnostic.plot=TRUE)
dev.off()

## Saving object
saveRDS(res.gr,file=grFile) 


### Export

## export segments to bed file
methSeg2bed(segments = res.gr,
            trackLine = paste0("track name='",sampleid,"' ",
                               "description='meth segments of ",
                               methRawDB@sample.id,
                               " mapped to ",
                               methRawDB@assembly,
                               "' itemRgb=On"),
            colramp=colorRamp(c("gray","green", "darkgreen")),
            filename = output)

