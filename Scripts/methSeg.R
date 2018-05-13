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

input     <- args[1]#argsL$rds
output    <- args[3]#argsL$outBed
grFile    <- args[2]#argsL$grds
pngFile   <- args[4]#argsL$png

print(args)

## read input methylRaw
methRaw <- readRDS(input)

## convert to GRanges
methRaw.gr= as(methRaw,"GRanges")
## calculate methylation score 
mcols(methRaw.gr)$meth=100*methRaw.gr$numCs/methRaw.gr$coverage
##destrand
strand(methRaw.gr) <- "*"
##sort 
methRaw.gr <- sort(methRaw.gr[,"meth"]) 


### Segmentation of methylation profile

png(filename = pngFile,units = "in",width = 8,height = 4.5,res=300)
res.gr = methSeg(methRaw.gr,diagnostic.plot=TRUE)
dev.off()

## Saving object
saveRDS(res.gr,file=grFile) 


### Export

## export segments to bed file
methSeg2bed(segments = res.gr,
            trackLine = paste0("track name='meth segments ' ",
                               "description='meth segments of ",
                               methRaw@sample.id,
                               " mapped to ",
                               methRaw@assembly,
                               "' itemRgb=On"),
            colramp=colorRamp(c("gray","green", "darkgreen")),
            filename = output)

