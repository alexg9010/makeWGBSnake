# PiGx BSseq Pipeline.
#
# Copyright © 2017, 2018 Bren Osberg <Brendan.Osberg@mdc-berlin.de>
# Copyright © 2017 Alexander Gosdschan <alexander.gosdschan@mdc-berlin.de>
# Copyright © 2017, 2018 Katarzyna Wreczycka <katwre@gmail.com>
# Copyright © 2017, 2018 Ricardo Wurmus <ricardo.wurmus@mdc-berlin.de>
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


suppressPackageStartupMessages(expr = {
  library(GenomicRanges)
  library(stringr)
  library(methylKit)
  library(rtracklayer)
  library(data.table)
  
})

args <- commandArgs(trailingOnly = TRUE)

tabix_filepath    <- args[1]
seqlengths_path <- args[2]
assembly        <- args[3]
sampleid        <- args[4]
out_path        <- args[5]


seqdat_temp = read.table(seqlengths_path, sep="\t", header=FALSE, col.names=c("chr", "len"), stringsAsFactors = FALSE)
Sinfo <- Seqinfo(seqnames   = seqdat_temp$chr,
                 seqlengths = seqdat_temp$len,
                 genome     = assembly)

## When working really big files GRanges consume too much memory and I stuck, it just doesnt work, 
## that's why I rewrote it to data.table,
## and after filtering, calculating a score make a GRanges out of it
# m1 = methRead(tabix_filepath,   #import the methylRaw object from tabix file.
#              sampleid, 
#              assembly, 
#              dbtype='tabix')
# 
# 
# 
# G1            <- as(m1 , "GRanges")            # convert it to a GRanges object
# 
# seqlevels(G1, pruning.mode="coarse") <- seqlevels(Sinfo, pruning.mode="coarse")              # ensure the full set of seqnames 
#                                                # from the ref-genome are included
#                                                # (even if this data set is low-
#                                                # coverage and missing chrom's)
# seqinfo(G1)   <- Sinfo
# G1$score = G1$numCs/G1$coverage
# 
# G1$coverage = NULL
# G1$numCs    = NULL
# G1$numTs    = NULL
# 
# print("Exporting")
# 
# export.bw( object = G1, con=out_path  )


m1.dt = fread(paste0("zcat ", tabix_filepath))
colnames(m1.dt) = c('chr', 'start',   'end', 'strand', 'coverage', 'numCs', 'numTs')
m1.dt <- m1.dt[chr %in%   seqnames(Sinfo) ] 

big_list <-
  lapply(1:length(seqlengths(Sinfo)), function(i){
    Sinfo_i <- seqlengths(Sinfo)[i]
    print(Sinfo_i)
    m1.dt_i <-  m1.dt[ (chr %in% names(Sinfo_i)) & end <= Sinfo_i ] 
    m1.dt_i
})
m1.dt.filtered = do.call("rbind", big_list)


m1.dt.filtered$score = m1.dt.filtered$numCs/m1.dt.filtered$coverage

m1.dt.filtered$coverage = NULL
m1.dt.filtered$numCs    = NULL
m1.dt.filtered$numTs    = NULL


G1            <- as(m1.dt.filtered , "GRanges")            # convert it to a GRanges object
seqinfo(G1)   <- Sinfo

print("Exporting")

export.bw( object = G1, con=out_path  )





