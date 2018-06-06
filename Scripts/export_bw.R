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
})

args <- commandArgs(trailingOnly = TRUE)

tabix_filepath    <- args[1]
seqlengths_path <- args[2]
assembly        <- args[3]
sampleid        <- args[4]
out_path        <- args[5]

m1 = methRead(tabix_filepath,   #import the methylRaw object from tabix file.
             sampleid, 
             assembly, 
             dbtype='tabix')

seqdat_temp = read.table(seqlengths_path, sep="\t", header=FALSE, col.names=c("chr", "len"), stringsAsFactors = FALSE)
Sinfo <- Seqinfo(seqnames   = seqdat_temp$chr,
                 seqlengths = seqdat_temp$len,
                 genome     = assembly)


G1            <- as(m1 , "GRanges")            # convert it to a GRanges object

seqlevels(G1) <- seqlevels(Sinfo)              # ensure the full set of seqnames 
                                               # from the ref-genome are included
                                               # (even if this data set is low-
                                               # coverage and missing chrom's)
seqinfo(G1)   <- Sinfo
G1$score = G1$numCs/G1$coverage

G1$coverage = NULL
G1$numCs    = NULL
G1$numTs    = NULL

export.bw( object = G1, con=out_path  )


