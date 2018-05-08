

suppressPackageStartupMessages(expr = {
  library(GenomicRanges)
  library(stringr)
  library(methylKit)
  library(rtracklayer)
})

args <- commandArgs(trailingOnly = TRUE)

RDS_filepath    <- args[1]
outputdir <- args[2]


methylBase.obj = readRDS(RDS_filepath)


makeMethylDB(methylBase.obj, dbdir=outputdir)



