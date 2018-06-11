
args <- commandArgs(TRUE)

# saveRDS(args, "~/args.RDS")
# args = readRDS("~/args.RDS")

inputfiles = args[1]
outfile = args[2]

inputs <- strsplit(inputfiles, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

# this function works only per chromosome and I want to merge 
# files that have different chromosome so it won't work for this case
# library(methylKit)
# methylKit:::mergeTabix(tabixList=inputs,
#                        dir=outdir,
#                        filename=outfilename,
#                        mc.cores=as.numeric("{params.cores}"),
#                        all=FALSE)

library(Rsamtools)

save.methylBase.per.chrom2onefile = function(inputfile, con1, chunk.size=1e6){
  tbx <- open(TabixFile(inputfile, yieldSize=chunk.size))
  while(length(res <- scanTabix(tbx)[[1]]))
    write.table(res,
                con1,
                append = TRUE,
                quote=FALSE,
                col.names=FALSE,
                row.names=FALSE,
                sep="\t")
  close(tbx)
}

con <- file(outfile, open="a")
for(infile in inputs)
  save.methylBase.per.chrom2onefile(infile, con, chunk.size=100)
close(con)

# bgzip and index
methylKit:::makeMethTabix(outfile)



