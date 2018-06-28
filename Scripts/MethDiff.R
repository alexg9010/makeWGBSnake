# Copyright © 2018 Katarzyna Wreczycka <Katarzyna.Wreczycka@mdc-berlin.de>

args <- commandArgs(TRUE)


library(methylKit)

input = args[1]
output = args[2]
t = as.numeric( args[3] )
cores = as.numeric( args[4] )
treatments = args[5]
sampleids = args[6]
context = args[7]
assembly = args[8]
outputdir = args[9]
suffix = args[10]
save.db = args[11]

treatments <- as.numeric(strsplit(treatments, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]])
sampleids <- strsplit(sampleids, " ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]]

# read RDS
#my.methylBase.obj = readRDS(input)
# read a tabix file
myMethylBaseDB <- methylKit:::readMethylBaseDB(dbpath = input,
                                    dbtype = "tabix",
                                    sample.ids = sampleids,
                                    assembly = assembly,
                                    context = context,
                                    resolution = "base",
                                    treatment = treatments,
                                    destranded = FALSE)


# change treatment vector
myMethylBaseDB@treatment  = ifelse(myMethylBaseDB@treatment==t,  1, 0)

myDiff<-calculateDiffMeth(myMethylBaseDB,
                          overdispersion="MN",
                          test="Chisq",
                          mc.cores=cores,
                          save.db = TRUE,
                          dbdir=outputdir,
                          suffix = suffix
                          )

saveRDS(myDiff, output)   
