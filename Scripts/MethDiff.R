
args <- commandArgs(TRUE)


library(methylKit)

input = args[1]
output = args[2]
t = as.numeric( args[3] )


my.methylBase.obj = readRDS(input)

# change treatment vector
my.methylBase.obj@treatment  = ifelse(my.methylBase.obj@treatment==t,  1, 0)

myDiff<-calculateDiffMeth(my.methylBase.obj,
                          overdispersion="MN",
                          test="Chisq",
                          mc.cores=1,
                          save.db = TRUE)

saveRDS(myDiff, output)   
