

args <- commandArgs(TRUE)

library("methylKit")

saveRDS(args, "~/arg.RDS")

output <- args[1]
treatment <- args[2]
mincov <- args[3]
input     <- args[4:length(args)]

print(output)
print(class(input))


#-------------------------------------------------------------------
#' Combine multiple methylRaw objects into a methylRawList object
#'
#' @param list.of.methylRaw a list of methylRaw objects from the methylKit package
#' @param treatment a numeric vector indicating treaments
#' @param min.cov a number indicating a minimum coverage
combine2methylRawList <- function(list.of.methylRaw, treatment, min.cov) {

  ## check if treatment has same length as number of inputs
  if(length(list.of.methylRaw)!=length(treatment))
    stop("Treatment vector doesnt have the same length as list of methylRaw objects.")

  ## check if input is really of type methylRaw
  if(!all(sapply(list.of.methylRaw, function(x) class(x)=="methylRaw")))
    stop("Input objects are not methylRaw objects.")

  ## remove data beyond min coverage
  list.of.methylRaw.mincov = lapply(1:length(list.of.methylRaw), function(i){
    #which.gtmincov = which( getData(rds.objs[[i]])$coverage >= min.cov )
    which.gtmincov = which( rds.objs[[i]][[5]] >= min.cov )
    rds.objs[[i]][which.gtmincov,]
    })

  ## merge
  mrl <- new("methylRawList", list.of.methylRaw.mincov)
  mrl@treatment <- treatment
  mrl
}


stich_list_of_methylRaw2methylRaw = function(list_methylRaw,
                             patientid, 
                             cores=20){
  require(methylKit)
  require(data.table)
  
  a=mclapply(rdsfiles, readRDS, mc.cores=cores)
  my.df.list = mclapply(a, function(x) data.table(getData(x)), mc.cores=cores)

  # adata = do.call("rbind", my.df.list) # too slow, use data.table instead
  dt <- rbindlist(my.df.list)
  
  dt.ordered=dt[order(chr,start,decreasing=FALSE),]
  
  obj=new("methylRaw", dt.ordered, 
                       sample.id=patientid,
                       assembly='hg19', ########### TODO:
                       context='CpG',
                       resolution='base')
  
  return( obj )
  
}


rds.objs = lapply(input, readRDS)
methRawList.obj = stich_list_of_methylRaw2methylRaw(rds.objs, treatment, mincov)


# TODOOOOOOOOOOO
saveRDS(methRawList.obj, file=output)








