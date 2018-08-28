
##################################### my.methSeg START

#' @param estimate.params.density a numeric value indicating percentage of regions to sample
#' and then estimate initial parameters (G and modelName) for a function to calculate density 
#' (mclust::densityMclust()). The value can be between 0 and 1, e.g. 0.1 means that 10% of
#' data will be used for sampling estimate initial parameters, otherwise it uses the whole
#' dataset (Default: 1) and it's time consuming on large datasets. If 0 or 1 then the function
#' will be executed without sampling. 
#' 
my.methSeg<-function(obj, diagnostic.plot=TRUE, join.neighbours=FALSE,
                     estimate.params.density=1, ...){
  
  
  require(fastseg)
  require(mclust)
  densityFind <- methylKit:::densityFind
  joinSegmentNeighbours <- methylKit:::joinSegmentNeighbours
  
  dots <- list(...)  
  
  
  ##coerce object to granges
  if(class(obj)=="methylRaw" | class(obj)=="methylRawDB") {
    obj= as(obj,"GRanges")
    ## calculate methylation score 
    mcols(obj)$meth=100*obj$numCs/obj$coverage
    ## select only required mcol
    obj = obj[,"meth"]
  }else if (class(obj)=="methylDiff" | class(obj)=="methylDiffDB") {
    obj = as(obj,"GRanges")
    ## use methylation difference as score
    obj = obj[,"meth.diff"]
  }else if (class(obj) != "GRanges"){
    stop("only methylRaw or methylDiff objects ", 
         "or GRanges objects can be used in this function")
  }
  
  # destrand
  strand(obj) <- "*"
  
  ## check wether obj contains at least one metacol 
  if(ncol(elementMetadata(obj))<1)
    stop("GRanges does not have any meta column.")
  
  ## check wether obj contains is sorted by position
  if(is.unsorted(obj,ignore.strand=TRUE)) {
    obj <- sort(obj,ignore.strand=TRUE)
    message("Object not sorted by position, sorting now.")
  }
  
  ## check wether obj contains at least two ranges else stop
  if(length(obj)<=1)
    stop("segmentation requires at least two ranges.")
  
  # match argument names to fastseg arguments
  args.fastseg=dots[names(dots) %in% names(formals(fastseg)[-1] ) ]  
  
  # match argument names to Mclust
  args.Mclust=dots[names(dots) %in% names(formals(Mclust)[-1])  ]
  
  args.fastseg[["x"]]=obj
  
  # do the segmentation
  #seg.res=fastseg(obj)
  seg.res <- do.call("fastseg", args.fastseg)
  #seg.res <- do.call("fastseg2", args.fastseg)
  
  # stop if segmentation produced only one range
  if(length(seg.res)==1) {
    warning("segmentation produced only one range, no mixture modeling possible.")
    seg.res$seg.group <- "1"
    return(seg.res)
  }
  
  # if joining, do not show first clustering
  if(join.neighbours) {
    diagnostic.plot.old = diagnostic.plot
    diagnostic.plot = FALSE
  }
  
  if(estimate.params.density>0 & estimate.params.density<1 ){
    
    nbr.sample = floor(length(seg.res) * estimate.params.density)
    
    # estimate parameters for mclust
    # finds the optimal number of componets from a sampled data
    # and then use it as starting points
    args.Mclust.esti = args.Mclust
    args.Mclust.esti[["score.gr"]]=seg.res[ sample(1:length(seg.res), nbr.sample) ]
    print(args.Mclust.esti[["score.gr"]])
    args.Mclust.esti[["diagnostic.plot"]]=FALSE
    dens_estimate=do.call("densityFind", args.Mclust.esti)
    
    args.Mclust[["modelName"]] = dens_estimate$modelName
    args.Mclust[["G"]] = 1:dens_estimate$G
    print(dens_estimate$G)
  }
  
   # decide on number of components/groups
    args.Mclust[["score.gr"]]=seg.res
    args.Mclust[["diagnostic.plot"]]=diagnostic.plot
    dens=do.call("densityFind", args.Mclust  )

  
  # add components/group ids 
  mcols(seg.res)$seg.group=as.character(dens$classification)
  
  # if joining, show clustering after joining
  if(join.neighbours) {
    message("joining neighbouring segments and repeating clustering.")
    seg.res <- joinSegmentNeighbours(seg.res)
    diagnostic.plot <- diagnostic.plot.old
    
    # get the new density
    args.Mclust[["score.gr"]]=seg.res
    args.Mclust[["diagnostic.plot"]]=diagnostic.plot
    # skip second progress bar
    args.Mclust[["verbose"]]=FALSE
    dens=do.call("densityFind", args.Mclust  )
    
  }
  
  list(seg.res=seg.res,
       dens=dens,
       dens_estimate=dens_estimate)
}

##################################### my.methSeg END



###################################### methylkit methSeg START

fastermethSeg<-function(obj, diagnostic.plot=TRUE, join.neighbours=FALSE,
                  initialize.on.subset=1, ...){
  
  
  require(fastseg)
  require(mclust)
  densityFind <- methylKit:::densityFind
  joinSegmentNeighbours <- methylKit:::joinSegmentNeighbours
  
  dots <- list(...)  
  
  
  ##coerce object to granges
  if(class(obj)=="methylRaw" | class(obj)=="methylRawDB") {
    obj= as(obj,"GRanges")
    ## calculate methylation score 
    mcols(obj)$meth=100*obj$numCs/obj$coverage
    ## select only required mcol
    obj = obj[,"meth"]
  }else if (class(obj)=="methylDiff" | class(obj)=="methylDiffDB") {
    obj = as(obj,"GRanges")
    ## use methylation difference as score
    obj = obj[,"meth.diff"]
  }else if (class(obj) != "GRanges"){
    stop("only methylRaw or methylDiff objects ", 
         "or GRanges objects can be used in this function")
  }
  
  # destrand
  strand(obj) <- "*"
  
  ## check wether obj contains at least one metacol 
  if(ncol(elementMetadata(obj))<1)
    stop("GRanges does not have any meta column.")
  
  ## check wether obj contains is sorted by position
  if(is.unsorted(obj,ignore.strand=TRUE)) {
    obj <- sort(obj,ignore.strand=TRUE)
    message("Object not sorted by position, sorting now.")
  }
  
  ## check wether obj contains at least two ranges else stop
  if(length(obj)<=1)
    stop("segmentation requires at least two ranges.")
  
  # match argument names to fastseg arguments
  args.fastseg=dots[names(dots) %in% names(formals(fastseg)[-1] ) ]  
  
  # match argument names to Mclust
  args.Mclust=dots[names(dots) %in% names(formals(Mclust)[-1])  ]
  
  args.fastseg[["x"]]=obj
  
  # do the segmentation
  #seg.res=fastseg(obj)
  seg.res <- do.call("fastseg", args.fastseg)
  #seg.res <- do.call("fastseg2", args.fastseg)
  
  # stop if segmentation produced only one range
  if(length(seg.res)==1) {
    warning("segmentation produced only one range, no mixture modeling possible.")
    seg.res$seg.group <- "1"
    return(seg.res)
  }
  
  # if joining, do not show first clustering
  if(join.neighbours) {
    diagnostic.plot.old = diagnostic.plot
    diagnostic.plot = FALSE
  }
  
  if("initialization" %in% names(args.Mclust)){
    if("subset" %in% names(args.Mclust[["initialization"]])) {
      if(length(args.Mclust[["initialization"]][["subset"]]) < 9 ){
        stop("too few samples, increase the size of subset.") 
      }
      message(paste("initializing clustering with",
                    length(args.Mclust[["initialization"]][["subset"]]),
                    "segments."))
      initialize.on.subset = 1
    }
  }
  
  if(initialize.on.subset != 1 && initialize.on.subset > 0 ) {
    
    if( initialize.on.subset > 0 & initialize.on.subset < 1 )
      nbr.sample = floor(length(seg.res) * initialize.on.subset)
    
    if( initialize.on.subset > 1) 
      nbr.sample = initialize.on.subset
    
    if( nbr.sample < 9 ){stop("too few samples, increase the size of subset.") }
    
    message(paste("initializing clustering with",nbr.sample,"out of",length(seg.res),"total segments."))
    # estimate parameters for mclust
    sub <- sample(1:length(seg.res), nbr.sample,replace = FALSE)
    args.Mclust[["initialization"]]=list(subset = sub)
  }  
  
  
  # decide on number of components/groups
  args.Mclust[["score.gr"]]=seg.res
  args.Mclust[["diagnostic.plot"]]=diagnostic.plot
  dens=do.call("densityFind", args.Mclust  )
  
  
  # add components/group ids 
  mcols(seg.res)$seg.group=as.character(dens$classification)
  
  # if joining, show clustering after joining
  if(join.neighbours) {
    message("joining neighbouring segments and repeating clustering.")
    seg.res <- joinSegmentNeighbours(seg.res)
    diagnostic.plot <- diagnostic.plot.old
    
    # get the new density
    args.Mclust[["score.gr"]]=seg.res
    args.Mclust[["diagnostic.plot"]]=diagnostic.plot
    # skip second progress bar
    args.Mclust[["verbose"]]=FALSE
    dens=do.call("densityFind", args.Mclust  )
    
  }
  
  seg.res
}

########################################################################## methylkit methSeg END

########### Benchamrk

options(scipen=999)
library(methylKit)

setwd("/fast/users/kwreczy_m/projects/Meth_analysis/")
source("./scripts/loadData.R")
source("./scripts/functions.R")

patient.ids = as.character(wgbs.table$DKFZ_PATIENT_ID)

### Read filtered methylation files in tabix
methyldir = '/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/subset_hg19/snakemake_test/06_methyl_calls/'

methylRawDB.list.obj_filtered = mclapply(patient.ids, function(x)
  methRead(paste0(methyldir,'Sampled_',x,'/','Sampled_',x,'_filtered.txt.bgz'),
           x, 
           'hg19',
           dbtype='tabix')
  ,mc.cores=24)
methylRawListDB.obj_filtered <- as(methylRawDB.list.obj_filtered, "methylRawListDB")

# Add treatment
methylRawListDB.obj_filtered@treatment = df[match(patient.ids, df$DKFZ_PATIENT_ID),]$CLIN_SUBGROUP_code # not yet implemented in methylkit

#methylRawDB.list.obj_filtered[[5]]
methylRawDB.list.obj_filtered[[14]]


mybench = mclapply(1:length(methylRawDB.list.obj_filtered), function(i){
  methRawDB = methylRawDB.list.obj_filtered[[i]]
  ## convert to GRanges
  methRaw.gr= as(methRawDB,"GRanges")
  ## calculate methylation score 
  mcols(methRaw.gr)$meth=100*methRaw.gr$numCs/methRaw.gr$coverage
  ##destrand
  strand(methRaw.gr) <- "*"
  ##sort 
  methRaw.gr <- sort(methRaw.gr[,"meth"]) 
  
  require(rbenchmark)
  benchmark(methSeg=methSeg(methRaw.gr.subset, diagnostic.plot=FALSE),
                 fastermethSeg=fastermethSeg(methRaw.gr.subset, diagnostic.plot=FALSE, initialize.on.subset=.2),
                 kasia.fastermethSeg=my.methSeg(methRaw.gr.subset, diagnostic.plot=FALSE, estimate.params.density=.2),
                 replications=100,
                 columns=c('test', 'replications', 'elapsed'))
  
}, mc.cores=24)

a=do.call("rbind",lapply(1:length(mybench), function(i) mybench[[i]]$elapsed))
colnames(a) = c( "fastermethSeg", "kasia.fastermethSeg","methSeg")

  
  
png("~/fastermethSeg.png")
boxplot(a,names=colnames(a), ylab="Seconds")
dev.off()



# require(rbenchmark)
# b <- benchmark(methSeg.subsetnbdata=methSeg(methRaw.gr.subset, diagnostic.plot=FALSE),
#                fastermethSeg.subsetnbdata=fastermethSeg(methRaw.gr.subset, diagnostic.plot=FALSE, initialize.on.subset=.2),
#                kasia.fastermethSeg.subsetnbdata=my.methSeg(methRaw.gr.subset, diagnostic.plot=FALSE, estimate.params.density=.2),
#                replications=100,
#                columns=c('test', 'replications', 'elapsed'))
# 
# e <- benchmark(methSeg.subsetnbdata=methSeg(methRaw.gr.subset, diagnostic.plot=FALSE),
#                fastermethSeg.subsetnbdata=fastermethSeg(methRaw.gr.subset, diagnostic.plot=FALSE, initialize.on.subset=.2),
#                kasia.fastermethSeg.subsetnbdata=my.methSeg(methRaw.gr.subset, diagnostic.plot=FALSE, estimate.params.density=.2),
#                replications=100,
#                columns=c('test', 'replications', 'elapsed'))
# > e
# test replications elapsed
# 2       fastermethSeg.subsetnbdata          100  84.886
# 3 kasia.fastermethSeg.subsetnbdata          100  79.418
# 1             methSeg.subsetnbdata          100  82.865
# dd <- benchmark(methSeg.subsetnbdata=methSeg(methRaw.gr1, diagnostic.plot=FALSE),
#                fastermethSeg.subsetnbdata=fastermethSeg(methRaw.gr1, diagnostic.plot=FALSE, initialize.on.subset=.2),
#                kasia.fastermethSeg.subsetnbdata=my.methSeg(methRaw.gr1, diagnostic.plot=FALSE, estimate.params.density=.2),
#                replications=100,
#                columns=c('test', 'replications', 'elapsed'))


library(rbenchmark)
d <- benchmark(methSeg.subsetnbdata=methSeg(methRaw.gr.wholedata, diagnostic.plot=FALSE),
               fastermethSeg.subsetnbdata=fastermethSeg(methRaw.gr.wholedata, diagnostic.plot=FALSE, initialize.on.subset=.2),
               kasia.fastermethSeg.subsetnbdata=my.methSeg(methRaw.gr.wholedata, diagnostic.plot=FALSE, estimate.params.density=.2),
               replications=1,
               columns=c('test', 'replications', 'elapsed'))



##########


## Collect arguments
args <- commandArgs(TRUE)


# Run Functions -----------------------------------------------------------

## Segmentation

## load methylKit
library("methylKit")

input     <- args[1]#DIR_methcall+"{sample}/{sample}_filtered.txt.bgz" or methylRaw_filtered_GRanges_",p,".RDS"
output    <- args[3]#argsL$outBed
grFile    <- args[2]#argsL$grds
pngFile   <- args[4]#argsL$png
assembly <- args[5]#argsL$png
sampleid <- args[6]#argsL$png
logfile<- args[7]#argsL$png

segdir = "/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/hg19/bsseq_pigx/08_segmentation/"
sampleid = "ZTFVN7"

input <- paste0("/fast/projects/peifer_wgs/work/2017-12-19_WGBS/Project/Results/RDS//methylRaw_filtered_GRanges_",sampleid,".RDS")
output <- paste0(segdir,sampleid,".deduped_meth_segments_gr.RDS")
grFile<- paste0(segdir,sampleid,".deduped_meth_segments_gr.RDS")
pngFile   <- paste0(segdir,sampleid,".deduped_meth_segments.png")
assembly <- "hg19"
sampleid <- sampleid
logfile<-paste0(segdir,sampleid,".deduped_meth_segments.log")


## catch output and messages into log file
out <- file(logfile, open = "wt")
sink(out,type = "output")
sink(out, type = "message")


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

if( substrRight(input, 7)==".txt.bgz" ){
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
  # if inut is a Granges object with a column `meth` that indicates % methylation
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
res.gr = my.methSeg(methRaw.gr, diagnostic.plot=TRUE)
dev.off()

## Saving object
saveRDS(res.gr,file=grFile) 


### Export

## export segments to bed file
methSeg2bed(segments = res.gr$seg.res,
            trackLine = paste0("track name='meth segments ' ",
                               "description='meth segments of ",
                               methRawDB@sample.id,
                               " mapped to ",
                               methRawDB@assembly,
                               "' itemRgb=On"),
            colramp=colorRamp(c("gray","green", "darkgreen")),
            filename = output)


###############################

mybench.length = mclapply(1:length(methylRawDB.list.obj_filtered), function(i){
  methRawDB = methylRawDB.list.obj_filtered[[i]]
  ## convert to GRanges
  methRaw.gr= as(methRawDB,"GRanges")
  ## calculate methylation score 
  mcols(methRaw.gr)$meth=100*methRaw.gr$numCs/methRaw.gr$coverage
  ##destrand
  strand(methRaw.gr) <- "*"
  ##sort 
  methRaw.gr <- sort(methRaw.gr[,"meth"]) 
  
  length(methRaw.gr)
  
  # require(rbenchmark)
  # benchmark(methSeg=methSeg(methRaw.gr.subset, diagnostic.plot=FALSE),
  #           fastermethSeg=fastermethSeg(methRaw.gr.subset, diagnostic.plot=FALSE, initialize.on.subset=.2),
  #           kasia.fastermethSeg=my.methSeg(methRaw.gr.subset, diagnostic.plot=FALSE, estimate.params.density=.2),
  #           replications=100,
  #           columns=c('test', 'replications', 'elapsed'))
  
}, mc.cores=24)

> mybench

> as.numeric(unlist(mybench.length))
[1]  1127    NA  1495  1490 21073  2186  1152  1682  1578  2404  2334   526
[13]  2589 13003  3807  1089  2084  1006 10013   874  1835  1859   750   520


[[1]]                                                                                                                                              [94/1812]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 120.197
3 kasia.fastermethSeg.subsetnbdata          100 114.069
1             methSeg.subsetnbdata          100 117.504

[[2]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100  79.189
3 kasia.fastermethSeg.subsetnbdata          100  76.014
1             methSeg.subsetnbdata          100  77.541

[[3]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 117.146
3 kasia.fastermethSeg.subsetnbdata          100 116.787
1             methSeg.subsetnbdata          100 121.180

[[4]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 120.372
3 kasia.fastermethSeg.subsetnbdata          100 116.981
1             methSeg.subsetnbdata          100 119.617

[[5]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 119.156
3 kasia.fastermethSeg.subsetnbdata          100 116.698
1             methSeg.subsetnbdata          100 121.571

[[6]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 120.859
3 kasia.fastermethSeg.subsetnbdata          100 115.658
1             methSeg.subsetnbdata          100 120.028

[[7]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 118.214
3 kasia.fastermethSeg.subsetnbdata          100 118.218
1             methSeg.subsetnbdata          100 118.998

[[8]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 116.458
3 kasia.fastermethSeg.subsetnbdata          100 116.242
1             methSeg.subsetnbdata          100 115.838

[[9]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.879
3 kasia.fastermethSeg.subsetnbdata          100 115.901
1             methSeg.subsetnbdata          100 116.943

[[10]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.169
3 kasia.fastermethSeg.subsetnbdata          100 116.747
1             methSeg.subsetnbdata          100 118.725

[[11]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 116.428
3 kasia.fastermethSeg.subsetnbdata          100 117.599
1             methSeg.subsetnbdata          100 117.193

[[12]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 119.771
3 kasia.fastermethSeg.subsetnbdata          100 116.913
1             methSeg.subsetnbdata          100 120.611

[[13]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 119.893
3 kasia.fastermethSeg.subsetnbdata          100 116.902
1             methSeg.subsetnbdata          100 117.231

[[14]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.275
3 kasia.fastermethSeg.subsetnbdata          100 115.392
1             methSeg.subsetnbdata          100 120.103

[[15]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 117.097
3 kasia.fastermethSeg.subsetnbdata          100 116.780
1             methSeg.subsetnbdata          100 116.537

[[16]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 118.714
3 kasia.fastermethSeg.subsetnbdata          100 116.790
1             methSeg.subsetnbdata          100 115.705


[[17]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 120.306
3 kasia.fastermethSeg.subsetnbdata          100 117.272
1             methSeg.subsetnbdata          100 119.660

[[18]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.575
3 kasia.fastermethSeg.subsetnbdata          100 118.130
1             methSeg.subsetnbdata          100 118.685

[[19]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 120.476
3 kasia.fastermethSeg.subsetnbdata          100 116.689
1             methSeg.subsetnbdata          100 119.181

[[20]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.107
3 kasia.fastermethSeg.subsetnbdata          100 116.795
1             methSeg.subsetnbdata          100 117.275

[[21]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 118.310
3 kasia.fastermethSeg.subsetnbdata          100 116.707
1             methSeg.subsetnbdata          100 119.764

[[22]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.267
3 kasia.fastermethSeg.subsetnbdata          100 117.338
1             methSeg.subsetnbdata          100 115.733

[[23]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 118.974
3 kasia.fastermethSeg.subsetnbdata          100 115.635
1             methSeg.subsetnbdata          100 120.401

[[24]]
test replications elapsed
2       fastermethSeg.subsetnbdata          100 121.855
3 kasia.fastermethSeg.subsetnbdata          100 113.872
1             methSeg.subsetnbdata          100 114.913

