

# It can be done only on beast right now, 
# ebcause I cant install the newest methylKit on BIH cluster right now
# It's due to the fact that even if I can install the newest R from conda
# I cant install the newwest methylKit via conda (they have only and old 1.4..sth version of methylKit)

BENDIR = "/fast/users/kwreczy_m/projects/makeWGBSnake/Test/output/benchmarks/"

benfiles = list.files(path = BENDIR, pattern = NULL, all.files = FALSE,
           full.names = TRUE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


rule.names = sapply(1:length(benfiles), function(i) strsplit(benfiles[i], "[.]")[[1]][2] )

max_vms = sapply(1:length(benfiles), function(i){
  tbl = read.table(benfiles[i], header=TRUE)
  tbl$max_vms
})
secs = sapply(1:length(benfiles), function(i){
  tbl = read.table(benfiles[i], header=TRUE)
  tbl$s
})

colnames(max_vms ) = rule.names
colnames(secs ) = rule.names

library(reshape)
library(ggplot2)

pdf("~/boxplots_VM.pdf")

plot.data=melt(max_vms)
q=ggplot(plot.data, aes(x=X2, y=value, fill=X2)) +  # This is the plot function
  geom_boxplot()      # This is the geom for box plot in ggplot.

q + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Rules") + ylab("VM [Mb]") + guides(fill=FALSE)

dev.off()
  


pdf("~/boxplots_Time.pdf")

plot.data=melt(secs)
q=ggplot(plot.data, aes(x=X2, y=value, fill=X2)) +  # This is the plot function
  geom_boxplot()      # This is the geom for box plot in ggplot.

q + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("Rules") + ylab("Execution time [sec]") + guides(fill=FALSE)

dev.off()





