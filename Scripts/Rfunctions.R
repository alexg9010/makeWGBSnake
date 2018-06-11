       save.methylBase.per.chrom2onefile = function(inputfile, con1, chunk.size=1e6){
          tbx <- open(TabixFile(inputfile, yieldSize=chunk.size))
          while(length(res <- scanTabix(tbx)[[1]]))
              write.table(res,
                        con1,
                        append = TRUE, ###
                        quote=FALSE,
                 col.names=FALSE,
                 row.names=FALSE,
                 sep="\t")
          close(tbx)
       }
