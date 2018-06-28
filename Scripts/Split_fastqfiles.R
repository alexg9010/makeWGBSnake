# Copyright © 2018 Katarzyna Wreczycka <Katarzyna.Wreczycka@mdc-berlin.de>

args <- commandArgs(TRUE)

in1 = args[1]
in2 = args[2]
DIR_inputdir_parts = args[3]
parts = as.numeric(args[4])

suffix = "_part"
cores=2

# calculate number of records
nbr.reads = as.numeric( system(paste0("zcat ",in1," | awk '{s++}END{print s/4}'"), intern = TRUE) )
# calcualte to how many parts divide the file
number_or_records_per_part = floor(nbr.reads / parts)

library("ShortRead")

split_perXM <- function(infile,
                        n=1000000, # 1M
                        outfile_basename,
                        outputdir,
                        suffix
                        ){
  
  outfile.prefix = tools::file_path_sans_ext(
                     tools::file_path_sans_ext(outfile_basename)
                     )
  number <- 0
  strm <- FastqStreamer(infile, n=n)
  repeat { 
    fq <- yield(strm) # The default size for both streams and samples is 1M records
    if (length(fq) == 0)
      break
    
    number = number+1
    ## as.character(id(fq)) # ids of reads
    outfile_final = paste0(outputdir, outfile.prefix, suffix, number, ".fq.gz")
    writeFastq(fq, outfile_final, "w")
  }
  close(strm)
}

mclapply(c(in1, in2), 
         function(x) split_perXM(x, 
                                 n=number_or_records_per_part, 
                                 basename(x), 
                                 DIR_inputdir_parts,
                                 suffix), 
         mc.cores=cores)

