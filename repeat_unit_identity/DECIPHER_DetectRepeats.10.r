# Usage: 
# Rscript DECIPHER_DetectRepeats.10.r K00001 10 1

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
minScore_h <- as.numeric(ARGS[2])
proc_n <- as.numeric(ARGS[3])

print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the seqs file is too big

######################## main ########################
# set filenames
seq_h <- readAAStringSet(paste(seqID, '.faa.gz', sep=''))

cat(seqID, length(seq_h), '\n')

result <- DetectRepeats(seq_h, minScore=minScore_h, processors=proc_n)
saveRDS(result, file=paste(seqID, '.tr', sep=''), compress = TRUE)