#!/usr/bin/env Rscript

# Usage: 
# Rscript DECIPHER_DetectRepeats.20.r infile out_wE out_woE nt

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
infile <- ARGS[1]
out_wE <- ARGS[2]
out_woE <- ARGS[3]
ntaa <- ARGS[4]

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# set filenames
if (ntaa == 'nt') {
    seq <- readDNAStringSet(infile)
    cat("# Input nt,",infile, ", number of seqs:", length(seq), '\n')
} else {
    seq <- readAAStringSet(infile)
    cat("# Input aa,",infile, ", number of seqs:", length(seq), '\n')
}

tr <- DetectRepeats(seq, minScore = 0, useEmpirical=TRUE)
saveRDS(tr, file=out_wE, compress = TRUE)

tr <- DetectRepeats(seq, minScore = 0, useEmpirical=FALSE)
saveRDS(tr, file=out_woE, compress = TRUE)