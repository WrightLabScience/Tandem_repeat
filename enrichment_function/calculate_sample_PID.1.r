#!/usr/bin/env Rscript

# Usage: 
# Rscript calculate_sample_PID.1.r infile KgroupID 1000

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
infile <- ARGS[1]
KgroupID <- ARGS[2]
n_seq_sample <- as.numeric(ARGS[3])

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# set filenames
total_seqs <- readAAStringSet(infile)
n_total_seqs <- length(total_seqs)
cat(n_total_seqs, n_seq_sample, n_total_seqs>n_seq_sample, '\n')
if (n_total_seqs > 1) {
    if (n_total_seqs > n_seq_sample) {
        sampled_seqs <- sample(total_seqs, n_seq_sample)
    } else {
        sampled_seqs <- total_seqs
    }
    cat('running MSA\n')
    msa <- AlignSeqs(sampled_seqs, verbose=TRUE, proc=NULL)
    cat('running DM\n')
    dm <- DistanceMatrix(msa, verbose=TRUE, proc=NULL)
    PID <- mean(1-dm[lower.tri(dm)], na.rm=TRUE)
} else {
    PID <- 1
}

# print input filename
cat("\n\n\n\nResult Report PID:", KgroupID, PID, '\n\n\n\n\n', sep='\t')
