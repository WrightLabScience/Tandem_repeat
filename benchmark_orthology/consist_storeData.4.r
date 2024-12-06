#!/usr/bin/env Rscript

# Usage: 
# consist_storeData.3.r K20276 0.5 0.6 ./

# re-install DECIPHER
LIB <- "./libraries/"
dir.create(LIB)
Sys.setenv(TMPDIR=getwd())
install.packages("DECIPHER_2.29.3_11272023.tar.gz", lib=LIB)
.libPaths(c(LIB, .libPaths()))

ARGS <- commandArgs(trailingOnly = TRUE)
inputID <- ARGS[1]
PID_min <- ARGS[2]
PID_max <- ARGS[3]
outdir <- ARGS[4]

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) 

######################## functions ########################
# matrix of MSA, gap as NA
get_mt_MSA <- function(MSA_h, MSA_width, n_seqs) {
	# matrix of MSA, gap as NA
	mt_MSA <- matrix(, nrow = MSA_width, ncol = n_seqs)
	for (index_h in seq(n_seqs)) {
		mt_MSA[,index_h] <- strsplit(toString(MSA_h[index_h]),"")[[1]]
	}
	mt_MSA[mt_MSA=='-'] <- NA
	mt_MSA[!is.na(mt_MSA)] <- 0
	return(mt_MSA)
}

# pairwise consistency for positive locations ("1"), mask gaps (NA)
get_pairs_in_PID_range <- function(MSA_h, n_seqs) {
	pairs_included_i <- c()
	pairs_included_j <- c()
	mt_PID <- 1-DistanceMatrix(MSA_h, includeTerminalGaps = TRUE, penalizeGapLetterMatches = TRUE, verbos=FALSE)
	# look at upper right triangle
	for (seq_i in seq(n_seqs-1)) { # row
		for (seq_j in c((seq_i+1):(n_seqs))) { # col
			PID_h <- mt_PID[seq_i,seq_j]
			if (PID_h >= PID_min && PID_h <= PID_max) { # save pair if in range
				pairs_included_i <- c(pairs_included_i, seq_i)
				pairs_included_j <- c(pairs_included_j, seq_j)
			}
		}
	}
	if (length(pairs_included_i) > 0) {
		mt_pairs_included <- rbind(pairs_included_i, pairs_included_j)
		return(mt_pairs_included)
	} else {
		return(NULL)
	}
}

######################## main ########################

# 'InputID', 'n_seqs', 'aa_MSA_width', 'n_pairs'

suppressMessages(library(DECIPHER))

# read aa MSA
input_msa <- readRDS(paste(inputID, '.aa.msa', sep=''))
MSA_width <- width(input_msa)[1]
n_seqs <- length(input_msa) # number of seq in InputID
mt_pairs_included <- get_pairs_in_PID_range(input_msa, n_seqs)	
if (is.null(mt_pairs_included)) {
    cat('##Report:', inputID, n_seqs, MSA_width, 0, '\n', sep='\t')
} else {
    cat('##Report:', inputID, n_seqs, MSA_width, ncol(mt_pairs_included), '\n', sep='\t')
    saveRDS(mt_pairs_included, file=paste(outdir, inputID, '_pairs.rds', sep=''), compress = TRUE)
    saveRDS(get_mt_MSA(input_msa, MSA_width, n_seqs), file=paste(outdir, inputID, '_mt_MSA_aa.rds', sep=''), compress = TRUE)
    # read nt MSA
    input_msa <- readRDS(paste(inputID, '.nt.msa', sep=''))
    MSA_width <- width(input_msa)[1]
    n_seqs <- length(input_msa) # number of seq in InputID
    saveRDS(get_mt_MSA(input_msa, MSA_width, n_seqs), file=paste(outdir, inputID, '_mt_MSA_nt.rds', sep=''), compress = TRUE)
}
cat('DONE', inputID, '\n')