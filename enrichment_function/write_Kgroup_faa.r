# Usage: 
# Rscript write_Kgroup_faa.r K00001

ARGS <- commandArgs(trailingOnly = TRUE)
Kgroup_h <- ARGS[1] # Kgroup

# load libraries
library(Biostrings)

######################## main ########################
count <- 0
dir_list <- list.dirs('Kgroup_seq',recursive = FALSE) 
for (genomeID_h in dir_list) {
	seq_h <- NULL
	filename_h <- paste(genomeID_h, '/', Kgroup_h, '.faa.gz', sep='')
	try(seq_h <- readAAStringSet(filename_h), silent = TRUE)
	if (!is.null(seq_h)) { # if Kgroup_h is in genomeID_h
		count <- count + 1
		if (count == 1) { # the first seq
			Kgroup_seq_h <- seq_h
		} else {
			Kgroup_seq_h <- c(Kgroup_seq_h, seq_h)
		}
	}
}
filename_h <- paste(Kgroup_h, '.faa.gz', sep='')
writeXStringSet(Kgroup_seq_h, filename_h, append=FALSE, compress=TRUE, format="fasta")
cat('\nDone write faa.\n')