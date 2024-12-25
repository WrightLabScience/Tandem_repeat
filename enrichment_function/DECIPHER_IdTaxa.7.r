# Usage: 
# Rscript DECIPHER_IdTaxa.7.r GCA_932276165 KEGG_Animals_r95.RData

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
genomeID_h <- ARGS[1]
TRS <- ARGS[2] # training set

print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# load training set
load(TRS, verbose = TRUE)
cat("\nInput training set:", TRS, '\n')

# set filenames
faa_filename <- paste(genomeID_h, '.faa.gz', sep='')
out_filename <- paste(genomeID_h, '.id', sep='')

record_len <- 100
cat('\n', 'read in at most', record_len,'seqs.\n')

end <- FALSE
start_i <- 1
while (end!=TRUE) {
	if (start_i==1) {
		# read in at most record_len seqs, starts from start_i (skip start_i-1)
		seq_h <- readAAStringSet(faa_filename, nrec=record_len)
		seq_h <- RemoveGaps(seq_h, "all")
		cat('\nRunning IdTaxa on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
		result_all <- IdTaxa(test = seq_h, trainingSet = trainingSet, threshold = 40, fullLength = 0.99, verbose = TRUE)
		start_i <- start_i + length(seq_h) # next round should starts from 101
	} else {
		# read in at most record_len seqs, starts from start_i
		seq_h <- readAAStringSet(faa_filename, nrec=record_len, skip=start_i-1)
		seq_h <- RemoveGaps(seq_h, "all")
		if (length(seq_h)!=0) { # if there is seq
			cat('\nRunning IdTaxa on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
			result_h <- IdTaxa(test = seq_h, trainingSet = trainingSet, threshold = 40, fullLength = 0.99, verbose = TRUE)
            result_all <- c(result_all, result_h)
			start_i <- start_i + length(seq_h) # next round should starts from
		} else {
			end <- TRUE
		}
	}
}
cat('\nDONE running IdTaxa on', length(result_all), 'CDS...\n\n')
# write output
saveRDS(result_all, file=out_filename, compress = TRUE)
cat('\nSave IdTaxa result as', out_filename, '\n')