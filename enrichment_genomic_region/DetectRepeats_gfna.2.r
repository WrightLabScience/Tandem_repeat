# Usage: 
# Rscript DetectRepeats_gfna.2.r GCF_028858775 3

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
contig_index <- as.numeric(ARGS[2])

print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

########################################## START OF SCRIPT ##########################################

# set filenames
in_filename <- paste(seqID, '.fna.gz', sep='')
out_filename <- paste(seqID, '_', contig_index, '.tr', sep='')

seq_h <- readDNAStringSet(in_filename, nrec=1, skip=contig_index-1)
if (length(seq_h)>0) { # if there is seq
	cat('\nRunning DetectRepeats on contig No.', contig_index, ":", names(seq_h), width(seq_h), '...\n\n')
	TR_h <- DetectRepeats(seq_h)
	if (nrow(TR_h)>0) {
		TR_h$Index <- contig_index
		saveRDS(TR_h, file=out_filename, compress = TRUE)
	}
	cat('#### DONE ####', nrow(TR_h), 'TR detected in', seqID, 'contig No.', contig_index, '\n')
} else {
	cat('!!!!!There is no contig No.', contig_index, 'in', seqID, '\n')
}