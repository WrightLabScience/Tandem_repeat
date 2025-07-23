# Usage: 
# Rscript DECIPHER_DetectRepeats_batch.1.r GCA_932276165.faa.gz 1 1 100

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--in_file", type="character", default="NA", help="input file (Required)"),
	make_option("--proc_n", type="integer", default=1, help="number of processors"),
	make_option("--batch_i", type="integer", default=NA, help="batch index"),
	make_option("--record_len", type="integer", default=100, help="number of sequences per batch"),
	make_option("--minScore", type="double", default=8, help="score cutoff for DetectRepeats' results")
);

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
in_file <-  opt$in_file
seqID <- gsub('([^\\.]+)\\..*','\\1',in_file)
proc_n <- as.numeric(opt$proc_n)
minScore <- as.numeric(opt$minScore)
cat(seqID, minScore, opt$batch_i, opt$record_len,'\n')

# load libraries
suppressMessages(library(DECIPHER, quietly = TRUE, verbose=FALSE))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
if (!is.na(opt$batch_i)) {
	batch_i <- as.numeric(opt$batch_i)
	record_len <- as.numeric(opt$record_len)
	out_filename <- paste(seqID, '_', batch_i, '.tr', sep='') # set out_filename
	cat('\n', 'read in at most', record_len,'seqs.\n')
	start_i <- (batch_i-1)*record_len + 1
	seq_h <- readAAStringSet(in_file, nrec=record_len, skip=start_i-1)
	seq_h <- RemoveGaps(seq_h, "all")
	if (length(seq_h)!=0) { # if there is seq
		cat('\nRunning DetectRepeats on sequences No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
		result_h <- DetectRepeats(seq_h, minScore=minScore, processors=proc_n)
		result_h$Index <- result_h$Index + start_i-1
		# write output
		saveRDS(result_h, file=out_filename, compress = TRUE)
		cat('\nSave DetectRepeats result as', out_filename, '\n')
		cat('###Result:', seqID, batch_i, length(seq_h), length(unique(result_h$Index)), '\n', sep='\t')
	} else {
		cat('\n!!!No sequence in batch', batch_i, '(No. ', start_i, 'to No.', start_i+length(seq_h)-1, ')...\n\n')
	}
} else {
	out_filename <- paste(seqID, '.tr', sep='') # set out_filename
	seq_h <- readAAStringSet(in_file)
	seq_h <- RemoveGaps(seq_h, "all")
	result_h <- DetectRepeats(seq_h, minScore=minScore, processors=proc_n)
	# write output
	saveRDS(result_h, file=out_filename, compress = TRUE)
	cat('\nSave DetectRepeats result as', out_filename, '\n')
	cat('###Result:', seqID, length(seq_h), length(unique(result_h$Index)), '\n', sep='\t')
}