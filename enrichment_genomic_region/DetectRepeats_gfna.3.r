# Usage: 
# Rscript DetectRepeats_gfna.3.r --in_file=GCF_964030795.fna.gz

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--in_file", type="character", default="NA", help="input file (Required)", metavar="filename"),
	make_option("--proc_n", type="integer", default=1, help="number of processors", metavar="number"),
	make_option("--contig_index", type="integer", default=NA, help="contig index", metavar="number")
);

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
in_filename <- opt$in_file
seqID <- gsub('([^\\.]+)\\..*','\\1', in_filename)
proc_n <- as.numeric(opt$proc_n)
contig_index <- as.numeric(opt$contig_index)

# load libraries
suppressMessages(library(DECIPHER, quietly = TRUE, verbose=FALSE))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big
########################################## START OF SCRIPT ##########################################
if (!is.na(opt$contig_index)) {
	seq_h <- readDNAStringSet(in_filename, nrec=1, skip=contig_index-1)
	out_filename <- paste(seqID, '_', contig_index, '.tr', sep='')
	if (length(seq_h)>0) { # if there is seq
		cat('\nRunning DetectRepeats on contig No.', contig_index, ":", names(seq_h), width(seq_h), '...\n\n')
		TR_h <- DetectRepeats(seq_h, processors=proc_n)
		if (nrow(TR_h)>0) {
			TR_h$Index <- contig_index
		}
		saveRDS(TR_h, file=out_filename)
		cat('#### DONE ####', nrow(TR_h), 'TR detected in', seqID, 'contig No.', contig_index, '\n')
	} else {
		cat('!!!!!There is no contig No.', contig_index, 'in', seqID, '\n')
	}
} else {
	seq_h <- readDNAStringSet(in_filename)
	out_filename <- paste(seqID, '.tr', sep='')
	TR_h <- DetectRepeats(seq_h, processors=proc_n)
	saveRDS(TR_h, file=out_filename)
	cat('#### DONE ####', nrow(TR_h), 'TR detected in', seqID, '\n')
}