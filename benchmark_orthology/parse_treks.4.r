#!/usr/bin/Rscript

# parse_treks.4.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read T-REKS's result and output tsv

# v1 2025/07/08

# Usage:
# Rscript ~/Documents/scripts_TR/parse_treks2Tally.1.r --in_file=/Users/shc167/Documents/project/TR02/TR02.102/CESymm_out_T-REKS/TPs_v5_rename.out --out_dir=/Users/shc167/Documents/project/TR02/TR02.102/Tally_treks/

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--in_file"), type="character", default="NA", help="tandem repeat result (Required)", metavar="filename"),
	make_option(c("-o", "--out_prefix"), type="character", default="NA", help="prefix of output file (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_file=DIRNAME] [--out_prefix=FILENAME]", option_list=option_list))

# title of .tsv
out_tsv <- paste0(opt$out_prefix, '_out_TREKS.tsv')
write(paste('seqID','repeatID','start','end','Psim', sep='\t'), file=out_tsv, append=FALSE)

# read line
lines <- readLines(opt$in_file)
lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
lines <- trimws(lines)
in_seqIDs <- c()
out_file_index=1
output_nchar=0
out_n_repeats=0
output_lines <- NULL
for (line in lines) {
	if (startsWith(line, ">") && substr(line, 2, 7)!='Repeat') { # new seq
		seqID <- gsub(">", "", line)
		in_seqIDs <- c(in_seqIDs, seqID)
		repeatID <- 0
		save_alignment <- FALSE
	} else if (startsWith(line, 'Length:')) { # save the repeat
		repeatID = repeatID + 1
		start <- gsub('.*from  (\\d+) to.*', '\\1', line)
		end <- gsub('.* to (\\d+) - Psim:.*', '\\1', line)
		score <- as.numeric(gsub('.*Psim:(\\S+) region.*', '\\1', line))
		write(paste(seqID, repeatID, start, end, score, sep='\t'), file=out_tsv, append=TRUE)
		save_alignment <- TRUE
		output_lines <- paste0('#', seqID, '_', repeatID)
	} else if (startsWith(line, '**********************')) { # repeat end
		save_alignment <- FALSE
		if (length(output_lines) > 0) { # write alignment to file
			if (output_nchar + nchar(output_lines) > 18000) {
				out_file_index=out_file_index+1
				output_nchar = 0
			}
			write(output_lines, file=paste0(opt$out_prefix, '.', out_file_index,'.fastally'), append=TRUE)
			output_nchar = output_nchar + nchar(output_lines)
			out_n_repeats = out_n_repeats + 1
			output_lines <- NULL
		}
	} else if (save_alignment==TRUE) {
		output_lines <- paste0(output_lines, '\n', gsub('[Xx.]','-',line, ignore.case = TRUE))
	} 
}
cat("n_seq_in = ", length(unique(in_seqIDs)), ", out_n_repeats = ", out_n_repeats,'\n', sep='')