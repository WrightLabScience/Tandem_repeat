#!/usr/bin/Rscript

# parse_trust.5.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read TRUST's result and output tsv and fastally

# v1 2022/04/21
# v2 2022/03/02
# v3 2022/05/17 Use argparse
# v5 2025/07/07

# Usage:
# Rscript parse_trust.5.r --in_file=/Users/shc167/Documents/project/TR02/TR02.102/CESymm_out_TRUST/TPs_v5_rename.out --out_prefix=/Users/shc167/Documents/project/TR02/TR02.102/TP_out_TRUST.tsv

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
out_tsv <- paste0(opt$out_prefix, '_out_TRUST.tsv')
write(paste('seqID', 'type', 'repeat_length', 'start', 'end', 'score_unit', 'score_sum', sep='\t'), file=out_tsv, append=FALSE)

# read line
lines <- readLines(opt$in_file)
lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
lines <- trimws(lines)
in_seqIDs <- c()
out_file_index=0
output_nchar=20000
out_n_repeats=0
output_lines <- ''
for (line in lines) {
	if (startsWith(line, ">") && substr(line, 2, 7)!='Repeat') { # new seq
		seqID <- gsub(">", "", line)
		in_seqIDs <- c(in_seqIDs, seqID)
		df_repeat_candidates <- data.frame(matrix(nrow=0,ncol=3))
		colnames(df_repeat_candidates) <- c('start', 'end', 'score')
		save_alignment <- FALSE
	} else if (startsWith(line, '# For protein ') && endsWith(line, 'included')) { # store candidate units
		df_repeat_candidates[nrow(df_repeat_candidates)+1,] <- c(gsub('.*repeat (\\d+)\\..*', '\\1', line), gsub('.*\\.\\.(\\d+) scored.*', '\\1', line), gsub('.*scored (\\S+) scored.*', '\\1', line))
	} else if (startsWith(line, "REPEAT_TYPE")) { # save the repeat type index
		typeID <- gsub('REPEAT_TYPE ', '', line)
		out_n_repeats = out_n_repeats + 1
		starts <- c()
		ends <- c()
		scores <- c()
	} else if (startsWith(line, "REPEAT_LENGTH")) { # save the repeat length
		repeat_length <- gsub('REPEAT_LENGTH ', '', line)
	} else if (grepl("# Repeat \\d+$", line)) { # save selected repeat units
		selected_start <- gsub('(\\d+) \\d+\\t.*', '\\1', line)
		starts <- c(starts, selected_start)
		selected_length <- gsub('\\d+ (\\d+)\\t.*', '\\1', line)
		selected_end=as.character(as.numeric(selected_start)+as.numeric(selected_length)-1)
		ends <- c(ends, selected_end)
		scores <- c(scores, as.numeric(df_repeat_candidates[(df_repeat_candidates$start==selected_start & df_repeat_candidates$end==selected_end), 'score']))
	} else if (startsWith(line, "# The multiple alignment of repeats")) { # save selected repeat info to tsv
		score_sum <- sum(scores)
		for (unit_i in seq_along(starts)) { # write each unit to .tsv
			write(paste(seqID, typeID, repeat_length, starts[unit_i], ends[unit_i], scores[unit_i], score_sum, sep='\t'), file=out_tsv, append=TRUE)
		}
		if (nchar(output_lines) > 0) {
			output_lines <- paste0(output_lines, '\n')
		}
		output_lines <- paste0(output_lines, '#', seqID, '_', typeID) # alignment name
	} else if (startsWith(line, ">Repeat ")) {
		save_alignment <- TRUE
	} else if (save_alignment==TRUE) {
		output_lines <- paste0(output_lines, '\n', gsub('[Xx.]','-',line, ignore.case = TRUE))
		save_alignment <- FALSE
	} else if (startsWith(line, "# end of the protein ") && nchar(output_lines) > 0) { # write alignment to file
		if (output_nchar + nchar(output_lines) > 18000) {
			out_file_index=out_file_index+1
			output_nchar = 0
			write(output_lines, file=paste0(opt$out_prefix, '.', out_file_index,'.fastally'), append=FALSE)
		} else {
			write(output_lines, file=paste0(opt$out_prefix, '.', out_file_index,'.fastally'), append=TRUE)
		}
		output_nchar = output_nchar + nchar(output_lines)
		output_lines <- ''
	}
}
cat("n_seq_in = ", length(unique(in_seqIDs)), ", out_n_repeats = ", out_n_repeats,'\n', sep='')

