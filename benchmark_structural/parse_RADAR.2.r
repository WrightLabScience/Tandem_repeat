#!/usr/bin/Rscript

# parse_RADAR.2.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read RADAR's result and output tsv

# v2 2025/07/09
# v1 2023/11/16

# Usage:
# Rscript parse_RADAR.2.r --in_dir=/net/dali/home/mscbio/shc167/project/TR02/TR02.102/TP_out_RADAR/ --out_prefix=/net/dali/home/mscbio/shc167/project/TR02/TR02.102/CESymm_out_RADAR/TP

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--in_dir"), type="character", default="NA", help="tandem repeat result (Required)", metavar="filename"),
	make_option(c("-o", "--out_prefix"), type="character", default="NA", help="prefix of output file (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dir=DIRNAME] [--out_file=FILENAME]", option_list=option_list))

# write heading to output file
out_tsv <- paste0(opt$out_prefix, '_out_RADAR.tsv')
write(paste('seqID','start','end','score','No_of_Repeats','Total_Score','Length','Diagonal','BW_From','BW_To','Level', sep='\t'), file=out_tsv, append=FALSE)

# list of input files
in_file_names <- list.files(path=opt$in_dir, pattern=".out.txt", all.files=TRUE, full.names=TRUE)

# opening result file
in_seqIDs <- c()
out_file_index=0
output_nchar=20000
out_n_repeats=0
output_lines <- NULL
repeat_unit <- NULL
for (in_file_name in in_file_names) {
	seqID <- gsub('.out.txt', '', basename(in_file_name))
	in_seqIDs <- c(in_seqIDs, seqID)
	lines <- readLines(in_file_name)
	lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
	lines <- trimws(lines)
	save_repeat <- FALSE
	save_alignment <- FALSE
	for (line in lines) {
		if (startsWith(line, "No. of Repeats")) {
			save_repeat <- TRUE
			starts <- c()
			ends <- c()
		} else if (save_repeat==TRUE) {
			repeat_info <- unlist(strsplit(line,'\\|\\s+'))
			output_lines <- paste0('#', seqID, '_', repeat_info[length(repeat_info)]) # alignment name
			repeat_info <- paste(repeat_info, collapse ='\t')
			save_repeat <- FALSE
			save_alignment <- TRUE
			repeat_unit <- NULL
		} else if (save_alignment==TRUE & !startsWith(line, "-")) { # save alignment
			start <- gsub('\\s*(\\d+)-.*', '\\1', line)
			end <- gsub('\\s*\\d+-\\s*(\\d+) \\(.*', '\\1', line)
			write(paste(seqID, start, end, repeat_info, sep='\t'), file=out_tsv, append=TRUE)
			repeat_unit <- gsub('.*\\)\\s+(\\S+)$', '\\1', line)
			output_lines <- paste0(output_lines, '\n', gsub('[Xx.]','-',repeat_unit, ignore.case = TRUE))
		} else if (save_alignment==TRUE & !is.null(repeat_unit) & startsWith(line, "-")) { # write alignment to file
			if (output_nchar + nchar(output_lines) > 18000) {
				out_file_index=out_file_index+1
				output_nchar = 0
				write(output_lines, file=paste0(opt$out_prefix, '.', out_file_index,'.fastally'), append=FALSE)
			} else {
				write(output_lines, file=paste0(opt$out_prefix, '.', out_file_index,'.fastally'), append=TRUE)
			}
			output_nchar = output_nchar + nchar(output_lines)
			out_n_repeats = out_n_repeats + 1
			output_lines <- NULL
			save_alignment <- FALSE
		}
	}
}

cat("n_seq_in = ", length(unique(in_seqIDs)), ", out_n_repeats = ", out_n_repeats,'\n', sep='')