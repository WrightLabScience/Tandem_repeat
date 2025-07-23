#!/usr/bin/Rscript

# parse_HHrepID2Tally.1.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read HHrepID's result and output fastally

# v1 2025/07/07

# Usage:
# Rscript parse_HHrepID.2.r --in_dir=~/project/TR02/TR02.72/HHrepID.out/7PN0A.aa_T9E-0.out --out_file=~/project/TR02/TR02.72/HHrepID/aa_tsv/7PN0A.tsv

library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-i", "--in_dir"), type="character", default="NA", help="predict result dir containing *.out (Required)", metavar="dirname"),
	make_option(c("-o", "--out_dir"), type="character", default="NA", help="prefix of output .fastally file (Required)", metavar="dirname")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dir=DIRNAME] [--out_dir=FILENAME]", option_list=option_list))

# list of input files
in_file_names <- list.files(path=opt$in_dir, pattern=".out", all.files=TRUE, full.names=TRUE)

# function to write alignment to file
write_alignment <- function(Index_h,typeID,alignment_h,output_nchar,out_file_index, out_dir) {
	output_lines <- paste0('#', Index_h, '_', typeID) # title
	for (unit_ID_h in names(alignment_h)) {
		output_lines <- paste0(output_lines, '\n', alignment_h[[unit_ID_h]])
	}
	if (output_nchar + nchar(output_lines) > 18000) {
		out_file_index=out_file_index+1
		output_nchar = 0
	}
	write(output_lines, file=paste0(out_dir, out_file_index,'.fastally'), append=TRUE)
	output_nchar = output_nchar + nchar(output_lines)
	return(c(out_file_index, output_nchar))
}

# loop through result file
in_seqIDs <- c()
out_n_repeats=0
alignment_h <- list()
out_file_index=1
output_nchar=0
for (in_file_name in in_file_names) {
	Index_h <- gsub('.out', '', basename(in_file_name))
	in_seqIDs <- c(in_seqIDs, Index_h)
	lines <- readLines(in_file_name)
	lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
	lines <- trimws(lines)
	for (line in lines) {
		if (startsWith(line, "Results")) {
			if (length(alignment_h) > 0) { # write previous repeat alignment
				out_index_nchar <- write_alignment(Index_h,typeID,alignment_h,output_nchar,out_file_index, opt$out_dir)
				alignment_h <- list()
				out_file_index <- out_index_nchar[1]
				output_nchar <- out_index_nchar[2]
				out_n_repeats=out_n_repeats+1
			}
			typeID <- gsub("Results for repeats type |:", "", line)
		} else if (grepl(" \\+", line)) {
			words_list <- unlist(strsplit(line, " "))
			unit_ID_h <- words_list[1]
			alignment_unit_h <- gsub('[Xx.]','-',words_list[length(words_list)], ignore.case = TRUE)
			alignment_h[[unit_ID_h]] <- paste0(alignment_h[[unit_ID_h]], alignment_unit_h)
		}
	}
	if (length(alignment_h) > 0) { # write previous repeat alignment
		out_index_nchar <- write_alignment(Index_h,typeID,alignment_h,output_nchar,out_file_index, opt$out_dir)
		alignment_h <- list()
		out_file_index <- out_index_nchar[1]
		output_nchar <- out_index_nchar[2]
		out_n_repeats=out_n_repeats+1
	}
}

cat("n_seq_in = ", length(unique(in_seqIDs)), ", out_n_repeats = ", out_n_repeats,'\n', sep='')