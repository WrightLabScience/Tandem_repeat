#!/usr/bin/Rscript

# parse_HHrepID.2.r

# Shu-Ting Cho <vivianlily6@hotmail.com>
# Read HHrepID's result and output tsv

# v2 2024/02/25
# v1 2023/06/06

# Usage:
# Rscript parse_HHrepID.2.r --in_dir=~/project/TR02/TR02.72/HHrepID.out/7PN0A.aa_T9E-0.out --out_file=~/project/TR02/TR02.72/HHrepID/aa_tsv/7PN0A.tsv

library(optparse)

# specify options in a list
option_list = list(
    make_option(c("-i", "--in_dir"), type="character", default="NA", help="predict result dir containing *.out (Required)", metavar="dirname"),
    make_option(c("-o", "--out_file"), type="character", default="NA", help="output .tsv file (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--in_dir=DIRNAME] [--out_file=FILENAME]", option_list=option_list))

# write heading to output file
write(paste('seqID', 'Index', 'start','end','P-value','Type','Repeats','Length','Offset', sep='\t'), file=opt$out_file, append=FALSE)

# list of input files
in_file_names <- list.files(path=opt$in_dir, pattern=".out", all.files=TRUE, full.names=TRUE)

# opening result file
in_seqIDs <- c()
out_seqIDs <- c()
for (in_file_name in in_file_names) {
    Index_h <- gsub('.out', '', basename(in_file_name))
    in_seqIDs <- c(in_seqIDs, Index_h)
    lines <- readLines(in_file_name)
    lines <- iconv(lines, "UTF-8", "UTF-8", sub = "byte")
    lines <- trimws(lines)
    for (line in lines) {
        if (startsWith(line, "Results")) {
            typeID <- gsub("Results for repeats type |:", "", line)
        } else if (startsWith(line, "Repeats")) {
            repeats <- gsub("Repeats| ", "", line)
        } else if (startsWith(line, "P-value")) {
            Pvalue <-        gsub("P-value| ", "", line)
        } else if (startsWith(line, "Length")) {
            Length <- gsub("Length| ", "", line)
        } else if (startsWith(line, "Offset")) {
            Offset <- gsub("Offset| ", "", line)
        } else if (grepl(" \\+", line)) {
            words_list <- unlist(strsplit(line, " +"))
            seqID <- words_list[2]
            positions <- unlist(strsplit(words_list[3], "-"))
            start <- positions[1]
            end <- positions[2]
            write(paste(seqID, Index_h, start, end, Pvalue, typeID, repeats, Length, Offset, sep="\t"), file=opt$out_file, append=TRUE)
            out_seqIDs <- c(out_seqIDs, Index_h)
        }
    }
}

cat("n_seq_in = ", length(unique(in_seqIDs)), ", n_seq_out = ", length(unique(out_seqIDs)),'\n', sep='')