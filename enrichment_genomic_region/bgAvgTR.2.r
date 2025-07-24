# Usage: 
# Rscript bgAvgTR.1.r GCA_932276165

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
print(ARGS)

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

########################################## START OF SCRIPT ##########################################
suppressMessages(library(DECIPHER))

############## Process input files ##############
# read genomic fna for contig names and lengths
filename_gfna <- paste(seqID,'.fna.gz',sep='')
contig_lengths <- c()
contig_index <- 1
end <- FALSE
while (end!=TRUE) {
	seq_h <- readDNAStringSet(filename_gfna, nrec=1, skip=contig_index-1)
	if (length(seq_h)>0) { # if there is seq
		contig_lengths <- c(contig_lengths, setNames(width(seq_h), names(seq_h)))
		contig_index <- contig_index + 1 # next contig
	} else {
		end <- TRUE
	}
}
rm(seq_h)
names(contig_lengths) <- gsub('(\\S+).*', '\\1', names(contig_lengths))

# read TR and get TR ranges
TR_all <- readRDS(paste(seqID,'.tr',sep=''))
TR_covs <- c()
for (i in seq_along(contig_lengths)) {
	TR_h <- TR_all[TR_all$Index == i, ]
	if (nrow(TR_h) == 0) {
		TR_cov <- 0
	} else {
		TR_ranges <- IRanges(start = TR_h$Begin, end = TR_h$End)
		TR_cov <- sum(width(reduce(TR_ranges)))
	}
	cat(sprintf("Contig: %s\tLength: %d\tTR bases (non-overlapping): %d\n",
	    names(contig_lengths)[i], contig_lengths[i], TR_cov))
	TR_covs <- c(TR_covs, TR_cov)
}

cat('# chromosome lengths:', contig_lengths, '\n')
cat('# TR position number:', TR_covs, '\n')
cat('# background_average_TR:', seqID, '\t', sum(TR_covs)/sum(contig_lengths), '\n', sep = "")
cat('############################################################\n\n')

[SUPP DATA] Fig S. TR unit PID distributions. (eCDF curves)