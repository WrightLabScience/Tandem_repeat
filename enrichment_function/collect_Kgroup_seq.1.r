# Usage: 
# Rscript collect_Kgroup_seq_1.r GCF_904067545

ARGS <- commandArgs(trailingOnly = TRUE)
genomeID_h <- ARGS[1] # in_dataframe

# load libraries
library(stringr)
library(DECIPHER)

######################## main ########################
# read data frame that contains FTP address and genome names
faa_file <- paste(genomeID_h, '.faa.gz', sep='')
IdTaxa_rds <- paste(genomeID_h, '.id', sep='')

cat('\n', genomeID_h, '\n')

# Id for all genes
IdTaxa_h <- readRDS(IdTaxa_rds)
assignment <- sapply(IdTaxa_h, function(x) paste(x$taxon, collapse=";"))
assignment_df <- as.data.frame(str_split_fixed(assignment, ';',5))
assignment_df <- assignment_df[,c(2:5)] # remove root
colnames(assignment_df) <- c('lv1', 'lv2', 'lv3', 'Kgroup')
cat('\nDone read IdTaxa.\n')

# loop through all genes
for (idx_h in seq(nrow(assignment_df))) {
	Kgroup_h <- assignment_df[idx_h,'Kgroup']
	if ((Kgroup_h != "") & (!startsWith(Kgroup_h, "unclassified_"))) { # if gene has assigned Kgroup
		Kgroup_h <- substr(Kgroup_h, 1, 6)
        cat(idx_h, Kgroup_h, Sys.time(), '\n')
		seq_h <- readAAStringSet(faa_file, nrec=1, skip=idx_h-1)
		filename_h <- paste(genomeID_h, '/', Kgroup_h, '.faa.gz', sep='')
		writeXStringSet(seq_h, filename_h, append=TRUE, compress=TRUE, format="fasta")
	}
}
cat('\nDone write faa.\n')
list.files()