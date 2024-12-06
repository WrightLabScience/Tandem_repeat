# Usage: 
# Rscript DECIPHER_MSA.4.r ${Kgroup_h}.fas.gz ${Kgroup_h}.msa

# re-install DECIPHER
LIB <- "./libraries/"
dir.create(LIB)
Sys.setenv(TMPDIR=getwd())
install.packages("DECIPHER_2.29.3_11272023.tar.gz", lib=LIB)
.libPaths(c(LIB, .libPaths()))


ARGS <- commandArgs(trailingOnly = TRUE)
in_file <- ARGS[1]
out_file <- ARGS[2]
ntaa <- ARGS[3]

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# set filenames
if (ntaa == 'nt') {
    seq <- readDNAStringSet(in_file)
    cat("# Input nt,",in_file, ", number of seqs:", length(seq), '\n')
    msa <- AlignTranslation(seq)
} else {
    seq <- readAAStringSet(in_file)
    cat("# Input aa,",in_file, ", number of seqs:", length(seq), '\n')
    msa <- AlignSeqs(seq)
}
saveRDS(msa, file=out_file, compress = TRUE)