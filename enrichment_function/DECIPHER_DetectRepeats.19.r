# Usage: 
# Rscript DECIPHER_DetectRepeats.19.r GCA_932276165 1 100

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
batch_i <- as.numeric(ARGS[2])
record_len <- as.numeric(ARGS[3])

print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# set filenames
in_filename <- paste(seqID, '.faa.gz', sep='')
out_filename <- paste(seqID, '_', batch_i, '.tr', sep='')

cat('\n', 'read in at most', record_len,'seqs.\n')

start_i <- (batch_i-1)*record_len + 1

seq_h <- readAAStringSet(in_filename, nrec=record_len, skip=start_i-1)
seq_h <- RemoveGaps(seq_h, "all")

if (length(seq_h)!=0) { # if there is seq
    cat('\nRunning DetectRepeats on sequences No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
    result_h <- DetectRepeats(seq_h)
    result_h$Index <- result_h$Index + start_i-1
    # write output
    saveRDS(result_h, file=out_filename, compress = TRUE)
    cat('\nSave DetectRepeats result as', out_filename, '\n')
    cat('###Result:', seqID, batch_i, length(seq_h), length(unique(result_h$Index)), '\n', sep='\t')    
} else {
    cat('\n!!!No sequence in batch', batch_i, '(No. ', start_i, 'to No.', start_i+length(seq_h)-1, ')...\n\n')
}