# Usage: 
# Rscript DECIPHER_DetectRepeats.18.r GCA_932276165

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# set filenames
in_filename <- paste(seqID, '.faa.gz', sep='')
out_filename <- paste(seqID, '.tr', sep='')

record_len <- 100
cat('\n', 'read in at most', record_len,'seqs.\n')

result_all <- NULL # create space holder for rbind
end <- FALSE
start_i <- 1
while (end!=TRUE) {
    # read in at most record_len seqs, starts from start_i
    seq_h <- readAAStringSet(in_filename, nrec=record_len, skip=start_i-1)
    seq_h <- RemoveGaps(seq_h, "all")
    if (length(seq_h)!=0) { # if there is seq
        cat('\nRunning DetectRepeats on CDS No.', start_i, 'to No.', start_i+length(seq_h)-1, '...\n\n')
        result_h <- DetectRepeats(seq_h)
        result_h$Index <- result_h$Index + start_i-1
        result_all <- rbind(result_all, result_h)
        start_i <- start_i + length(seq_h) # next round should starts from
    } else {
        end <- TRUE
    }
}
cat('\nDONE running DetectRepeats on', (start_i-1), 'CDS...\n\n')
# write output
saveRDS(result_all, file=out_filename, compress = TRUE)
cat('\nSave DetectRepeats result as', out_filename, '\n')

frac_seqs_w_TR <- length(unique(result_all$Index))/(start_i-1)
cat('###Result:', seqID, frac_seqs_w_TR, '\n', sep='\t')