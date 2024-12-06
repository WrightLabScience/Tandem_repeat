#!/usr/bin/Rscript

# load libraries
library(optparse, quietly = TRUE, verbose=FALSE)

# specify options in a list
option_list = list(
	make_option("--path_answer", type="character", default="/Users/shc167/Documents/project/TR02/TR02.98/answer.tsv", help="answer", metavar="filename"),
	make_option("--path_aa", type="character", default="NA", help="aa predict *.tsv (Required)", metavar="filename"),
	make_option("--path_decoy", type="character", default="NA", help="decoy predict *.tsv (Required)", metavar="filename"),
	make_option("--path_out", type="character", default="NA", help="ROC result *.rds (Required)", metavar="filename"),
    make_option("--num_cutoffs", type="integer", default=1000, help="Maximum number of scores for plotting", metavar="number"),
    make_option("--le_score", action = "store_true", default = FALSE, help = "lower than or equal to score cutoff (for P-value or err)"),
    make_option("--default_score", type="double", default=10, help="score cutoff", metavar="number"),
    make_option("--colname_seqID", type="character", default="seqID", help="colname_seqID", metavar="colname"),
    make_option("--colname_start", type="character", default="start", help="colname_start", metavar="colname"),
    make_option("--colname_end", type="character", default="end", help="colname_end", metavar="colname"),
    make_option("--colname_score", type="character", default="score", help="colname_score", metavar="colname")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "usage: %prog [options]", option_list=option_list))

# input
path_answer <- opt$path_answer #'/Users/shc167/Documents/project/TR02/TR02.98/answer.tsv'
path_aa <- opt$path_aa #'/Users/shc167/Documents/project/TR02/TR02.98/DetectRepeats_wE_aa_05102024.tsv'
path_decoy <- opt$path_decoy #'/Users/shc167/Documents/project/TR02/TR02.98/DetectRepeats_wE_decoy_05102024.tsv'
path_out <- opt$path_out #'/Users/shc167/Documents/project/TR02/TR02.98/ROC_DetectRepeats_wE.2.rds'
num_cutoffs <- opt$num_cutoffs #1000
#le_score <- opt$le_score #FALSE
default_score <- opt$default_score #10
colname_seqID <- opt$colname_seqID #'seqID'
colname_start <- opt$colname_start #'Begin'
colname_end <- opt$colname_end #'End'
colname_score <- opt$colname_score #'Score'

library(DECIPHER, quietly = TRUE, verbose=FALSE)


################# functions #################
# vector of TR positions
get_vector_TR <- function(seqID_h, seq_len, TR, colname_seqID, colname_start, colname_end) {
    vector_TR <- numeric(seq_len[seqID_h])
    for (row_i in rownames(TR[TR[,colname_seqID]==seqID_h,])) {
        vector_TR[TR[row_i,colname_start]:TR[row_i,colname_end]] <- 1
    }
    return(vector_TR)
}

################# main #################
# read input
answer <- read.table(path_answer, header = TRUE, sep= "\t")
predict_aa <- read.table(path_aa, header = TRUE, sep= "\t")
predict_decoy <- read.table(path_decoy, header = TRUE, sep= "\t")

# seq_len
seq_len <- answer[,c('ID1', 'seq_len')]
seq_len <- seq_len[!duplicated(seq_len),]
seq_len <- setNames(seq_len$seq_len, seq_len$ID1)

# get list of score cutoff
scores <- unique(sort(c(predict_aa[,colname_score], predict_decoy[,colname_score])))
if (length(scores) > num_cutoffs) {
    cutoffs <- quantile(scores, seq(0, 1, length.out=num_cutoffs+1))
} else {
    cutoffs <- scores
}
if (!is.element(default_score, cutoffs)) {
    cutoffs <- c(cutoffs, default_score)
    cutoffs <- sort(cutoffs)
}

cat('Number of cutoffs = ', length(cutoffs), '\n', sep='')

# TR answer
list_answer <- list()
for (seqID_h in names(seq_len)) {
    list_answer[[seqID_h]] <- get_vector_TR(seqID_h, seq_len, answer, 'ID1', 'begin_label', 'end_label')
}


# calculate TP and FP for each cutoff
ROC_df <- as.data.frame(matrix(, nrow = length(cutoffs), ncol = 3))
colnames(ROC_df) <- c('cutoff', 'TPR', 'FPR')

i=1
for (cutoff_h in cutoffs) {
    if (opt$le_score) { # filter TR result by score cutoff
        predict_aa_filtered <- predict_aa[predict_aa[,colname_score]<=cutoff_h,]
        predict_decoy_filtered <- predict_decoy[predict_decoy[,colname_score]<=cutoff_h,]
    } else {
        predict_aa_filtered <- predict_aa[predict_aa[,colname_score]>=cutoff_h,]
        predict_decoy_filtered <- predict_decoy[predict_decoy[,colname_score]>=cutoff_h,]
    }
    # vector of TPR or FPR for each seq
    TPRs <- c()
    FPRs <- c()
    for (seqID_h in names(seq_len)) {
        # TPR = overlapped TR / full length of answer TR
        vector_predict_aa_h <- get_vector_TR(seqID_h, seq_len, predict_aa_filtered, colname_seqID, colname_start, colname_end)
        TPRs <- c(TPRs, as.numeric(list_answer[[seqID_h]] %*% vector_predict_aa_h) / sum(list_answer[[seqID_h]]))
        # FPR = TR / full length of chain
        vector_predict_decoy_h <- get_vector_TR(seqID_h, seq_len, predict_decoy_filtered, colname_seqID, colname_start, colname_end)
        FPRs <- c(FPRs, sum(vector_predict_decoy_h)/seq_len[seqID_h])
    }
    TPR_mean <- mean(TPRs)
    FPR_mean <- mean(FPRs)
    cat(i, cutoff_h, TPR_mean, FPR_mean, '\n',sep='\t')
    ROC_df[i,] <- c(cutoff_h, TPR_mean, FPR_mean)
    i=i+1
}

saveRDS(ROC_df, path_out, compress = TRUE)