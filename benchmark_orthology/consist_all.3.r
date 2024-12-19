#!/usr/bin/env Rscript

# Usage: 
# Rscript consist_all.3.r --Kgroup_list=Kgroup_ge500.rds --mt_MSA_dir=mt_MSA_aa/ --pairs_dir=pairs/ --DetectRepeat_dir=DetectRepeats_wE_aa/ --tr_ext='.wEtr' --pair_ext='_pairs.rds' --msa_ext='_mt_MSA_aa.rds' --score=${score} --out_tsv=${score}.tsv



# load libraries
library(optparse)

# specify options in a list
option_list = list(
	make_option(c("-k", "--Kgroup_list"), type="character", default="NA", help="Kgroup_list.rds (Required)", metavar="filename"),
	make_option(c("-m", "--mt_MSA_dir"), type="character", default="NA", help="dir for *_mt_MSA_*.rds (Required)", metavar="filename"),
	make_option(c("-p", "--pairs_dir"), type="character", default="NA", help="dir for *_pairs.rds (Required)", metavar="filename"),
	make_option(c("-r", "--DetectRepeat_dir"), type="character", default="NA", help="Directory path for *.tr (Required)", metavar="dirname"),
    make_option("--notDECIPHER", action = "store_true", default = FALSE, help = "input tr are not from DECIPHER"),
    make_option("--fromRADAR", action = "store_true", default = FALSE, help = "input tr are from RADAR"),
	make_option("--fromTRUST", action = "store_true", default = FALSE, help = "input tr are from TRUST"),
    make_option("--in_fasta_dir", type="character", default="NA", help="to get seq Index", metavar="dirname"),
    make_option("--fasta_ext", type="character", default=".faa.gz", help="extension of fasta filename", metavar="ext"),
    make_option("--name_seqID", type="character", default="Index", help="name_seqID", metavar="filename"),
    make_option("--name_start", type="character", default="Begin", help="name_start", metavar="filename"),
    make_option("--name_end", type="character", default="End", help="name_end", metavar="filename"),
    make_option("--name_score", type="character", default="Score", help="name_score", metavar="filename"),
    make_option("--le_score", action = "store_true", default = FALSE, help = "lower than or equal to score cutoff (for P-value or err)"),
	make_option("--tr_ext", type="character", default=".tr", help="extension of TR filename", metavar="ext"),
	make_option("--pair_ext", type="character", default="_pairs.rds", help="extension of pair filename", metavar="ext"),
	make_option("--msa_ext", type="character", default="_mt_MSA_aa.rds", help="extension of mt_msa filename", metavar="ext"),
	make_option(c("-s", "--score"), type="double", default=10, help="score cutoff (Required)", metavar="number"),
	make_option(c("-o", "--out_tsv"), type="character", default="NA", help="output tsv (Required)", metavar="filename")
); 

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
opt <- parse_args(OptionParser(usage = "%prog [options] [--Kgroup_list=FILENAME] [--DetectRepeat_dir=DIRNAME] [--score=FLOAT] [--out_tsv=FILENAME]", option_list=option_list))

Kgroup_list <- readRDS(opt$Kgroup_list)
mt_MSA_dir <- opt$mt_MSA_dir
pairs_dir <- opt$pairs_dir
DetectRepeat_dir <- opt$DetectRepeat_dir

tr_ext <- opt$tr_ext
pair_ext <- opt$pair_ext
msa_ext <- opt$msa_ext
score <- opt$score
out_tsv <- opt$out_tsv

######################## functions ########################
# matrix of repeat loci in MSA
get_mt_TR <- function(DRresult_h, mt_MSA, MSA_width) {
	mt_TR <- mt_MSA
	DRresult_h <- DRresult_h[order(DRresult_h$Index), ]
	Index_cur <- 'none'
	for (row_h in seq(nrow(DRresult_h))) {
		index_h <- DRresult_h[row_h, 'Index']
		begin_h <- DRresult_h[row_h, 'Begin']
		end_h <- DRresult_h[row_h, 'End']
		if (index_h!=Index_cur) {
			if (Index_cur!='none') {
				# new seq, store old seq
				mt_TR[!is.na(mt_MSA[,Index_cur]),Index_cur] <- seq_v
			}
			# get new seq vector
			Index_cur <- index_h
			seq_v <- mt_MSA[,Index_cur][!is.na(mt_MSA[,Index_cur])]
		}
        seq_v[c(begin_h:end_h)] <- 1
	}
	mt_TR[!is.na(mt_MSA[,Index_cur]),Index_cur] <- seq_v
	return(mt_TR)
}

# pairwise consistency for positive locations ("1"), mask gaps (NA)
get_pair_consistency <- function(pairs, mt_TR) {
	pair_consistency <- c()
	for (pair_h in seq(ncol(pairs))) {
		pair_mt <- mt_TR[,pairs[,pair_h]] # 2 cols
		pair_mt <- pair_mt[complete.cases(pair_mt), ] # mask gaps
		pair_consistency_h <- (sum(pair_mt[,1] == pair_mt[,2]))/nrow(pair_mt) # count same / all non-gap positions
		pair_consistency <- c(pair_consistency, pair_consistency_h)
    }
	return(pair_consistency)
}

# fraction of seq is in repeat region
get_frac_repeat_h <- function(mt_TR) {
	n_repeat_loci <- colSums(mt_TR==1,na.rm=TRUE)
	n_not_repeat_loci <- colSums(mt_TR==0,na.rm=TRUE)
	return(n_repeat_loci / (n_repeat_loci+n_not_repeat_loci))
}

# get weighted mean
func_w_mean <- function(d, i) {
	d2 <- d[i,]
	return(weighted.mean(d2$mean_cons, d2$n_seqs_inpairs))
}

# convert seqID to Index
format2DECIPHER <- function(DRresult_h, in_fasta, name_seqID, name_start, name_end, name_score) {
    if (opt$fromRADAR) {
        DRresult_h$Index <- DRresult_h[[name_seqID]]
    } else {
		seq_names <- gsub('(\\S+)\\|.*', '\\1', names(readAAStringSet(in_fasta)))
        if (!opt$fromTRUST) {
			DRresult_h[[name_seqID]] <- gsub('(\\S+)\\|.*', '\\1', DRresult_h[[name_seqID]])
		}
        DRresult_h$Index <- match(DRresult_h[[name_seqID]], seq_names)
    }
    DRresult_h[,c('Begin', 'End', 'Score')] <- DRresult_h[,c(name_start, name_end, name_score)]
    return(DRresult_h)
}

######################## main ########################
suppressMessages(library(DECIPHER))
suppressMessages(library(dplyr))
#suppressMessages(library(stringr))
suppressMessages(library(boot))

# Y: Consistency:
# X: mean(fraction of seq is in repeat region) %

line <- paste('KgroupID', 'n_seqs', 'n_seqs_r', 'n_pairs', 'n_seqs_inpairs', 'mean_frac_repeat_all', 'mean_cons', sep='\t')
write(line, file=out_tsv, append=FALSE)

df_h <- data.frame(matrix(ncol = 6, nrow = length(Kgroup_list)))
rownames(df_h) <- Kgroup_list
colnames(df_h) <- c('n_seqs', 'n_seqs_r', 'n_pairs', 'n_seqs_inpairs', 'mean_frac_repeat_all', 'mean_cons')

for (Kgroup_i in seq_along(Kgroup_list)){
	Kgroup_h <- Kgroup_list[Kgroup_i]
    # get the filtered results
    TR_filename <- paste(DetectRepeat_dir, Kgroup_h, tr_ext, sep='')
    if (opt$notDECIPHER) {
        in_fasta <- paste(opt$in_fasta_dir, Kgroup_h, opt$fasta_ext, sep='')
        DRresult_h <- read.delim(TR_filename, header=TRUE, sep = "\t", quote="")
        if (nrow(DRresult_h)>0) {
            DRresult_h <- format2DECIPHER(DRresult_h, in_fasta, opt$name_seqID, opt$name_start, opt$name_end, opt$name_score)
            DRresult_h <- DRresult_h[!is.na(DRresult_h$Score),]
            if (opt$le_score) {
                DRresult_h <- DRresult_h[DRresult_h$Score <= score,]
            } else {
                DRresult_h <- DRresult_h[DRresult_h$Score >= score,]
            }
        }
    } else {
        DRresult_h <- readRDS(TR_filename)
        DRresult_h <- DRresult_h[DRresult_h$Score >= score,]
    }
	n_seqs_r <- n_distinct(DRresult_h$Index) # seq have repeat
    # read MSA # (MSA width) rows x (n seqs) cols
	mt_MSA <- readRDS(paste(mt_MSA_dir, Kgroup_h, msa_ext, sep=''))
    n_seqs <- ncol(mt_MSA) # number of seq in Kgroup
    MSA_width <- nrow(mt_MSA)
	# read pairs # 2 rows x (n pairs) cols 
	mt_pairs <- readRDS(paste(pairs_dir, Kgroup_h, pair_ext, sep=''))
	n_pairs <- ncol(mt_pairs)
	n_seqs_inpairs <- n_distinct(c(mt_pairs)) # flatten then count unique
	if (nrow(DRresult_h)==0) {
		# if no repeat detect, add 0 to frac_repeat, consistency = 100%
		cat(Kgroup_i, Kgroup_h, n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, 0, 1, '\n', sep='\t')
		line <- paste(Kgroup_h, n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, 0, 1, sep='\t')
		write(line, file=out_tsv, append=TRUE)
        df_h[Kgroup_h,] <- c(n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, 0, 1)
	} else {
		mt_TR <- get_mt_TR(DRresult_h, mt_MSA, MSA_width)
		# pairwise consistency for all locations: 1 & 0, mask gaps (NA)
		mean_cons <- mean(get_pair_consistency(mt_pairs, mt_TR))
		mean_frac_repeat_all <- mean(get_frac_repeat_h(mt_TR)) # fraction of seq is in repeat region
		# report result
		cat(Kgroup_i, Kgroup_h, n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, mean_frac_repeat_all, mean_cons, '\n', sep='\t')
		line <- paste(Kgroup_h, n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, mean_frac_repeat_all, mean_cons, sep='\t')
		write(line, file=out_tsv, append=TRUE)
		df_h[Kgroup_h,] <- c(n_seqs, n_seqs_r, n_pairs, n_seqs_inpairs, mean_frac_repeat_all, mean_cons)
	}
}


# Weight the mean of a Kgroup base on # of seqs be considered

cat('### Calculating values for plotting...\n')
# calculate values for plotting
frac_repeat_all_h <- weighted.mean(df_h$mean_frac_repeat_all, df_h$n_seqs) # x
w_mean_h <- weighted.mean(df_h$mean_cons, df_h$n_seqs_inpairs) # y at solid line
if (n_distinct(df_h$mean_cons)>1) {
    # calculate ci
    bootmean <- boot(df_h, func_w_mean, R=1000)
    ci_h <- boot.ci(boot.out = bootmean, conf = 0.95, type = 'perc')
    # store values
    ci_lb <- ci_h[["percent"]][4]
    ci_ub <- ci_h[["percent"]][5]
} else {
    # store values
    ci_lb <- w_mean_h
    ci_ub <- w_mean_h
}
cat('###Result:', score, frac_repeat_all_h, ci_lb, w_mean_h, ci_ub, '\n', sep='\t')