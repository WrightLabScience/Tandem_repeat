# Usage: 
# Rscript TR_dist.15.r K00001

# read inputs
ARGS <- commandArgs(trailingOnly = TRUE)
Kgroup_h <- ARGS[1]

print(ARGS)

MSA_h <- readRDS(paste(Kgroup_h, '.msa', sep=''))
TR_h <- readRDS(paste(Kgroup_h, '.tr', sep=''))

start_time <- Sys.time()

# load libraries
suppressMessages(library(DECIPHER))
packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## functions ########################
# matrix of MSA
get_mt_MSA <- function(MSA_h) {
	mt_MSA <- matrix(, nrow = width(MSA_h)[1], ncol = length(MSA_h))
	for (index_h in seq_along(MSA_h)) {
		mt_MSA[,index_h] <- strsplit(toString(MSA_h[index_h]),"")[[1]]
	}
	return(mt_MSA)
}

# matrix of repeat loci in MSA
get_mt_rep_loci <- function(TR_h, mt_MSA) {
	mt_rep_loci <- mt_MSA
	mt_rep_loci[mt_rep_loci!='-'] <- 0
	for (row_h in seq(nrow(TR_h))) {
		index_h <- TR_h[row_h, 'Index']
		begin_h <- TR_h[row_h, 'Begin']
		end_h <- TR_h[row_h, 'End']
		# convert location in ori_seq to location in alignment(with "-")
		# marked as:
		# 1: if location is in repeat region
		# 0: if location is not in repeat region
		# 1: if is gap inside TR
		# 0: if is gap outside TR
		pos_ori <- 0
		for (pos_align in seq(nrow(mt_MSA))) {
			if (mt_rep_loci[pos_align,index_h]!='-') {
				pos_ori <- pos_ori + 1
				if (pos_ori >= begin_h & pos_ori <= end_h ) {
					mt_rep_loci[pos_align, index_h] <- 1
				}
			} else if (pos_ori >= begin_h-1 & pos_ori <= end_h+1 ) { # 1: if is gap inside TR
				mt_rep_loci[pos_align, index_h] <- 1
			}
		}
	}
	mt_rep_loci[mt_rep_loci=='-'] <- 0 # 0: if is gap outside TR
	mt_rep_loci <- matrix(as.numeric(mt_rep_loci),ncol = ncol(mt_rep_loci))
	return(mt_rep_loci)
}

# get TR region for 1 seq 1 TR
get_vec_rep_loci <- function(TR_e_h, MSA_i_h) {
	vec_rep_loci <- MSA_i_h
	vec_rep_loci[vec_rep_loci!='-'] <- 0
	begin_h <- TR_e_h$Begin
	end_h <- TR_e_h$End
	# convert location in ori_seq to location in alignment(with "-")
	# marked as:
	# 1: if location is in repeat region
	# 0: if location is not in repeat region
	# 1: if is gap inside TR
	# 0: if is gap outside TR
	pos_ori <- 0
	for (pos_align in seq_along(MSA_i_h)) {
		if (vec_rep_loci[pos_align]!='-') {
			pos_ori <- pos_ori + 1
			if (pos_ori >= begin_h & pos_ori <= end_h ) {
				vec_rep_loci[pos_align] <- 1
			}
		} else if (pos_ori >= begin_h-1 & pos_ori <= end_h+1 ) { # 1: if is gap inside TR
			vec_rep_loci[pos_align] <- 1
		}
	}
	vec_rep_loci[vec_rep_loci=='-'] <- 0 # 0: if is gap outside TR
	vec_rep_loci <- as.numeric(vec_rep_loci)
	return(vec_rep_loci)
}

# get unit length
get_unit_len <- function(TR_h, output){
	tryCatch(
		{
		# accessing elements from third column
		units_h <- IRanges(TR_h$Left, TR_h$Right)
		#unit_len_h <- mean(width(units_h))
		unit_len_h <- median(width(units_h))
		return(unit_len_h)
		},
		#if an error occurs, tell me the error
		error=function(e) {
			message('An Error Occurred')
			print(e)
			print(TR_h)
		},
		#if a warning occurs, tell me the warning
		warning=function(w) {
			message('A Warning Occurred')
			print(w)
			return(NA)
		}
	)
}


######################## main ########################

# vector to store distances
dists <- c()
copy_numbers <- c()
count_all_TR <- 0

if (!is.null(TR_h) & nrow(TR_h)>0) { # ignore K group if no repeat detect
    #cat('\n', Kgroup_h, 'has TR.\n')
    count_seq <- length(MSA_h)
    count_all_TR <- length(unique(TR_h$Index))
    # filter TR by unit length median >= 10
    unit_len_h <- apply(TR_h,1,get_unit_len)
    TR_h <- TR_h[unit_len_h >=10,]
    count_10 <- length(unique(TR_h$Index))
    if (count_10 > 1) {
        # MSA
        mt_MSA <- get_mt_MSA(MSA_h)
        # matrix of repeat loci in MSA
        mt_rep_loci <- get_mt_rep_loci(TR_h, mt_MSA)
        # get main TR region
        row_sums_h <- rowSums(mt_rep_loci)
        max_val <- max(row_sums_h)
        if (max_val > count_seq/2) { # there should be more than 1/2 seqs have the main TR
            max_i <- which.max(row_sums_h)
            threshold_h <- max_val/2
            # get start pos
            pos_i <- max_i
            while ((row_sums_h[pos_i] >= threshold_h) & (pos_i >= 1) & (pos_i <= length(row_sums_h))) {
                TR_start_i <- pos_i
                pos_i <- pos_i-1
                if (pos_i < 1) {
                    break
                }
            }
            # get end pos
            pos_i <- max_i
            while ((row_sums_h[pos_i] >= threshold_h) & (pos_i >= 1) & (pos_i <= length(row_sums_h))) {
                TR_end_i <- pos_i
                pos_i <- pos_i+1
                if (pos_i < 1) {
                    break
                }
            }
            #cat('\n Main TR region', TR_start_i, '-', TR_end_i, '\n')
            # make a vector for main TR region
            TR_region_main <- integer(length(row_sums_h))
            TR_region_main[c(TR_start_i:TR_end_i)] <- 1

            # check for each seq and each TR, get mode of unit len
            count_overlap <- 0
            # unit length
            unit_len_l <- c()
            for (seq_index_h in seq_along(MSA_h)) {
                TR_region_h <- mt_rep_loci[,seq_index_h]
                if (sum(TR_region_h * TR_region_main) > 0) { # has main TR
                    TR_i_h <- TR_h[TR_h$Index==seq_index_h,]
                    for (row_h in seq(nrow(TR_i_h))) {
                        TR_e_h <- TR_i_h[row_h,] # get one TR element
                        MSA_i_h <- mt_MSA[,seq_index_h]
                        vec_rep_loci <- get_vec_rep_loci(TR_e_h, MSA_i_h)
                        if (sum(vec_rep_loci * TR_region_main) > 0) { # this TR element is main TR
                            seq_h <- RemoveGaps(MSA_h[seq_index_h], "all")
                            # get dist
                            reps <- extractAt(seq_h[[1]], IRanges(TR_e_h[[1, "Left"]], TR_e_h[[1, "Right"]]))
                            unit_len_l <- c(unit_len_l, width(reps))
                            count_overlap <- count_overlap + 1
                        }
                    }
                }
            }
            mode_unit_len <- as.numeric(names(sort(-table(unit_len_l))[1]))

            # only look at conserved TR
            count_cons <- 0
            for (seq_index_h in seq_along(MSA_h)) {
                TR_region_h <- mt_rep_loci[,seq_index_h]
                if (sum(TR_region_h * TR_region_main) > 0) { # has main TR
                    TR_i_h <- TR_h[TR_h$Index==seq_index_h,]
                    for (row_h in seq(nrow(TR_i_h))) {
                        TR_e_h <- TR_i_h[row_h,] # get one TR element
                        MSA_i_h <- mt_MSA[,seq_index_h]
                        vec_rep_loci <- get_vec_rep_loci(TR_e_h, MSA_i_h)
                        if (sum(vec_rep_loci * TR_region_main) > 0) { # this TR element is main TR
                            seq_h <- RemoveGaps(MSA_h[seq_index_h], "all")
                            # get dist
                            reps <- extractAt(seq_h[[1]], IRanges(TR_e_h[[1, "Left"]], TR_e_h[[1, "Right"]]))
                            reps <- AlignSeqs(reps, verbose=FALSE) # align the repeats
                            #reps <- MaskAlignment(reps) # mask alignment
                            reps <- as(MaskAlignment(reps, includeTerminalGaps=TRUE), "AAStringSet")
                            mask_aln_len <- width(reps)[1]
                            if ((mask_aln_len <= mode_unit_len*1.1) & (mask_aln_len >= mode_unit_len*0.9)) {
                                distances_m <- DistanceMatrix(reps, includeTerminalGaps = TRUE, penalizeGapLetterMatches = NA, verbose=FALSE)
                                distances_h <- distances_m[lower.tri(distances_m)] # distance values
                                count_cons <- count_cons + 1
                                dists <- c(dists, distances_h)
                                copy_numbers <- c(copy_numbers, length(reps))
                            }
                        }
                    }
                }
            }
            #cat('\nDONE calculating distances on', count_overlap, 'TR elements (len>=10 main TR),', count_seq, 'genes in total...\n\n')
            if (length(dists) > 0) {
                mode_copy_number <- as.numeric(names(sort(-table(copy_numbers))[1]))
                cat('### Report counts include:', Kgroup_h, count_seq, count_all_TR, count_10, count_overlap, count_cons, mode_unit_len, mode_copy_number, '\n', sep = '\t')
                saveRDS(dists, file=paste(Kgroup_h, '.dist', sep=''), compress = TRUE)
                saveRDS(copy_numbers, file=paste(Kgroup_h, '.copyn', sep=''), compress = TRUE)
            } else {
                cat('Report counts exclude:', Kgroup_h, count_seq, count_all_TR, count_10, count_overlap, count_cons, mode_unit_len, NA, '\n', sep = '\t')
            }
        } else {
            cat('Report counts exclude:', Kgroup_h, count_seq, count_all_TR, count_10, 0, NA, NA, NA, '\n', sep = '\t')
        }
    } else {
        cat('\nTR (unit len >10 ) < 2 in', Kgroup_h, '\n')
    }
}
if (count_all_TR == 0) {
	cat('\nTR not detected in', Kgroup_h, '\n')
}
end_time <- Sys.time()
cat('Done TR_dist', end_time - start_time, 'sec.\n')