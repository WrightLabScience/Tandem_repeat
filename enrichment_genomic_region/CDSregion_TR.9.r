# Usage: 
# Rscript CDSregion_TR.9.r GCA_932276165

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]

print(ARGS)

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

########################################## START OF SCRIPT ##########################################
suppressMessages(library(DECIPHER))
suppressMessages(library(ape))
packageVersion("DECIPHER")

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
# list of vectors for each contig in the genome
l_tr_pos <- list()
for (i in seq_along(contig_lengths)) {
	l_tr_pos[[names(contig_lengths)[i]]] <- numeric(contig_lengths[i]) # 0 to all positions
}
l_CDS_pos <- l_tr_pos

# read gff
filename_gff <- paste(seqID,'.gff.gz',sep='')
gff <- read.gff(filename_gff, na.strings = c(".", "?"), GFF3 = TRUE)
CDS <- gff[gff$type=='CDS' & grepl("protein_id=", gff$attributes, fixed = TRUE),] # non-pseudogene CDS
rm(gff)
if (nrow(CDS)<1) {
	cat('!!!!!!!!!!!!!!! NO CDS in GFF.\n')
	quit(save = "no", status = 1, runLast = FALSE) # terminate script here
}
CDS$protein_id <- gsub(".*protein_id=([^;]+).*", "\\1", CDS$attributes)
CDS$GeneID <- gsub(".*GeneID:([0-9]+).*", "\\1", CDS$attributes)
cat('# protein_id count:', length(unique(CDS$protein_id)), '\n')
cat('# GeneID count:', length(unique(CDS$GeneID)), '\n')
cat('# CDS count:', nrow(CDS), '\n')

# read CDS and mark positions as 1
for (row_i in rownames(CDS)) {
	l_CDS_pos[[as.character(CDS[row_i, 'seqid'])]][c(as.numeric(CDS[row_i, 'start']):as.numeric(CDS[row_i, 'end']))] <- 1
}
cat('# CDS position count:', sapply(l_CDS_pos, sum), '\n') # list of each chromosome

# read TR and mark positions as 1
TR_all <- readRDS(paste(seqID,'.tr',sep=''))
for (row_i in rownames(TR_all)) {
	l_tr_pos[[TR_all[row_i, 'Index']]][c(TR_all[row_i, 'Begin']:TR_all[row_i, 'End'])] <- 1
}
cat('# TR position number:', sapply(l_tr_pos, sum), '\n')
cat('############################################################\n\n')



############## generate matrix for each CDS ##############
num_protein_coding_gene_proc <- 0
sum_TR_199_start <- numeric(199)
sum_pos_199_start <- numeric(199)
sum_TR_199_end <- numeric(199)
sum_pos_199_end <- numeric(199)

# loop through GeneID
GeneIDs <- unique(CDS$GeneID)
for (GeneID_h in GeneIDs) {
	#### CDS info ####
	CDSs_h <- CDS[CDS$GeneID==GeneID_h,]
	CDSs_h <- CDSs_h[order(CDSs_h$start), ] # sort by position
	rownames(CDSs_h) <- seq(nrow(CDSs_h))
	left_h <- CDSs_h[1,'start']
	right_h <- CDSs_h[nrow(CDSs_h),'end']
	contig_name_h <- as.character(CDSs_h[1, 'seqid'])
	strand <- as.character(CDSs_h[1, 'strand'])
	contig_len <- contig_lengths[contig_name_h]
	for (check_end in c(FALSE, TRUE)) {
		# left or right border
		center_position <- ifelse(((strand=='+' & check_end==FALSE) | (strand=='-' & check_end==TRUE)), left_h, right_h) 
		TR_region_199 <- numeric(199)
		CDS_region_199 <- rep(NA, 199)
		if (center_position > contig_len) { # skip this CDS
			cat('# center_position exceed contig length:', GeneID_h, center_position, '>', contig_len, '\n')            
		} else {
			g_start_index <- ifelse(center_position > 99, center_position-99, 1)
			v_start_index <- ifelse(center_position > 99, 1, 101-center_position)
			g_end_index <- ifelse(contig_len-center_position >= 99, center_position+99, contig_len)
			v_end_index <- ifelse(contig_len-center_position >= 99, 199, contig_len-center_position+100)
			if (length(g_start_index:g_end_index)==length(v_start_index:v_end_index)) {
				TR_region_199[v_start_index:v_end_index] <- l_tr_pos[[contig_name_h]][g_start_index:g_end_index]
				CDS_region_199[v_start_index:v_end_index] <- l_CDS_pos[[contig_name_h]][g_start_index:g_end_index]
				num_protein_coding_gene_proc <- num_protein_coding_gene_proc + 1
			} else {
				cat('!!!!!!!!!!something wrong with row', GeneID_h, '\n')
				quit(save = "no", status = 1, runLast = FALSE)
			}
			# flip if strand is reverse "-"
			if (strand == '-') {
				TR_region_199 <- rev(TR_region_199)
				CDS_region_199 <- rev(CDS_region_199)
			}
			# calculate possible region: invert upstream postions, ignore NA
			possible_region_199 <- CDS_region_199
			if (check_end==FALSE) { # invert upstream/downstream 1 -> 0 ; 0 -> 1
				possible_region_199[c(1:99)] <- (possible_region_199[c(1:99)]-1)*-1
				possible_region_199[is.na(possible_region_199)] <- 0
				trueTR_region_199 <- possible_region_199*TR_region_199
				sum_TR_199_start <- sum_TR_199_start + trueTR_region_199
				sum_pos_199_start <- sum_pos_199_start + possible_region_199
			} else {
				possible_region_199[c(101:199)] <- (possible_region_199[c(101:199)]-1)*-1
				possible_region_199[is.na(possible_region_199)] <- 0
				trueTR_region_199 <- possible_region_199*TR_region_199
				sum_TR_199_end <- sum_TR_199_end + trueTR_region_199
				sum_pos_199_end <- sum_pos_199_end + possible_region_199
			}
		}
	}
}
cat('# Number of protein coding gene processed:', num_protein_coding_gene_proc, '\n') 


############## Write output ##############
out_plot <- paste(seqID,'_out.png',sep='')

norm_TR_start <- sum_TR_199_start/sum_pos_199_start
norm_TR_end <- sum_TR_199_end/sum_pos_199_end
saveRDS(norm_TR_start, file=paste(seqID,'_start.rds',sep=''), compress = TRUE)
saveRDS(norm_TR_end, file=paste(seqID,'_end.rds',sep=''), compress = TRUE)
cat('###norm_TR_start:', seqID, norm_TR_start, '\n', sep=',')
cat('###norm_TR_end:', seqID, norm_TR_end, '\n', sep=',')

# save plots as png
show_labels <- c(1,50, 65, 100, 149,199)
png(filename=out_plot)

par(mfrow = c(6, 1), mar=c(2,5,2,1))

barplot(height=sum_pos_199_start, ylab='Possibles',space=0, col="gray", border="gray", main='start')
bar_center <- barplot(sum_pos_199_start,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

barplot(height=sum_TR_199_start, ylab='TR',space=0, col="gray", border="gray")
bar_center <- barplot(sum_TR_199_start,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

barplot(height=norm_TR_start, ylab='Normalized TR',space=0, col="gray", border="gray")
bar_center <- barplot(norm_TR_start,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

barplot(height=sum_pos_199_end, ylab='Possibles',space=0, col="gray", border="gray", main='end')
bar_center <- barplot(sum_pos_199_end,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

barplot(height=sum_TR_199_end, ylab='TR',space=0, col="gray", border="gray")
bar_center <- barplot(sum_TR_199_end,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

barplot(height=norm_TR_end, ylab='Normalized TR',space=0, col="gray", border="gray")
bar_center <- barplot(norm_TR_end,space=0, plot = FALSE)
axis(side = 1, at = bar_center[show_labels], labels = c('-99','-50','-35','+1','+50','+100'))

mtext(paste(seqID, ' (', nrow(GeneIDs),' GeneIDs)', sep=''), side = 3, line = -2, outer = TRUE)
dev.off()

cat('Save plots as', out_plot, '\n') 