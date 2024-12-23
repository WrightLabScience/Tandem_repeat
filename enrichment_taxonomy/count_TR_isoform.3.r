# Usage: 
# Rscript count_TR_isoform.3.r GCA_932276165 GCA_932276165.tr protein.faa genomic.gff

# read arguments
ARGS <- commandArgs(trailingOnly = TRUE)
seqID <- ARGS[1]
TR_filename <- ARGS[2]
seqs_filename <- ARGS[3]
gff_filename <- ARGS[4]

print(ARGS)

# load libraries
suppressMessages(library(DECIPHER))
suppressMessages(library(ape))

packageVersion("DECIPHER")

options(timeout=999999999) # default is 60 sec, will fail if the genome is too big

######################## main ########################
# read files
seqs <- readAAStringSet(seqs_filename)
TR_all <- readRDS(TR_filename)
gff <- read.gff(gff_filename, na.strings = c(".", "?"), GFF3 = TRUE)

# keep protein_id
seq_names <- names(seqs)
seq_names_clean <- gsub("(^[^ ]+) .*", "\\1", seq_names)

# convert gff to protein_id, GeneID
CDS <- gff[gff$type=='CDS',]
attributes <- CDS$attributes[grepl("protein_id=", CDS$attributes, fixed = TRUE)] # exclude pseudogenes
protein_id_gff <- gsub(".*protein_id=([^;]+).*", "\\1", attributes)
GeneID_gff <- gsub(".*GeneID:([0-9]+).*", "\\1", attributes)
# lookup table
df <- data.frame(protein_id_gff, GeneID_gff) 
df <- df[!duplicated(df), ]
df <- df[order(df$protein_id_gff),]
rownames(df) <- seq(nrow(df))
df$TR <- 0

vec_uniq_GeneID <- unique(df$GeneID_gff)

cat('### faa seqs count: ', length(seqs), '\n')
cat('### gff protein_id count: ', length(unique(df$protein_id_gff)), '\n')
cat('### gff GeneID count: ', length(vec_uniq_GeneID), '\n')
cat('### gff not in faa: ', length(df$protein_id_gff[!is.element(df$protein_id_gff, seq_names_clean)]), '\n')
cat('### faa not in gff: ', length(seq_names_clean[!is.element(seq_names_clean, df$protein_id_gff)]), '\n')

GeneID_in_seqs_n <- 0
uniq_seq_w_TR_n <- 0

for (GeneID_h in vec_uniq_GeneID) {
    protein_ids <- df[df$GeneID_gff==GeneID_h,'protein_id_gff']
    protein_ids_Index <- which(is.element(seq_names_clean, protein_ids))
    if (length(protein_ids_Index) > 0) {
        GeneID_in_seqs_n <- GeneID_in_seqs_n + 1
        TR_h <- TR_all[is.element(TR_all$Index, protein_ids_Index),]
        if (nrow(TR_h)>0) {
            uniq_seq_w_TR_n <- uniq_seq_w_TR_n + 1
        }
    }
}

cat('### Report uniq:  seqID =', seqID, ', uniq_seq_n =', GeneID_in_seqs_n, ', uniq_TR_n =', uniq_seq_w_TR_n,', frac_uniq_TR =', uniq_seq_w_TR_n/GeneID_in_seqs_n, '\n')
cat('### Report total: seqID =', seqID, ', total_seq_n =', length(seqs), ', total_TR_n =', length(unique((TR_all$Index))), ', frac_total_TR =', length(unique((TR_all$Index)))/length(seqs), '\n')

saveRDS(df, paste(seqID, '.df_TR', sep=''))