Analysis_enrichment_function.txt

Shu-Ting Cho <shutingcho@pitt.edu>

# Tandem repeat enrichment analysis for functional groups
# Fig4B

# Step0: Prepare source data
# Step1: Classify protein function
# Step2: Run DetectRepeats
# Step4: Compute values and plot

# load library
library(DECIPHER)

# Step 1.1: Classify protein function
# Training data for IdTaxa downloaded from DECIPHER website
# Training data was selected based on genome's taxonomy
# Use GCA_932276165 genome as an example
# http://www2.decipher.codes/Classification/TrainingSets/KEGG_Animals_r95.RData
Rscript DECIPHER_IdTaxa.7.r GCA_932276165 KEGG_Animals_r95.RData

# collect Kgroup aa fasta from each genome
Rscript collect_Kgroup_seq.1.r GCA_932276165

mkdir Kgroup_seq
mv GCA_*.tar.gz Kgroup_seq
tar -zcvf Kgroup_seq.tar.gz Kgroup_seq/

# concat Kgroup seq from all genomes
Rscript write_Kgroup_faa.r K00001
# result files: euk_Kgroup_faa/*.faa.gz prok_Kgroup_faa/*.faa.gz


# Step2: Run DetectRepeats
Rscript DECIPHER_DetectRepeats.19.r K00001 1 100



# Step 3: Calculate values for plotting
# values needed:
# x-axis: Fraction of genomes have the Kgroup (perc_genomes)
# y-axis: Fraction of Kgroup have tandem repeat (perc_TR_genes)
# dot size: DetectRepeats average score (Scores_highest_ave) 
# shade: Average sequence identity of genes in Kgroup (PID)

# Use prokaryote as an example

# Step3.1: x-axis: Fraction of genomes have the Kgroup (perc_genomes)
load('Kgroup_mt.RData')
prok_Kgroup_Brite_TR <- readRDS('prok_Kgroup_Brite_TR.rds')
inDir <- 'pro_Kgroup_faa/'
Kgroups <- rownames(prok_Kgroup_Brite_TR)

prok4028_mt_Kgroups_1 <- pro4028_mt_Kgroups
colnames(prok4028_mt_Kgroups_1) <- gsub("^(K[0-9]{5}) .*", "\\1", colnames(prok4028_mt_Kgroups_1))

prok_Kgroup_Brite_TR$n_genomes <- 0
prok_Kgroup_Brite_TR$perc_genomes <- 0
for (i in seq_along(Kgroups)) {
	Kgroup <- Kgroups[i]
	filename <- paste(inDir, Kgroup, '.faa.gz', sep='')
	prok_Kgroup_Brite_TR[Kgroup,'n_genomes'] <- sum(prok4028_mt_Kgroups_1[,Kgroup]>0, na.rm=TRUE)
	prok_Kgroup_Brite_TR[Kgroup,'perc_genomes'] <- mean(prok4028_mt_Kgroups_1[,Kgroup]>0, na.rm=TRUE)
	cat(i, Kgroup, prok_Kgroup_Brite_TR[Kgroup,'n_genomes'], prok_Kgroup_Brite_TR[Kgroup,'perc_genomes'], '\n', sep='\t')
}

# Step3.2: shade: Average sequence identity of genes in Kgroup (PID)
# We randomly sampled up to 1000 sequences in each KO group to determine the group's average percent identity
Rscript calculate_sample_PID.1.r K20276.faa.gz K20276 1000
grep 'Result\ Report\ PID' *.out | sed 's/.*:Result Report PID:\t//g' > prok_PID.txt

PID_tab <- read.table(file = 'prok_PID.txt', header = FALSE, row.names = 1, sep= "\t")
prok_Kgroup_Brite_TR <- prok_Kgroup_Brite_TR[order(row.names(prok_Kgroup_Brite_TR)), ]
PID_tab <- PID_tab[order(row.names(PID_tab)), ]
prok_Kgroup_Brite_TR$PID <- PID_tab$V2


# Step3.3:
# y-axis: Fraction of Kgroup have tandem repeat (perc_TR_genes)
# dot size: DetectRepeats average score (Scores_highest_ave) 
prok_Kgroup_Brite_TR <- prok_Kgroup_Brite_TR[,c(1:6, 17, 20, 21, 22)]
prok_Kgroup_Brite_TR$n_TR_genes <- 0
prok_Kgroup_Brite_TR$Scores_highest_ave <- 0
for (Kgroup_h in prok_Kgroup_Brite_TR$KgroupID) {
	print(Kgroup_h)
	TR_h <- readRDS(paste('DetectRepeats_Kgroup_prok/', Kgroup_h, '.tr', sep=''))
	if (nrow(TR_h) > 0) {
		# for perc_TR_genes
		prok_Kgroup_Brite_TR[Kgroup_h,'n_TR_genes'] <- length(unique(TR_h$Index))
		# Scores_highest_ave
		scores_h <- c()
		for (seq_index in unique(TR_h$Index)) {
			scores_h <- c(scores_h, max(TR_h[which(TR_h$Index == seq_index), 'Score']))
		}
		prok_Kgroup_Brite_TR[Kgroup_h,'Scores_highest_ave'] <- mean(scores_h)
	}
}
prok_Kgroup_Brite_TR$perc_TR_genes <- prok_Kgroup_Brite_TR$n_TR_genes/prok_Kgroup_Brite_TR$n_total_genes

# label text
prok_Kgroup_Brite_TR$label_text <- prok_Kgroup_Brite_TR$KgroupID
prok_Kgroup_Brite_TR[prok_Kgroup_Brite_TR$perc_TR_genes<0.5, 'label_text'] <- ''
prok_Kgroup_Brite_TR[prok_Kgroup_Brite_TR$perc_genomes<0.8, 'label_text'] <- ''



# Step3.4: Scatter plot
library(ggplot2)
library(ggrepel)

# sort data points by PID: darker points (high PID) will be in the foreground
prok_Kgroup_Brite_TR_sortPID <- prok_Kgroup_Brite_TR[order(prok_Kgroup_Brite_TR$PID),]
euk_Kgroup_Brite_TR_sortPID <- euk_Kgroup_Brite_TR[order(euk_Kgroup_Brite_TR$PID),]

Brite_lv1_names <- c("09100 Metabolism","09120 Genetic Information Processing","09130 Environmental Information Processing","09140 Cellular Processes","09150 Organismal Systems","09160 Human Diseases") 

pdf(file = "BriteLv1_scatter.5.pdf", width = 2.43448819, height = 7.30346457)
par(mfrow = c(6, 2), mar=c(0,0,0,0))
for (i in seq_along(Brite_lv1_names)) {
	group_name <- Brite_lv1_names[i]
	# Scatter plot Prokaryotes
    df_h <- prok_Kgroup_Brite_TR_sortPID[prok_Kgroup_Brite_TR_sortPID$lv1==group_name,]
	df_h <- df_h[df_h$Scores_highest_ave>=10,]
	df_h$perc_TR_genes <- df_h$perc_TR_genes
	df_h$perc_genomes <- df_h$perc_genomes
	if (i==6) {
		plot(y=df_h$perc_TR_genes, x=df_h$perc_genomes, 
		pch = 16,
		col = df_h$PID_colors,
		cex = sqrt(df_h$Scores_highest_ave)/15,
		xlim=c(0, 1), 
		ylim=c(0, 1),
		ann=FALSE)
	} else {
		plot(y=df_h$perc_TR_genes, x=df_h$perc_genomes, 
		xaxt='n',
		pch = 16,
		col = df_h$PID_colors,
		cex = sqrt(df_h$Scores_highest_ave)/15,
		xlim=c(0, 1), 
		ylim=c(0, 1),
		ann=FALSE)
	}
	low <- lowess(df_h$perc_genomes, df_h$perc_TR_genes, f=1)
	lines(low, lty = 1, col = 3, lwd = 2)

	# Scatter plot Eukaryotes
	df_h <- euk_Kgroup_Brite_TR_sortPID[euk_Kgroup_Brite_TR_sortPID$lv1==group_name,]
	df_h <- df_h[df_h$Scores_highest_ave>=10,]
	df_h$perc_TR_genes <- df_h$perc_TR_genes
	df_h$perc_genomes <- df_h$perc_genomes
	if (i==6) {
		plot(y=df_h$perc_TR_genes, x=df_h$perc_genomes, 
		yaxt='n',
		pch = 16,
		col = df_h$PID_colors,
		cex = sqrt(df_h$Scores_highest_ave)/15,
		xlim=c(0, 1), 
		ylim=c(0, 1),
		ann=FALSE)
	} else {
		plot(y=df_h$perc_TR_genes, x=df_h$perc_genomes, 
		xaxt='n',
		yaxt='n',
		pch = 16,
		col = df_h$PID_colors,
		cex = sqrt(df_h$Scores_highest_ave)/15,
		xlim=c(0, 1), 
		ylim=c(0, 1),
		ann=FALSE)
	}
	low <- lowess(df_h$perc_genomes, df_h$perc_TR_genes, f=1)
	lines(low, lty = 1, col = 3, lwd = 2)
}
dev.off()