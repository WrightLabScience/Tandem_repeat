library(DECIPHER)

path <- "~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/repeat_alignments_v2/"

ls <- list.files(path)

# subset to those sets not used in benchmarking
r <- readLines("~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/CESYMM_seqID_UsedInBenchmarking_v2.txt")
W <- !(substring(ls, 1, nchar(ls) - 4L) %in% r) # set not used in benchmarking

# subset to those sets with non-zero TM-scores
load('~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/ResultPDBs_v2.RData')
res <- unlist(res)
res <- strsplit(res, "\t")
w <- which(lengths(res) == 16) # missing repeats
res[w] <- lapply(res[w], c, NA_character_)
res <- do.call(cbind, res)
res <- t(res)
res <- as.data.frame(res)
res <- res[res[, "V11"] >= 0.4 & res[, "V11"] < 0.6,]
W[!(substring(ls, 1, nchar(ls) - 4L) %in% rownames(res))] <- FALSE

ls <- ls[W]

seqs <- readAAStringSet("~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/cullpdb_pc50.0_res0.0-2.0_noBrks_len40-10000_R0.25_Xray_d2025_02_20_chains16323.fasta")
m <- match(sapply(strsplit(ls, ".", fixed=TRUE), `[`, 1L),
	sapply(strsplit(names(seqs), " ", fixed=TRUE), `[`, 1L))
seqs <- seqs[m]

L <- 101L # number of points in interpolation
res1 <- numeric(L)
pBar <- txtProgressBar(style=3)
for (i in seq_along(ls)) {
	setTxtProgressBar(pBar, i/length(ls))
	aa <- readBStringSet(paste0(path, ls[i]))
	AA <- AAStringSet(toupper(aa))
	AA <- RemoveGaps(AA)
	
	# add intervening region between units back to sequences if possible
	B <- gsub(paste0("^.*(", paste0("(", AA, "", collapse=".*?)"), ")).*$"),
		"\\1",
		seqs[i])
	BB <- character(length(AA))
	for (j in seq_along(AA)) {
		BB[j] <- gsub(paste0("^(.*", AA[j], ").*$"), "\\1", B)
		B <- gsub(paste0("^(.*", AA[j], ")(.*)$"), "\\2", B)
	}
	BB <- AAStringSet(BB)
	if (length(BB) == length(AA) && all(width(BB) >= width(AA)))
		AA <- BB # replace sequences with full region
	
	AA <- AlignSeqs(AA, verbose=FALSE, anchor=NA, iterations=0, refinements=0) # realign
	pwm <- .Call("consensusProfileAA",
		AA,
		rep(1, length(AA)),
		NULL,
		PACKAGE="DECIPHER")
	a <- pwm[27,] # column occupancy
	res1 <- res1 + a[round(seq(1, length(a), length.out=L))]
}
close(pBar)

r1 <- res1/length(ls) # convert to average

# weight function
w <- 1 + -0.3*abs(seq(-1, 1, length.out=length(res1)))^3 # optimized function
w <- w/mean(w)

off <- 0.12
dev.new(height=3.5, width=4)
par(mar=c(4, 4, 2, 4))
plot(seq(0, 100, length.out=L),
	r1,
	ylim=range(r1) + c(-0.2, 0.15),
	type="l",
	xlab="",
	ylab="",
	yaxt="n")
axis(2, at=axTicks(2), labels=100*axTicks(2))
mtext(side=1, "Relative alignment position (%)", line=2)
mtext(side=2, "Average column occupancy (%)", line=2)
axis(4, at=axTicks(2) + (round(off, 1) - off), labels=axTicks(2) + round(off, 1), col="gray", col.axis="gray")
mtext(side=4, "Relative column weight", line=2, col="gray")
lines(w - off, col="gray", lty=2)
legend(x=-15,
	y=1.3,
	ncol=2,
	xpd=TRUE,
	bty="n",
	legend=c("Observed", "Optimized"),
	lty=1:2,
	text.col=c("black", "darkgray"),
	col=c("black", "darkgray"))
