library(DECIPHER)

path <- "~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/repeat_alignments_v2/"

ls <- list.files(path)

sM <- matrix(c(0.951, -0.266, -0.305, -0.327, 0.099, -0.126, -0.151, 0.041, -0.383, -0.26, -0.274, -0.236, -0.131, -0.451, -0.14, 0.191, 0.008, -0.586, -0.501, 0.015, -2.542, -0.266, 1.465, 0.013, -0.153, -0.741, 0.37, 0.117, -0.453, 0.178, -0.81, -0.701, 0.679, -0.46, -0.874, -0.311, -0.097, -0.129, -0.5, -0.415, -0.679, -2.542, -0.305, 0.013, 1.494, 0.532, -0.582, 0.189, 0.129, 0.037, 0.255, -0.962, -0.947, 0.202, -0.606, -0.879, -0.252, 0.261, 0.075, -0.873, -0.433, -0.834, -2.542, -0.327, -0.153, 0.532, 1.575, -1.007, 0.155, 0.597, -0.131, -0.045, -1.265, -1.194, 0.052, -0.915, -1.235, -0.131, 0.099, -0.121, -1.102, -0.8, -1.051, -2.542, 0.099, -0.741, -0.582, -1.007, 3.127, -0.777, -0.995, -0.499, -0.438, -0.174, -0.218, -0.884, -0.137, -0.189, -0.832, -0.091, -0.185, -0.446, -0.268, 0.016, -2.542, -0.126, 0.37, 0.189, 0.155, -0.777, 1.289, 0.494, -0.368, 0.251, -0.762, -0.636, 0.433, -0.259, -0.846, -0.241, 0.046, -0.01, -0.707, -0.444, -0.624, -2.542, -0.151, 0.117, 0.129, 0.597, -0.995, 0.494, 1.287, -0.38, -0.057, -0.967, -0.931, 0.342, -0.646, -1.116, -0.162, 0.007, -0.072, -0.97, -0.681, -0.758, -2.542, 0.041, -0.453, 0.037, -0.131, -0.499, -0.368, -0.38, 1.768, -0.42, -1.087, -1.022, -0.369, -0.758, -0.924, -0.333, 0.043, -0.342, -0.886, -0.863, -0.861, -2.542, -0.383, 0.178, 0.255, -0.045, -0.438, 0.251, -0.057, -0.42, 2.254, -0.781, -0.663, 0.033, -0.433, -0.29, -0.354, -0.099, -0.206, -0.217, 0.381, -0.664, -2.542, -0.26, -0.81, -0.962, -1.265, -0.174, -0.762, -0.967, -1.087, -0.781, 1.184, 0.585, -0.819, 0.423, 0.216, -0.799, -0.716, -0.29, -0.347, -0.258, 0.785, -2.542, -0.274, -0.701, -0.947, -1.194, -0.218, -0.636, -0.931, -1.022, -0.663, 0.585, 1.087, -0.799, 0.585, 0.39, -0.777, -0.73, -0.43, -0.123, -0.159, 0.343, -2.542, -0.236, 0.679, 0.202, 0.052, -0.884, 0.433, 0.342, -0.369, 0.033, -0.819, -0.799, 1.282, -0.5, -1.005, -0.175, 0.006, -0.035, -0.829, -0.564, -0.704, -2.542, -0.131, -0.46, -0.606, -0.915, -0.137, -0.259, -0.646, -0.758, -0.433, 0.423, 0.585, -0.5, 1.637, 0.285, -0.712, -0.406, -0.171, -0.135, -0.091, 0.219, -2.542, -0.451, -0.874, -0.879, -1.235, -0.189, -0.846, -1.116, -0.924, -0.29, 0.216, 0.39, -1.005, 0.285, 1.717, -0.837, -0.7, -0.528, 0.608, 0.885, 0.045, -2.542, -0.14, -0.311, -0.252, -0.131, -0.832, -0.241, -0.162, -0.333, -0.354, -0.799, -0.777, -0.175, -0.712, -0.837, 2.121, -0.015, -0.198, -0.777, -0.763, -0.588, -2.542, 0.191, -0.097, 0.261, 0.099, -0.091, 0.046, 0.007, 0.043, -0.099, -0.716, -0.73, 0.006, -0.406, -0.7, -0.015, 0.979, 0.427, -0.727, -0.505, -0.505, -2.542, 0.008, -0.129, 0.075, -0.121, -0.185, -0.01, -0.072, -0.342, -0.206, -0.29, -0.43, -0.035, -0.171, -0.528, -0.198, 0.427, 1.128, -0.659, -0.439, -0.062, -2.542, -0.586, -0.5, -0.873, -1.102, -0.446, -0.707, -0.97, -0.886, -0.217, -0.347, -0.123, -0.829, -0.135, 0.608, -0.777, -0.727, -0.659, 3.153, 0.763, -0.428, -2.542, -0.501, -0.415, -0.433, -0.8, -0.268, -0.444, -0.681, -0.863, 0.381, -0.258, -0.159, -0.564, -0.091, 0.885, -0.763, -0.505, -0.439, 0.763, 2.023, -0.287, -2.542, 0.015, -0.679, -0.834, -1.051, 0.016, -0.624, -0.758, -0.861, -0.664, 0.785, 0.343, -0.704, 0.219, 0.045, -0.588, -0.505, -0.062, -0.428, -0.287, 1.084, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, -2.542, 3.235),
	nrow=21,
	ncol=21,
	dimnames=list(c(AA_STANDARD, "*"), c(AA_STANDARD, "*")))

scores <- scores2 <- d <- d2 <- widths <- numeric(length(ls))
l <- integer(length(ls))
AAs <- matrix(0, length(ls), 20)
colnames(AAs) <- AA_STANDARD
pBar <- txtProgressBar(style=3)
for (i in seq_along(scores)) {
	setTxtProgressBar(pBar, i/length(scores))
	aa <- readAAStringSet(paste0(path, ls[i]))
	l[i] <- length(aa)
	aa <- AAStringSet(toupper(aa))
	AAs[i,] <- alphabetFrequency(unlist(aa))[1:20]
	d[i] <- mean(DistanceMatrix(aa, verbose=FALSE, type="dist"))
	scores[i] <- ScoreAlignment(aa,
		method="adjacent",
		gapOpening=-0.1,
		gapExtension=-2.6,
		substitutionMatrix=sM,
		includeTerminalGaps=TRUE)
	# realign
	aa <- RemoveGaps(aa)
	widths[i] <- mean(width(aa))
	AA <- AlignSeqs(aa, verbose=FALSE)
	d2[i] <- mean(DistanceMatrix(AA, verbose=FALSE, type="dist"))
	scores2[i] <- ScoreAlignment(AA,
		method="adjacent",
		gapOpening=-0.1,
		gapExtension=-2.6,
		substitutionMatrix=sM,
		includeTerminalGaps=TRUE)
}
close(pBar)

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

writeLines(substring(ls, 1, nchar(ls) - 4L)[W], "~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/PDB_TrainingSet_v4.txt")

### Single Amino Acids

dev.new()
z <- readAAStringSet("~/Desktop/Sequencing/Assembly/TandemRepeats/cesymm-2.2.0-SNAPSHOT/cullpdb_pc50.0_res0.0-2.0_noBrks_len40-10000_R0.25_Xray_d2025_02_20_chains16323.fasta")
m <- match(substring(ls, 1, nchar(ls) - 4L), gsub("^(.+?) .+$", "\\1", names(z)))
a1 <- alphabetFrequency(z[m])[, 1:20]
a1 <- colSums(a1[W,]) # background
a2 <- colSums(AAs[W,]) # foreground
a1 <- a1 - a2 # approximately subtract repeat region
a1 <- a1/sum(a1)
a2 <- a2/sum(a2)
S <- log(a2/a1)
plot(S, pch=NA, xlab="Amino acid", ylab="Log-odds in repeat")
text(S, colnames(AAs))
mean(colSums(S*t(AAs)))
cat(paste(colnames(AAs), round(S, 2), sep="="), sep=", ")

dev.new()
matplot(t(AAs/rowSums(AAs)), type="l", col="#00000001", lty=1, xlab="Amino acid", ylab="Relative frequency", main="Compositional bias")

### Periodicities

dev.new()
t <- tabulate(widths)/length(widths)
plot(t, xlab="Periodicity", ylab="Relative frequency", log="y") # power-law tail

dev.new()
Y <- cumsum(t)
plot(Y, log="y")
SAE <- function(p, plot=FALSE) {
	y <- (p[1] + p[2]*p[6]^(-p[3]*seq_along(t - p[5])))^(-1/p[4])
	if (plot)
		lines(y, col="red")
	sum(abs(log(y) - log(Y)))
}
o <- optim(c(1, 10, 0.01, 1, 0, exp(1)), SAE) # fit sigmoid to CDF
SAE(o$par, TRUE)

dev.new()
layout(matrix(1:2))
freqs <- colMeans(AAs/rowSums(AAs))
set.seed(123L)
rand <- AAStringSet(paste(sample(names(freqs), 1e7, replace=TRUE, prob=freqs), collapse=""))
x <- DetectRepeats(rand, useEmpirical=FALSE, minScore=0)
trand <- mapply(function(x, y) {
		mean(y - x + 1)
	},
	x[, "Left"],
	x[, "Right"])
breaks <- c(0, 7:30, seq(32, 50, 2), seq(55, 80, 5), seq(100, 200, 40), 500)
fg <- tabulate(.bincode(widths[W], breaks))
bg <- tabulate(.bincode(trand, breaks))
fg <- fg/sum(fg)
bg <- bg/sum(bg)
logodds <- log(fg/bg)

SAE <- function(p, plot=FALSE) {
	est <- p[1] + (p[2] - p[1])/(p[3] + p[4]*exp(p[5]*X))^(1/p[6])
	est[is.infinite(est)] <- p[1]
	if (plot)
		lines(X, est, col="red")
	sum(abs(est - c(logodds[1], logodds)))
}

X <- breaks
o <- optim(c(0.5, -4, 1, 0.01, 1, 10), SAE) # fit sigmoid to CDF
plot(breaks[-1], logodds, xlab="Periodicity", ylab="Log-odds", xlim=c(0, 100))
SAE(o$par, TRUE)
cat(o$par, sep=", ")

X <- breaks*3 # convert to nucleotides
o <- optim(c(0.5, -4, 1, 0.01, 1, 10), SAE) # fit sigmoid to CDF
plot(3*breaks[-1], logodds, xlab="Periodicity", ylab="Log-odds", xlim=c(0, 100))
SAE(o$par, TRUE)
cat(o$par, sep=", ")

### Copy Number of Repeat

dev.new()
l1 <- tabulate(l[W])
l1 <- l1/sum(l1)
l2 <- tabulate(lengths(x[, "Left"]), length(l1))
l2 <- l2/sum(l2)
subset <- 1:9 # NOTE: repeats start at 2 copies
l1 <- c(l1[subset], sum(l1[-subset]))
l2 <- c(l2[subset], sum(l2[-subset]))
lens <- log(l1/l2)
plot(lens, xlab="Number of repeats", ylab="Log-odds", pch=NA)
text(seq_along(lens), lens)
abline(lm(lens~seq_along(lens)), lty=2)
cat(round(lens, 1), sep=", ")
