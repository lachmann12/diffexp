setwd("~/code/diffexpression")


if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")

BiocManager::install("biomaRt")
BiocManager::install("Mus.musculus")

install.packages("RColorBrewer")



library("edgeR")
library("limma") 
library("RColorBrewer")
library("Mus.musculus")


# http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file
files <- c("data/GSM1545535_10_6_5_11.txt", "data/GSM1545536_9_6_5_11.txt", 
   "data/GSM1545538_purep53.txt", "data/GSM1545539_JMS8-2.txt", 
   "data/GSM1545540_JMS8-3.txt", "data/GSM1545541_JMS8-4.txt", 
   "data/GSM1545542_JMS8-5.txt", "data/GSM1545544_JMS9-P7c.txt", 
   "data/GSM1545545_JMS9-P8c.txt") 
read.delim(files[1], nrow=5)

x <- readDGE(files, columns=c(1,3)) 
dim(x)

samplenames <- substring(colnames(x), 12, nchar(colnames(x))) 

colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP")) 
x$samples$group <- group 
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2))) 
x$samples$lane <- lane 
x$samples

geneid <- rownames(x) 
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), keytype="ENTREZID")
dim(genes)
head(genes)
genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes



cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)

L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6

pdf("fig/cut_clean.pdf", 10, 6)
lcpm.cutoff <- log2(10/M + 2/L)

nsamples <- ncol(x)
col <- brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
lcpm <- cpm(x, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
legend("topright", samplenames, text.col=col, bty="n")
dev.off()


keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)

x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors

x2 <- x
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5

pdf("fig/barplots.pdf", 10, 6)
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data", ylab="Log-cpm")
x2 <- calcNormFactors(x2)
x2$samples$norm.factors

lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data", ylab="Log-cpm")
dev.off()

lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")

plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")

design <- model.matrix(~0+group+lane)
colnames(design) <- gsub("group", "", colnames(design))
design

contr.matrix <- makeContrasts(
   BasalvsLP = Basal-LP,
   BasalvsML = Basal - ML, 
   LPvsML = LP - ML,
   levels = colnames(design))
contr.matrix

v <- voom(x, design, plot=TRUE)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit)

summary(decideTests(efit))

de.common <- which(dt[,1]!=0 & dt[,2]!=0)
length(de.common)

pdf("fig/venn_diag.pdf")
vennDiagram(dt[,1:2], circle.col=c("turquoise", "salmon"))
dev.off()

write.fit(tfit, dt, file="results.txt")

basal.vs.lp <- topTreat(tfit, coef=1, n=Inf)
basal.vs.ml <- topTreat(tfit, coef=2, n=Inf)
head(basal.vs.lp)

pdf("fig/md_plot.pdf")
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
dev.off()

