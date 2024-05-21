library(recount3)
library(recount)
library(edgeR)

# Save the 3 files inside a RSE data object
rse_brain <- readRDS("rse_brain.RDS")
rse_colon <- readRDS("rse_colon.RDS")
rse_muscle <- readRDS("rse_muscle.RDS")

# Check on RIN
colData(rse_brain)$gtex.smrin[33]
colData(rse_brain)$gtex.smrin[34]
colData(rse_brain)$gtex.smrin[38]

colData(rse_colon)$gtex.smrin[32]
colData(rse_colon)$gtex.smrin[33]
colData(rse_colon)$gtex.smrin[34]

colData(rse_muscle)$gtex.smrin[32]
colData(rse_muscle)$gtex.smrin[33]
colData(rse_muscle)$gtex.smrin[34]

# Check on % of rRNA
colData(rse_brain)$gtex.smrrnart[33]
colData(rse_brain)$gtex.smrrnart[34]
colData(rse_brain)$gtex.smrrnart[38]

colData(rse_colon)$gtex.smrrnart[32]
colData(rse_colon)$gtex.smrrnart[33]
colData(rse_colon)$gtex.smrrnart[34]

colData(rse_muscle)$gtex.smrrnart[32]
colData(rse_muscle)$gtex.smrrnart[33]
colData(rse_muscle)$gtex.smrrnart[34]

# Check on % of uniquely mapped reads
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[38]

colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[32]
colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]

colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[32]
colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]

# Conversion of "raw counts" into the actual raw counts
assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_colon)$counts <- transform_counts(rse_colon)
assays(rse_muscle)$counts <- transform_counts(rse_muscle)

# Selection of the 3 replicates from each condition
rse_brain_selected <- rse_brain[,c(33,34,38)]
counts_brain_selected <- assays(rse_brain_selected)$counts

rse_colon_selected <- rse_colon[,c(32,33,34)]
counts_colon_selected <- assays(rse_colon_selected)$counts

rse_muscle_selected <- rse_muscle[,c(32,33,34)]
counts_muscle_selected <- assays(rse_muscle_selected)$counts


x <- cbind(counts_brain_selected, counts_colon_selected, counts_muscle_selected)

colnames(x) <- c("Brain33", "Brain34","Brain38","Colon32", "Colon33","Colon34",
                 "Muscle32","Muscle33","Muscle34")

# Since all the three tables contain exactly the same genes, 
# this can be done only once on the final table
rownames(x) <- rowData(rse_brain_selected)$gene_name

y <- DGEList(counts=x)


group <- as.factor(c("Brain","Brain","Brain","Colon","Colon","Colon",
                     "Muscle","Muscle","Muscle"))
y$samples$group <- group

# Adding interesting informations to 'y'
y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_colon_selected)$gtex.smrin,colData(rse_muscle_selected)$gtex.smrin))

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_colon_selected)$gtex.smtsd,colData(rse_muscle_selected)$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_colon_selected)$gtex.sex,colData(rse_muscle_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_colon_selected)$gtex.age,colData(rse_muscle_selected)$gtex.age))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_colon_selected)$gtex.smrrnart,colData(rse_muscle_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_colon_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_muscle_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_colon_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_muscle_selected)$"recount_qc.aligned_reads%.chrm"))

y

# Let's take a look at the library size: usually we want at least 30M/50M as library size
# but it can happen to see just 10M --> the problem is in the 'transform_counts' step,
# and the number I see is just an artifact, is NOT the real library size of each sample!
# (so we ignore this column)

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

logcpm_before <- cpm(y, log=TRUE)
y <- calcNormFactors(y, method = "TMM")
logcpm_after <- cpm(y, log=TRUE)

par(mfrow = c(1,2))
boxplot(logcpm_before, notch=T)
boxplot(logcpm_after, notch=T)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

# Sometimes a colon sample can be found very far from the others colon samples (this
# is NOT our case). In our case a brain sample is strange.
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)

# We look at those plots to try to understand why that brain sample is far from the
# others brain samples (we think that's because of the region of the brain
# from where it comes from --> the other 2 are cortex, it is not).
plotMDS(logcpm, labels=y$samples$rRNA)
plotMDS(logcpm, labels=y$samples$chrm)
plotMDS(logcpm, labels=y$samples$age)
plotMDS(logcpm, labels=y$samples$mapped)
plotMDS(logcpm, labels=y$samples$rin)
plotMDS(logcpm, labels=y$samples$slice)
plotMDS(logcpm, labels=y$samples$sex)
plotMDS(logcpm, labels=rownames(y$samples)) # To have the names of the samples displayed

# The most restrictive thresholds are 0.01 for FDR and 1 for logFC.
# If I think that I obtained too many genes DE, I consider the first 100
# and I look which kind of genes they are, then I look at the first 1000 and so on.

y <- estimateDisp(y, design)
plotBCV(y) # We don't need to comment this plot (just take it as it is)

fit <- glmQLFit(y, design)

# Colon (top) vs Brain (bottom)
qlfCB <- glmQLFTest(fit, contrast=c(-1,1,0))
# Muscle (top) vs Brain (bottom)
qlfMB <- glmQLFTest(fit, contrast=c(-1,0,1))
# Muscle (top) vs Colon (bottom)
qlfMC <- glmQLFTest(fit, contrast=c(0,-1,1))

resultsCB <- topTags(qlfCB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsCB, "resultsCB.txt")

resultsMB <- topTags(qlfMB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMB, "resultsMB.txt")

resultsMC <- topTags(qlfMC, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMC, "resultsMC.txt")

# Summary for the up- and down- regulated genes for the chosen thresholds of p-value and logFC
summary(decideTests(qlfCB, p.value=0.05, adjust.method = "BH", lfc=0))


# !!! Update R to the latest version to execute the commands below
# (or compute the TMP values by hand)
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_colon)$TPM <- recount::getTPM(rse_colon)
assays(rse_muscle)$TPM <- recount::getTPM(rse_muscle)
which(rowData(rse_brain)$gene_name == "PCDH10")

# Boxplot of the chosen gene among the 3 tissues considering all replicates
boxplot(assays(rse_brain)$TPM[39995,],assays(rse_colon)$TPM[39995,], 
        assays(rse_muscle)$TPM[39995,], outline=F )
