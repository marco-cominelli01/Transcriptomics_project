
###################

library(recount3)
library(recount)
library(edgeR)


colData(rse_brain)$gtex.smrin[38]
colData(rse_brain)$gtex.smrin[33]
colData(rse_brain)$gtex.smrin[34]

colData(rse_colon)$gtex.smrin[32]
colData(rse_colon)$gtex.smrin[33]
colData(rse_colon)$gtex.smrin[34]

colData(rse_muscle)$gtex.smrin[32]
colData(rse_muscle)$gtex.smrin[33]
colData(rse_muscle)$gtex.smrin[34]



colData(rse_brain)$gtex.smrrnart[38]
colData(rse_brain)$gtex.smrrnart[33]
colData(rse_brain)$gtex.smrrnart[34]

colData(rse_colon)$gtex.smrrnart[32]
colData(rse_colon)$gtex.smrrnart[33]
colData(rse_colon)$gtex.smrrnart[34]

colData(rse_muscle)$gtex.smrrnart[32]
colData(rse_muscle)$gtex.smrrnart[33]
colData(rse_muscle)$gtex.smrrnart[34]



colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[38]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]

colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[32]
colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_colon)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]

colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[32]
colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[33]
colData(rse_muscle)$"recount_qc.star.uniquely_mapped_reads_%_both"[34]


assays(rse_brain)$counts <- transform_counts(rse_brain)
assays(rse_colon)$counts <- transform_counts(rse_colon)
assays(rse_muscle)$counts <- transform_counts(rse_muscle)

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
# this can be done only once on the final table.
rownames(x) <- rowData(rse_brain_selected)$gene_name

y <- DGEList(counts=x)


group <- as.factor(c("Brain","Brain","Brain","Colon","Colon","Colon",
                     "Muscle","Muscle","Muscle"))
y$samples$group <- group


y$samples$rin <- as.factor(c(colData(rse_brain_selected)$gtex.smrin,colData(rse_colon_selected)$gtex.smrin,colData(rse_muscle_selected)$gtex.smrin))

y$samples$slice <- as.factor(c(colData(rse_brain_selected)$gtex.smtsd,colData(rse_colon_selected)$gtex.smtsd,colData(rse_muscle_selected)$gtex.smtsd))

y$samples$sex <- as.factor(c(colData(rse_brain_selected)$gtex.sex,colData(rse_colon_selected)$gtex.sex,colData(rse_muscle_selected)$gtex.sex))

y$samples$age <- as.factor(c(colData(rse_brain_selected)$gtex.age,colData(rse_colon_selected)$gtex.age,colData(rse_muscle_selected)$gtex.age))

y$samples$rRNA <- as.factor(c(colData(rse_brain_selected)$gtex.smrrnart,colData(rse_colon_selected)$gtex.smrrnart,colData(rse_muscle_selected)$gtex.smrrnart))

y$samples$mapped <- as.factor(c(colData(rse_brain_selected)$"recount_qc.star.uniquely_mapped_reads_%_both", colData(rse_colon_selected)$"recount_qc.star.uniquely_mapped_reads_%_both",colData(rse_muscle_selected)$"recount_qc.star.uniquely_mapped_reads_%_both"))

y$samples$chrm <- as.factor(c(colData(rse_brain_selected)$"recount_qc.aligned_reads%.chrm", colData(rse_colon_selected)$"recount_qc.aligned_reads%.chrm",colData(rse_muscle_selected)$"recount_qc.aligned_reads%.chrm"))

recount_qc.star.%_reads_unmapped:_too_short_both

y
# Library size: at least 30M/50M reads; perchè in alcni casi ne ho solo 10? 
# il problema è nel transform_counts --> ignoro quindi la colonna lib.size
# non è la real library size, è side effect della conversion they made to give me the read counts
# è un artefatto, non è la real library size of each sample

table(rowSums(y$counts==0)==9)

keep.exprs <- filterByExpr(y, group=group)
y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

logcpm_before <- cpm(y, log=TRUE)
y <- calcNormFactors(y, method = "TMM")
y

logcpm_after <- cpm(y, log=TRUE)
par(mfrow = c(1,2))
boxplot(logcpm_before, notch=T)
boxplot(logcpm_after, notch=T)

design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

# Colon: potrebbe andare closer to another tissue instead of staying close to other replicates of colon
# Unu replicate di Brain è un po' lontano dagli altri: 
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)

# Guardo queste cose per capire perch questo sample è un po' diverso dagli altri
# Aggiungi altre colonne per vedere di chi è la colpa
plotMDS(logcpm, labels=y$samples$rRNA)
plotMDS(logcpm, labels=y$samples$chrm)
plotMDS(logcpm, labels=y$samples$age)
plotMDS(logcpm, labels=y$samples$mapped)
plotMDS(logcpm, labels=y$samples$rin)
plotMDS(logcpm, labels=y$samples$slice)
plotMDS(logcpm, labels=y$samples$sex)
plotMDS(logcpm, labels=rownames(y$samples)) # To have the names of the samples displayed


# le threshold più restrittive sono 0.01 for FDR e logFC almeno = 1 
# se penso che i geni DE siano troppi, li ordino per increasing p-value
# e poi considero i primi 100, che geni trovo? Poi guardo ai primi 1000, che geni trovo?
# quando controllo per sequencing lane devo avere 0 geni DE
# 

# Always save txt files.

y <- estimateDisp(y, design)
plotBCV(y) # We don't need to comment this plot (just take it as it is)


fit <- glmQLFit(y, design)

#colon (top) vs brain (bottom)
qlfCB <- glmQLFTest(fit, contrast=c(-1,1,0))
#muscle (top) vs brain (bottom)
qlfMB <- glmQLFTest(fit, contrast=c(-1,0,1))
#muscle (top) vs colon (bottom)
qlfMC <- glmQLFTest(fit, contrast=c(0,-1,1))

resultsCB <- topTags(qlfCB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsCB, "resultsCB.txt")

resultsMB <- topTags(qlfMB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMB, "resultsMB.txt")

resultsMC <- topTags(qlfMC, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMC, "resultsMC.txt")

# Summary for the up- and down- regulated genes for the chosen thresholds of p-value and logFC
summary(decideTests(qlfCB, p.value=0.05, adjust.method = "BH", lfc=0))


####### AGGIORNARE R O CALCOLARE I TPM A MANO
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_heart)$TPM <- recount::getTPM(rse_heart)
assays(rse_pancreas)$TPM <- recount::getTPM(rse_pancreas)
which(rowData(rse_brain)$gene_name == "PCDH10")

####### BOXPLOT OF THE CHOSEN GENE AMONG THE 3 TISSUES CONSIDERING ALL THE REPLICATES
boxplot(assays(rse_brain)$TPM[39995,],assays(rse_heart)$TPM[39995,], 
        assays(rse_pancreas)$TPM[39995,], outline=F )
