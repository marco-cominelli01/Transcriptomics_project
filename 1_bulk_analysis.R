library(recount3)
library(recount)
library(edgeR)


# Save the 3 files inside a RSE data object
rse_brain <- readRDS("rse_brain.RDS")
rse_colon <- readRDS("rse_colon.RDS")
rse_muscle <- readRDS("rse_muscle.RDS")


original_rse_brain <- rse_brain
original_rse_colon <- rse_colon
original_rse_muscle <- rse_muscle


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

##############
##############
##############

# Challenge 1

## BRAIN
chrM_mask <- rowRanges(rse_brain)@seqnames != "chrM"

# under the assumption that null biotypes would not be pseudogenes
pseudogene_mask <- rowRanges(rse_brain)$gene_biotype != "pseudogene" | is.na(rowRanges(rse_brain)$gene_biotype )

# under the assumption that null gbkeys would not be rRNA
rRNA_mask <- rowRanges(rse_brain)$gbkey != "rRNA" | is.na(rowRanges(rse_brain)$gbkey)

ncRNA_mask <- rowRanges(rse_brain)$bp_length >= 200

filter <- chrM_mask & pseudogene_mask & rRNA_mask & ncRNA_mask
filtered_rse_brain <- rse_brain[filter, ] # 40903 rows out of 54042

colSums(assays(rse_brain)$counts[ , c(33, 34, 38)]) # 28510344 30433721 27988026
colSums(assays(rse_brain[chrM_mask, ])$counts[, c(33, 34, 38)]) # remaining rows after removing only chrM genes: 27513103 28788076 25370943
colSums(assays(rse_brain[pseudogene_mask, ])$counts[, c(33, 34, 38)]) # remaining rows after removing only pseudogenes: 28004851 29944047 27471401
colSums(assays(rse_brain[rRNA_mask, ])$counts[, c(33, 34, 38)]) # remaining rows after removing only rRNA genes: 27316214 28695874 25385419
colSums(assays(rse_brain[ncRNA_mask, ])$counts[, c(33, 34, 38)]) # remaining rows after removing only short ncRNA: 28457517 30367047 27829919
colSums(assays(filtered_rse_brain)$counts[, c(33, 34, 38)]) # all of the four removed: 26764344 28140060 24714616

### COLON

chrM_mask <- rowRanges(rse_colon)@seqnames != "chrM"

# under the assumption that null biotypes would not be pseudogenes
pseudogene_mask <- rowRanges(rse_colon)$gene_biotype != "pseudogene" | is.na(rowRanges(rse_colon)$gene_biotype )

# under the assumption that null gbkeys would not be rRNA
rRNA_mask <- rowRanges(rse_colon)$gbkey != "rRNA" | is.na(rowRanges(rse_colon)$gbkey)

ncRNA_mask <- rowRanges(rse_colon)$bp_length >= 200

filter <- chrM_mask & pseudogene_mask & rRNA_mask & ncRNA_mask
filtered_rse_colon <- rse_colon[filter, ] # 40903 rows out of 54042

colSums(assays(rse_colon)$counts[ , c(32, 33, 34)]) # 28510344 30433721 27988026
colSums(assays(rse_colon[chrM_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only chrM genes: 27513103 28788076 25370943
colSums(assays(rse_colon[pseudogene_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only pseudogenes: 28004851 29944047 27471401
colSums(assays(rse_colon[rRNA_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only rRNA genes: 27316214 28695874 25385419
colSums(assays(rse_colon[ncRNA_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only short ncRNA: 28457517 30367047 27829919
colSums(assays(filtered_rse_colon)$counts[, c(32, 33, 34)]) # all of the four removed: 26764344 28140060 24714616



### MUSCLE

chrM_mask <- rowRanges(rse_muscle)@seqnames != "chrM"

# under the assumption that null biotypes would not be pseudogenes
pseudogene_mask <- rowRanges(rse_muscle)$gene_biotype != "pseudogene" | is.na(rowRanges(rse_muscle)$gene_biotype )

# under the assumption that null gbkeys would not be rRNA
rRNA_mask <- rowRanges(rse_muscle)$gbkey != "rRNA" | is.na(rowRanges(rse_muscle)$gbkey)

ncRNA_mask <- rowRanges(rse_muscle)$bp_length >= 200

filter <- chrM_mask & pseudogene_mask & rRNA_mask & ncRNA_mask
filtered_rse_muscle <- rse_muscle[filter, ] # 40903 rows out of 54042

colSums(assays(rse_muscle)$counts[ , c(32, 33, 34)]) # 28510344 30433721 27988026
colSums(assays(rse_muscle[chrM_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only chrM genes: 27513103 28788076 25370943
colSums(assays(rse_muscle[pseudogene_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only pseudogenes: 28004851 29944047 27471401
colSums(assays(rse_muscle[rRNA_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only rRNA genes: 27316214 28695874 25385419
colSums(assays(rse_muscle[ncRNA_mask, ])$counts[, c(32, 33, 34)]) # remaining rows after removing only short ncRNA: 28457517 30367047 27829919
colSums(assays(filtered_rse_muscle)$counts[, c(32, 33, 34)]) # all of the four removed: 26764344 28140060 24714616


rse_brain <- filtered_rse_brain
rse_colon <- filtered_rse_colon
rse_muscle <- filtered_rse_muscle

rm(filtered_rse_brain, filtered_rse_colon, filtered_rse_muscle)

##### Challenge 2

## BRAIN

# Duplicates - or unique_duplicated_gene_names
duplicates_names <- names(table(rowRanges(rse_brain)$Name)[table(rowRanges(rse_brain)$Name) > 1]) # group by name
length(duplicates_names) # 737

canonchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
              "chrX", "chrY", "chrM")

duplicated_genes_index <- rowRanges(rse_brain)$Name %in% duplicates_names
duplicated_entries <- rse_brain[duplicated_genes_index, ] # 2415 entries corresponding to genes that have duplicate $Name

# Some of the entries that correspond to duplicated genes fall into cannonchr and some not
duplicated_canonical_entries <- duplicated_entries[rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 671 entries
duplicated_non_canonical_entries <- duplicated_entries[!rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 1744 entries

# unique Names
duplicated_canonical_names <- names(table(rowRanges(duplicated_canonical_entries)$Name)) # 1167 unique genes, thus: 1180 - 1167 = 13 repeats that fall on other canonicals;
table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1] # shows that there are actually 13 genes (those 13 repeats are evenly distributed among 13 genes), each with 1 other repeat falling in canonicals
duplicated_non_canonical_names <- names(table(rowRanges(duplicated_non_canonical_entries)$Name)) # 1228 unique genes; way more "non canonical repeats" falling in other non canonicals

# 1167 + 1228 - 1240 = 1155 genes that are "repeated" from canonical to non canonical
# Observation: since there were 1167 unique genes with repeats in the "canonical group" and 13 of them had also repeats falling in other canonicals, knowing
# that 1155 must have a repeat both in canoncial and non canonical, that leaves 1167 - 1155 = 12 "canonicals" NOT having repeats falling in non canonicals.
# Which are those genes? for sure not the ones with repeat number of 1 in canonicals; if they have made it so far, they for sure are duplicates, so if they are "repeated"
# only once inside canonicals, there must be at least one repeat falling in non canonicals. They must be the ones from the group of 13 that had repeats falling
# in other canonicals. Therefore, only one of them is repeated both in other canoncials as well as non canonicals, while the 12, are repeated only on other
# canonicals.
table(duplicated_non_canonical_names %in% duplicated_canonical_names) # 1155

# Indeed:
canonical_12 <- duplicated_canonical_names[(!(duplicated_canonical_names %in% duplicated_non_canonical_names))] # the 12 that fall ONLY in other canonicals
group_of13 <- names(table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1]) # the 13 that fall in canonicals
# from the original group of 13 that fall in other canonicals, only one falls also in non canonicals
# namely:
group_of13[!(group_of13 %in% canonical_12)] # FABP5P13 is repeated both in canonicals as well as in non canonicals

# Rse with canonicals only
canonical_rse_brain <- rse_brain[rowRanges(rse_brain)@seqnames %in% canonchr, ]

## COLON

# Duplicates - or unique_duplicated_gene_names
duplicates_names <- names(table(rowRanges(rse_colon)$Name)[table(rowRanges(rse_colon)$Name) > 1]) # group by name
length(duplicates_names) # 737

canonchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
              "chrX", "chrY", "chrM")

duplicated_genes_index <- rowRanges(rse_colon)$Name %in% duplicates_names
duplicated_entries <- rse_colon[duplicated_genes_index, ] # 2415 entries corresponding to genes that have duplicate $Name

# Some of the entries that correspond to duplicated genes fall into cannonchr and some not
duplicated_canonical_entries <- duplicated_entries[rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 671 entries
duplicated_non_canonical_entries <- duplicated_entries[!rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 1744 entries

# unique Names
duplicated_canonical_names <- names(table(rowRanges(duplicated_canonical_entries)$Name)) # 1167 unique genes, thus: 1180 - 1167 = 13 repeats that fall on other canonicals;
table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1] # shows that there are actually 13 genes (those 13 repeats are evenly distributed among 13 genes), each with 1 other repeat falling in canonicals
duplicated_non_canonical_names <- names(table(rowRanges(duplicated_non_canonical_entries)$Name)) # 1228 unique genes; way more "non canonical repeats" falling in other non canonicals

# 1167 + 1228 - 1240 = 1155 genes that are "repeated" from canonical to non canonical
# Observation: since there were 1167 unique genes with repeats in the "canonical group" and 13 of them had also repeats falling in other canonicals, knowing
# that 1155 must have a repeat both in canoncial and non canonical, that leaves 1167 - 1155 = 12 "canonicals" NOT having repeats falling in non canonicals.
# Which are those genes? for sure not the ones with repeat number of 1 in canonicals; if they have made it so far, they for sure are duplicates, so if they are "repeated"
# only once inside canonicals, there must be at least one repeat falling in non canonicals. They must be the ones from the group of 13 that had repeats falling
# in other canonicals. Therefore, only one of them is repeated both in other canoncials as well as non canonicals, while the 12, are repeated only on other
# canonicals.
table(duplicated_non_canonical_names %in% duplicated_canonical_names) # 1155

# Indeed:
canonical_12 <- duplicated_canonical_names[(!(duplicated_canonical_names %in% duplicated_non_canonical_names))] # the 12 that fall ONLY in other canonicals
group_of13 <- names(table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1]) # the 13 that fall in canonicals
# from the original group of 13 that fall in other canonicals, only one falls also in non canonicals
# namely:
group_of13[!(group_of13 %in% canonical_12)] # FABP5P13 is repeated both in canonicals as well as in non canonicals

# Rse with canonicals only
canonical_rse_colon <- rse_colon[rowRanges(rse_colon)@seqnames %in% canonchr, ]

## MUSCLE

# Duplicates - or unique_duplicated_gene_names
duplicates_names <- names(table(rowRanges(rse_muscle)$Name)[table(rowRanges(rse_muscle)$Name) > 1]) # group by name
length(duplicates_names) # 737

canonchr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
              "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15",
              "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22",
              "chrX", "chrY", "chrM")

duplicated_genes_index <- rowRanges(rse_muscle)$Name %in% duplicates_names
duplicated_entries <- rse_muscle[duplicated_genes_index, ] # 2415 entries corresponding to genes that have duplicate $Name

# Some of the entries that correspond to duplicated genes fall into cannonchr and some not
duplicated_canonical_entries <- duplicated_entries[rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 671 entries
duplicated_non_canonical_entries <- duplicated_entries[!rowRanges(duplicated_entries)@seqnames %in% canonchr, ] # 1744 entries

# unique Names
duplicated_canonical_names <- names(table(rowRanges(duplicated_canonical_entries)$Name)) # 1167 unique genes, thus: 1180 - 1167 = 13 repeats that fall on other canonicals;
table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1] # shows that there are actually 13 genes (those 13 repeats are evenly distributed among 13 genes), each with 1 other repeat falling in canonicals
duplicated_non_canonical_names <- names(table(rowRanges(duplicated_non_canonical_entries)$Name)) # 1228 unique genes; way more "non canonical repeats" falling in other non canonicals

# 1167 + 1228 - 1240 = 1155 genes that are "repeated" from canonical to non canonical
# Observation: since there were 1167 unique genes with repeats in the "canonical group" and 13 of them had also repeats falling in other canonicals, knowing
# that 1155 must have a repeat both in canoncial and non canonical, that leaves 1167 - 1155 = 12 "canonicals" NOT having repeats falling in non canonicals.
# Which are those genes? for sure not the ones with repeat number of 1 in canonicals; if they have made it so far, they for sure are duplicates, so if they are "repeated"
# only once inside canonicals, there must be at least one repeat falling in non canonicals. They must be the ones from the group of 13 that had repeats falling
# in other canonicals. Therefore, only one of them is repeated both in other canoncials as well as non canonicals, while the 12, are repeated only on other
# canonicals.
table(duplicated_non_canonical_names %in% duplicated_canonical_names) # 1155

# Indeed:
canonical_12 <- duplicated_canonical_names[(!(duplicated_canonical_names %in% duplicated_non_canonical_names))] # the 12 that fall ONLY in other canonicals
group_of13 <- names(table(rowRanges(duplicated_canonical_entries)$Name)[table(rowRanges(duplicated_canonical_entries)$Name) > 1]) # the 13 that fall in canonicals
# from the original group of 13 that fall in other canonicals, only one falls also in non canonicals
# namely:
group_of13[!(group_of13 %in% canonical_12)] # FABP5P13 is repeated both in canonicals as well as in non canonicals


# Rse with canonicals only
canonical_rse_muscle <- rse_muscle[rowRanges(rse_muscle)@seqnames %in% canonchr, ]

##############
##############
##############

rse_brain <- canonical_rse_brain
rse_colon <- canonical_rse_colon
rse_muscle <- canonical_rse_muscle


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
qlfBC <- glmQLFTest(fit, contrast=c(1,-1,0))

# Colon (top) vs Muscle (bottom)
qlfCM <- glmQLFTest(fit, contrast=c(0,1,-1))

# Muscle (top) vs Brain (bottom)
qlfMB <- glmQLFTest(fit, contrast=c(-1,0,1))
qlfBM <- glmQLFTest(fit, contrast=c(1,0,-1))

# Muscle (top) vs Colon (bottom)
qlfMC <- glmQLFTest(fit, contrast=c(0,-1,1))

resultsCB <- topTags(qlfCB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsCB, "resultsCB.txt")

resultsMB <- topTags(qlfMB, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMB, "resultsMB.txt")

resultsMC <- topTags(qlfMC, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)
write.table(resultsMC, "resultsMC.txt")

resultsCM <- topTags(qlfCM, n = 10000000, adjust.method = "BH", sort.by = "PValue", p.value = 1)


# Summary for the up- and down- regulated genes for the chosen thresholds of p-value and logFC
summary(decideTests(qlfCB, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfMB, p.value=0.05, adjust.method = "BH", lfc=0))
summary(decideTests(qlfMC, p.value=0.05, adjust.method = "BH", lfc=0))


## Colon

genes_up_colon_vs_brain <- decideTests(qlfCB, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_colon_vs_brain <- rownames(genes_up_colon_vs_brain)[genes_up_colon_vs_brain == 1]

genes_up_colon_vs_muscle <- decideTests(qlfCM, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_colon_vs_muscle <- rownames(genes_up_colon_vs_muscle)[genes_up_colon_vs_muscle == 1]

up_colon_both <- intersect(up_colon_vs_brain, up_colon_vs_muscle)

### Muscle

genes_up_muscle_vs_colon <- decideTests(qlfMC, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_muscle_vs_colon <- rownames(genes_up_muscle_vs_colon)[genes_up_muscle_vs_colon == 1]

genes_up_muscle_vs_brain <- decideTests(qlfMB, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_muscle_vs_brain <- rownames(genes_up_muscle_vs_brain)[genes_up_muscle_vs_brain == 1]

up_muscle_both <- intersect(up_muscle_vs_colon, up_muscle_vs_brain)

### Brain

genes_up_brain_vs_colon <- decideTests(qlfBC, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_brain_vs_colon <- rownames(genes_up_brain_vs_colon)[genes_up_brain_vs_colon == 1]

genes_up_brain_vs_muscle <- decideTests(qlfBM, p.value=0.05, adjust.method = "BH", lfc=0)@.Data
up_brain_vs_muscle <- rownames(genes_up_brain_vs_muscle)[genes_up_brain_vs_muscle == 1]

up_brain_both <- intersect(up_brain_vs_colon, up_brain_vs_muscle)





#### Scelta di un gene over espresso in muscolo rispetto agli altri due
filtered_results <- resultsMB[rownames(resultsMB) %in% up_muscle_both, ]

ordered_results <- filtered_results[order(filtered_results$table$logFC, decreasing = TRUE), ]

filtered_results_muscle_vs_colon <- resultsMC[rownames(resultsMC) %in% up_muscle_both, ]
filtered_results_muscle_vs_colon[rownames(filtered_results_muscle_vs_colon) == 'CLCN1',]

### Colon
filtered_results <- resultsCB[rownames(resultsCB) %in% up_colon_both, ]
ordered_results <- filtered_results[order(filtered_results$table$logFC, decreasing = TRUE), ]
top_200_genes_CB <- rownames(ordered_results[1:200, ])

filtered_results <- resultsCM[rownames(resultsCM) %in% up_colon_both, ]
ordered_results <- filtered_results[order(filtered_results$table$logFC, decreasing = TRUE), ]
top_200_genes_CM <- rownames(ordered_results[1:200, ])

write.table(top_200_genes_CB, "top200_CB.txt", row.names = F, quote = F)
write.table(top_200_genes_CM, "top200_CM.txt", row.names = F, quote = F)




# !!! Update R to the latest version to execute the commands below
# (or compute the TMP values by hand)
assays(rse_brain)$TPM <- recount::getTPM(rse_brain)
assays(rse_colon)$TPM <- recount::getTPM(rse_colon)
assays(rse_muscle)$TPM <- recount::getTPM(rse_muscle)
which(rowData(rse_brain)$gene_name == "CLCN1") # 31616

# Boxplot of the chosen gene among the 3 tissues considering all replicates
boxplot(assays(rse_brain)$TPM[31616,],assays(rse_colon)$TPM[31616,], 
        assays(rse_muscle)$TPM[31616,], outline=F )
