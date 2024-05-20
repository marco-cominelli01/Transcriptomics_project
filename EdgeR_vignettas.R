library(edgeR)

x <- read.table("mouse.exercise.txt",header=T,row.names=1)
head(x)

size <- colSums(x)   # Library size of each sample
size

y <- DGEList(counts=x) # Special data object that the package use to fill it 
                       # with all the data and results of the analysis.
                       # Now it contains only the count table and some information about the samples. 
                       # 9 samples, their names are the one found in the count table, then it computes 
                       # the library size by itself and it adds two additional vectors, 
                       # group and normalization factor, by default initialized with only 1; 
                       # 'group' is for defining which are the replicates of what in our conditions, 
                       # normalization factors are the normalization factors that will be computed when we'll 
                       # perform the normalization of the original raw counts
                       # The starting count table contains the raw counts.
y

# Here we are in the scenario of a starting cell line that differentiates into 2 others cell lines.
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", "Basal", "ML", "LP"))
y$samples$group <- group
y

# Sequencing lanes: the samples were loaded on the sequencer in the same run but splitted in different
# sequencing lanes.
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
y$samples$lane <- lane
y

head(y$counts)

# The following are preparation steps:
table(rowSums(y$counts==0)==9)  # Those are the genes with 0 counts in all the samples.
                                # The statistical estimation of the parameteres uses all the genes; 
                                # if I consider also the genes with mean = 0 and variance = 0, the 
                                # estimation will be confused.


# Default function employed by edgeR: it removes from the table the genes which have at 
# least one '0' in each condition
# Example1: if a gene have 0 in condition1 and 10 in condition2 I keep it
# Example2: if a gene have values (for replicates) in condition1 of 0,10,5 and 0,20,7 in condition2, 
# I remove it.
keep.exprs <- filterByExpr(y, group=group) 


y <- y[keep.exprs,, keep.lib.sizes=FALSE]
dim(y)

# It converts the raw counts into CPM --> 'log = TRUE' means that in the table it will not put
# CPM but log2(CPM)
logcpm_before <- cpm(y, log=TRUE)

# TMM normalization: in this function I have also other way to do normalization.
y <- calcNormFactors(y, method = "TMM")
y

# Now in fact the norm.factors are different from 1.
# We apply this normalization to obtain the situation in which most of the genes 
# don't change their expression.
# These norm.factors are the values used to multiply the original raw counts values.

# When I have the normalization factors, everything from that moment on will be done
# on normalized values, not on raw counts.

# So now I compute log(CPM), but after the normalization
logcpm <- cpm(y, log=TRUE)

# After the normalization, the medians are almost the same (each boxplot is shifted up or down).
# The effect of this shifting is that samples tend to have more or less the same median value
# (here the shift is not so much visible).
par(mfrow=c(1,2))
boxplot(logcpm)
boxplot(logcpm_before)

# The next step is to tell what are the dependencies among the different conditions 
# and replicates in my experiment.
design <- model.matrix(~group, data=y$samples)

# If we do not want to use an intercept, the syntax is model.matrix(~0+group, data=y$samples)
colnames(design) <- levels(y$samples$group)
design

# LP is beta1, ML is beta2, Basal is beta0.
# The model is: beta0 + beta1*X1 + beta2*X2
# Why the function choose by itself 'basal' as baseline? Because it's the first in alphabetical order,
# and also the others two are placed in alphabetical order.

# We don't model an explicit relationship between LP and ML.

# Quality control
logcpm <- cpm(y, log=TRUE)
plotMDS(logcpm, labels=group)

# Let's suppose that each sample is a point in the space and its coordinates are 
# the expression values of each gene --> so each sample is a point in a 20.000
# dimensional space --> the plot is a projection into 2D space of all the samples
# (human readable).
# The more similar two samples are according to expression of their genes, 
# the closer the points will be in 2D space. 
# 9 samples, so 9 points: the closer 2 samples are from the expression of genes
# point of view, the closest they are here in 2D.
# It's a quality control because I see that replicates are indeed of the same condition,
# they are similar.
# If I see the points mixed, I might have some problem in the data (in these cases
# we go back to check fastQC, % of rRNA, RIN, % of mapped reads...).
# plotMDS stands for 'multidimensional scaling' --> I have a dataset 
# with too many dimensions, I reduce it to 2 dimensions.

# I look if there's some influence of the sequencing lane
plotMDS(logcpm, labels=lane)
# The plot is good because they are all mixed up; otherwise I might have some bias
# in the sequencing.

# Estimate of dispersion (gene by gene the dispersion of expression of each of the gene across all the conditions)
# X-axis: gene by gene its average expression in logCPM
# Y-axis: square root of dispersion
# BCV is the square root of dispersion, key parameter for the estimation of the variance.
# Trend line used to correct original values and shrink them: the final estimate is 
# composed by the red line (common dispersion) and a gene specific contribution 
# that is computed according to the position of the original point and the trend line (blue).
y <- estimateDisp(y, design)
plotBCV(y)

y

y$common.dispersion

# Quasi-Linear fit: this function computes the beta values for all genes in all 
# conditions.
fit <- glmQLFit(y, design)

# Now everything is ready for the statistical test
# The first column in the design matrix was the baseline condition; 
# now I want to compare LP with basal: basal is the beta0
# 'coef = 2' means --> take the second column and compare it with the first one; 
# the function already knows implicitly to compare 2 with the first column (the baseline)

qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1)
# In the logFC, the condition specified (in this case 2) will be the numerator 
# and the intercept the denominato (condition 1).
# logFC can be positive or negative: in this case LP is at numerator,
# logCPM is the average logCPM across this two conditions,
# F is the test-statistic value, for which it computes the p-value and the FDR.
# Typically the genes with FDR lower than 0.01 or 0.05 are the ones selected.

FDR <- p.adjust(qlf.2vs1$table$PValue, method="BH")
sum(FDR < 0.05)

summary(decideTests(qlf.2vs1))

summary(decideTests(qlf.2vs1, p.value=0.01, lfc=1)) 
# Here I also add a threshold on the logFC

#### This is ML vs basal
qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)
# I test 3 vs 2, ML vs LP, so ignoring the base condition.
# If logFC is negative, the expression is higher in LP, 
# If positive, it means that the expression is higher in ML.

# We now write explicitly what condition has to be compared to which. 
# In the logFC, the condition "1" will be the numerator and the "-1" the denominator.
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.3vs2)

# I can also perform a statistical test taking all the conditions together.
# I specify a range of coefficients: ANOVA test, which is a generalization of 
# the t-test to more than 2 samples. 
# Considering that here we have 3 conditions, I do the test considering all 3 
# conditions. The p-value tells me which are the most variable genes across
# the conditions I'm studying (but I don't know if the change is significant
# in each condition).
qlf <- glmQLFTest(fit, coef=2:3)
topTags(qlf)

# I try to call DE genes grouping by sequencing lanes (it's a sort of quality check)
design <- model.matrix(~0+lane,data=y$samples) 
# '0' means that there isn't an intercept (the 3 lanes are independent with one another)
# In this case the number of replicates is not the same across conditions
# (3 for lane4, 4 for lane6 and 2 for lane8) --> it's not advisable but it still works.
# With only 1 replicate in a condition would throw an error (unreliable estimate 
# of the parameters) but with 2 is ok.

design

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

qlfLane1.2 <- glmQLFTest(fit, contrast=c(-1,1,0))  
topTags(qlfLane1.2) 
# With a threshold of 0.05 on the FDR, I consider only 1 gene to be DE
# (I expected 0 but this is good anyway)

qlfLane1.3 <- glmQLFTest(fit, contrast=c(-1,0,1))
topTags(qlfLane1.3)

summary(decideTests(qlfLane1.2))
# Happily, I get only 1 DE genes in all the comparisons 
# (otherwise there would have been sequencing bias, batch effect)
# Also, this gene is borderline (with a threshold of 0.01 it wouldn't be considered DE).

summary(decideTests(qlfLane1.3))


