library(recount3)
library(recount)
library(edgeR)

rse_brain <- readRDS("rse_brain.RDS")

# Those are brain samples: 2931 columns (= different individuals) and 54042 rows (= genes).
# rowData: I have 25 additional annotations.
# colData: I have 198 additional pieces of information for each sample.
rse_brain

# The counts are in "coverage" format.
assays(rse_brain)$counts <- transform_counts(rse_brain)

# I can see whether or not I have enough information to convert the counts I find 
# in the table to obtain the CPM (I need the library size and the gene length)

# Here I see all the additional informations for each sample: sex, age,
# the one starting with 'recount' contain the results of their processes
# (of the people working for this package), % of reads mapping on chrX and chrY,
# 'qc.star' contain all the parameters and results obtained by the output of STAR
# used for the mapping
names(colData(rse_brain))

# RIN value of each sample (number that goes from 0 to 10)
table(colData(rse_brain)$gtex.smrin)

# It tells me the specific portion of the brain each of the sample comes from.
# So here we have different individuals and different portions, so very high variability.
table(colData(rse_brain)$gtex.smtsd)

table(colData(rse_brain)$gtex.sex) # Look at the presence/absence of chrY in data 

table(colData(rse_brain)$gtex.age) # They are written in range because of privacy issues.

# Number of reads before the mapping (paired-end reads)
head(colData(rse_brain)$"recount_qc.star.number_of_input_reads_both")

# % of uniquely mapped reads --> the ones employed for computing the read counts
head(colData(rse_brain)$'recount_qc.star.uniquely_mapped_reads_%_both')

inputreads <- colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"
boxplot(inputreads)

minlib <- min(inputreads)
maxlib <- max(inputreads)
minlib
maxlib

mapped <- colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped) # Distribution of % of mapped reads on the genome: the boxplot is
                # way over 80% but there are a lot of outliers which have a very
                # low mapping (the lowest is 6.6%)
min(mapped)

min_map <- which.min(mapped)
colData(rse_brain)[min_map,]


rinplot <- colData(rse_brain)$gtex.smrin
boxplot(rinplot) # Median is good (7) and there are few outliers, but some are at 5
                 # --> the RIN is checked before the sequencing, so the authors were
# very permissive in sequencing also the samples with RIN = 5

# These were total RNA-Seq: let's look at the estimated % of rRNA in my samples
head(colData(rse_brain)$gtex.smrrnart) # Under 10% is good
boxplot(colData(rse_brain)$gtex.smrrnart)

head(colData(rse_brain)$"recount_qc.aligned_reads%.chrm") 
boxplot(colData(rse_brain)$"recount_qc.aligned_reads%.chrm") # There are samples with
                                                             # 70% of mtDNA

plot(colData(rse_brain)$gtex.smrrnart, colData(rse_brain)$"recount_qc.aligned_reads%.chrm")
# % of rRNA va % of mtDNA, each point is a sample.
# --> perfect correlation: more ribosomal I have, more mitocondrial I have
# --> why? 

names(rowData(rse_brain)) # All info we have for each row (so for each gene): name,
                          # gene_name, gene_id...

table(rowData(rse_brain)$gbkey) # This tells if the gene produces mRNA, it's a
                                # pseudogene (called "Gene", which is weird), 
                                # rRNA, tRNA...
# We also see C_region, D_segment and so on --> those are the loci where there 
# is the most individual variability (immunoglobulin)
# "Gene" means that it's a pseudogene: pseudogenes are genomic regions that 
# look like a protein coding gene but without any split in exons and introns 
# --> even if they look like protein coding genes, at the time they didn't find 
# any trace of gene activity. 
# Two sources of pseudogenes: 
# 1) Mature RNA retrotranscribedand insert back into the genome
# 2) Gene duplication and one of the two copy was inactivated.
# Now we sequence a lot (and very in depth) so we start to see that also pseudogenes
# are transcribed (at low levels) --> instead of annotate them as novel genes, 
# they kept them as a separate category.

rowRanges(rse_brain) # bp_length is the gene lenght

table(rowRanges(rse_brain)@seqnames) # The mapping was done on canonical 
                                     # chromosomes, but also on the non-canonical 
                                     # chromosomes

# Several checks on the genes mapping on chrM
mito <- rse_brain[rowRanges(rse_brain)@seqnames == 'chrM']
dim(mito)
mito
rowRanges(mito)

