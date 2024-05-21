library(recount3)
library(recount)
library(edgeR)

rse_brain <- readRDS("rse_brain.RDS")
rse_colon <- readRDS("rse_colon.RDS")
rse_muscle <- readRDS("rse_muscle.RDS")

# Those are brain samples: 2931 columns (= different individuals) and 54042 rows (= genes)
# rowData: I have 25 additional annotations;
# colData: 198 additional piece of information for each sample.
rse_brain

# The counts are in "coverage" format
assays(rse_brain)$counts <- transform_counts(rse_brain)

# I can see wheter or not I have enough information to convert the counts I find in table
# to obtain CPM (I need library size and gene length informations)

names(colData(rse_brain))
# Tutte le infor aggiuntive per ogni sample; tutto ciò che inizia per gtex
# è related to the sample; ho il sesso e età del donor; 
# recount ho tutti i risultati coming from their processes;
# 77 e 78 %di reads mapping on chromosome x and Y
# qc.star --> all the parameters and results obtained by the output using STAR for mapping


table(colData(rse_brain)$gtex.smrin)
# Il RIN value di ogni sample (numero da 0 a 10)

table(colData(rse_brain)$gtex.smtsd)
#  mi dice la specifica portion del brain each of the sample came from.
# sono diversi individui e diverse parti del cervello, quindi la variability è mlto alta

table(colData(rse_brain)$gtex.sex) # Biologicaylly dovrebbe essere con o senza chrY nei dati

table(colData(rse_brain)$gtex.age)
# Devo mettere range per privacy issue

head(colData(rse_brain)$"recount_qc.star.number_of_input_reads_both")
# Number of reads before the mapping (paired end reads)

head(colData(rse_brain)$'recount_qc.star.uniquely_mapped_reads_%_both')
# and the percentage of uniquely mapped reads - the ones employed for computing the read counts:

inputreads <- colData(rse_brain)$"recount_qc.star.number_of_input_reads_both"
boxplot(inputreads)

minlib <- min(inputreads)
maxlib <- max(inputreads)
minlib
maxlib

mapped <- colData(rse_brain)$"recount_qc.star.uniquely_mapped_reads_%_both"
boxplot(mapped)
# distribuzione di % mapped reads on the genome; il boxplot è way over 80% 
# ma ci sono molti outliers, alcuni sample hanno % di mapping reads bassissime
# il minimo è 6.6%
min(mapped)

min_map <- which.min(mapped)
colData(rse_brain)[min_map,]


rinplot <- colData(rse_brain)$gtex.smrin
boxplot(rinplot)
# Median is on 7, few outliers very good, ma some at 5 --> controllo RIN prima di seuencing
# quindi loro sono stati molto permissivi perché il più basso è 5

#Questi erano total RNA seq: quindi guardiamo esimtated % of rRNA in my samples
head(colData(rse_brain)$gtex.smrrnart) # il consiglio è rimane sotto il 10% di rRNA
boxplot(colData(rse_brain)$gtex.smrrnart)

head(colData(rse_brain)$"recount_qc.aligned_reads%.chrm") 
boxplot(colData(rse_brain)$"recount_qc.aligned_reads%.chrm") # samples in cui 
# 70% of reads were due only to mtDNA



plot(colData(rse_brain)$gtex.smrrnart, colData(rse_brain)$"recount_qc.aligned_reads%.chrm")
# % of rRNA va % of mtDNA, ogni punto è un sample
# --> perfect correlation: più ribosomal ho, più mitocondrial
# --> why? 

names(rowData(rse_brain)) # Info for all genes; ho name, gene_name, gene_id, gene_id, gene_synonym...
# they put to each gene annotation to what the gene is:
# produce mRNA, o rRNA, o tRNA, o 'Gene' whch means pseudo-gene, 
# C_Region, D_segment, V_ssegment --> loci where there is individual variability (immunoglobulin)
# "Gene" is indeed a pseudogene --> genomic regions that look like a gene
# On the human genome there are some regions that look like protein coding regions but without
# any split in exons and introns --> even if they look like protein coding genes,
# at time they didnt find trace of gene activity. 2 sources: mature RNA retrotranscribed 
# and insert back into the genome; or gene duplication e una delle due copie fu inattivata.
# dato che ora sequenzimao molto, le persone hanno visto che gli pseudo gene si comportano da 
# geni, a low levels they find evidence of transcruption --> quindi invcee d iannotarli 
# come novel genes, they kept them as a separate category.

table(rowData(rse_brain)$gbkey)

rowRanges(rse_brain)
# bp_length is the gene lenght

table(rowRanges(rse_brain)@seqnames)
# Ho i chanonical chromosomes; sotto ho tutti glialternative chromosomes e via dicendo
# quindi ho mapping su anche i non-chanonical chromosomes

mito <- rse_brain[rowRanges(rse_brain)@seqnames == 'chrM']
dim(mito) # 

mito
rowRanges(mito)

