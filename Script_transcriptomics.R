### For bulk RNA-Seq
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("recount3")

BiocManager::install(c(
  "recount", "GenomicRanges", "limma", "edgeR", "DESeq2",
  "regionReport", "clusterProfiler", "org.Hs.eg.db", "gplots",
  "derfinder", "GenomicState", "bumphunter", "derfinderPlot", "sessioninfo"
))

### For single cell RNA_Seq
install.packages("tidyverse")
install.packages('Seurat')
library(Seurat)  # In case of Warning message, enter 'y'

install.packages("devtools")
devtools::install_github("thomasp85/patchwork")

