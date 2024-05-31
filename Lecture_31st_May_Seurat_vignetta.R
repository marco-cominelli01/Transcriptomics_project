library(dplyr)
library(Seurat)
library(patchwork)

# First thing to do in any case: study the type of cells I'm dealing with,
# to know which cell types (more or less) to expect.

# Load the PBMC dataset (you need to download the zip from the link in the Vignetta)
pbmc.data <- Read10X(data.dir = "./pbmc3k_filtered_gene_bc_matrices/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized) data.
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", 
                           min.cells = 3, min.features = 200)
# Notice that the parameters filter out all cells in which there are less than 
# 200 genes expressed (min.features) and genes expressed in less that 3 cells 
# (min.cells). The matrix will contain 13714 genes and 2700 cells. That is, an 
# initial cell/gene removal. The parameters used are here usually not changed.

# The 'features' are the genes!!!
pbmc

# Lets examine the counts for a few genes in the first 30 cells. 
# They are in pbmc.data
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# It's called sparse matrix because it has a lot of empty cells.
# It's normal to find a lot of 'low' values for the counts (like 1,2,3...)

# The label put to each column is the anonymous barcode assigned.
head(colnames(pbmc))

# Select rows with name starting with MT- and then compute the 
# % of reads of each cell going to these genes.
grep("^MT-",rownames(pbmc),value = TRUE)

# (Also remember that in mouse, gene names are usually in lowercase (Mt- or even mt-)).
# The [[ operator can add columns to object metadata. This is a great place 
# to store additional info/data
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# We can play around a little with this feature. We saw how ribosomal protein genes 
# “eat up” a lot of reads because they're highly expressed. 
# Their gene symbol usually starts by RPL or RPS:
grep("^RP[LS]",rownames(pbmc),value = TRUE)

pbmc[["percent.rbp"]] <- PercentageFeatureSet(pbmc, pattern = "^RP[LS]")


head(pbmc@meta.data, 5)
# In order we have: Cell barcode, how many reads associated to the cell (library size
# of that cell), % of reads mapping on MT genes and % of reads mapping on ribosomal
# protein genes.
# The numbers I see are normal in the context of single cell RNA-Seq.

# Visualize QC metrics as violin plots - also adding the RPL genes
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), 
        ncol = 4)
# Without dots:
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), 
        ncol = 4, pt.size=0)

# FeatureScatter is typically used to visualize feature-feature relationships, 
# but can be used for anything calculated by the object, i.e. columns in 
# object metadata, PC scores etc.
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# If I find too few reads coming from a cell, we might be looking at an empty droplet,
# instead if I find too many we might be looking at a doublet (or triplet etc...)
# From the plot on the right, I see a correlation between nCount and nFeature,
# so I can find a cut point to remove empty droplet and multiplet.
# To understand where to cut, I can look at the previous violin plots.

plot3 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.rbp")
plot3

# We keep the cells with at least 200 genes expressed (we remove empty droplets) 
# and with less than 2500 count reads (we remove multiplets) and % of reads
# mapping on MT genes lower than 5.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# This a rule of thumb, there is no 'mathematical formula' to choose 
# the above thresholds.

pbmc

### Normalization (we will obtain log2(CP10K) values)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# The original and normalized counts are buried inside the Seurat object pbmc. 
# Let us find them:
pbmc@assays
pbmc@assays$RNA

# Now that we have normalized the expression values, we can take a 
# look to the genes that have the highest mean expression across our cells:
apply(pbmc[["RNA"]]$data,1,mean) -> gene.expression

sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=50) 
# Those are the 50 genes with the highest average expression in our data.
# The record holder is usually MALAT1: it's a long non coding RNA but is polyadenilated 
# (overexpressed by lung cancer cells) and is expressed everywhere, not only
# by that specific metastasis of lung cancer.

# Distribution of the expression of MALAT1 and of another typical housekeeping gene:
VlnPlot(pbmc, features = c("MALAT1","GAPDH"))

# GAPDH is an housekeeping gene, but some cells have expression 0 for it
# (typical example of dropout): this could be because of:
# - Biological reasons: in that moment the gene was not expressed
# - Technical reasons: we had few reads coming from this gene ans we were not able
#   to capture any of them

cc.genes.updated.2019
# Precomputed list of genes expressed only (or expressed more) in specific cell
# cycle phases

CellCycleScoring(pbmc, s.features = cc.genes.updated.2019$s.genes, 
                 g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> pbmc
# The function tells, for each cell, the likelihood of being in any cell cycle phase

pbmc[[]]

# The default method -vst- computes (or better, estimates) the mean-variance 
# relationship of each gene, and chooses the 2000 genes with the highest variance 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

x11()
plot1 + plot2

# After choosing the 2000 most variable genes, there is the scaling
# in order to have each gene with mean 0 and variance 1 (across all the cells).
# Another table is added to the data structure, containing the scaled values.

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc@assays$RNA

# But before proceeding another scaling of the counts is advised. 
# The idea is to shift the expression of each gene, so that the mean 
# expression across cells is 0 and the variance across cells is 1 This step 
# gives equal weight in downstream analyses, so that highly-expressed genes 
# do not dominate. In practice, values are sort of “binarized”, or rather 
# “ternarized” - >0 “high expression”, 0 “average expression” <0 “under 
# expression, or no expression at all”. Notice that this is done for all the 
# genes - not only the most variable ones.

# We perform this step because in this way the PCA will perform better,
# but why the PCA will perform better?
# - Scaling means I correct to have the same mean but also the same variance so
# is not just a shifting.
# Also in this case, as in bulk, the most variable genes are the ones with the
# highest expression, so if I don't do the scaling the first PCs are driven by the
# genes with the higher expression, which are the RP genes --> with scaling, 
# I consider all most variable genes, regardless their expression level.
# In this way PCA performs better.

# The VariableFeatures are only the top 2000 genes
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# Many of these genes are those that at the end will be marker genes which will 
# tell us which cell type is each cluster.
# The one showned are the most relevant genes for each PC.

VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")

DimPlot(pbmc, reduction = "pca") 
# By default coloured by predicted cell cycle phase
# ==> the colours are mixed, so the cell cycle phase is not a key factor
# dividing my cells.

# with ndims we can choose how many PC to plot
ElbowPlot(pbmc, ndims=20)
# The elbow is more or less around 9/10 Principal Components.
# Before the plateaux, it's all biological variability; after, it's 
# just technical random background noise.

# In many cases we choose the number of PC that explains 75% of the variance
pc.touse <- (pbmc$pca@stdev)^2
pc.touse <- pc.touse/sum(pc.touse)
pc.touse <- cumsum(pc.touse)[1:50]
pc.touse <- min(which(pc.touse>=0.75))
pc.touse
# In this case the number of PCs chosen is 28.

# I can try both 10 and 28 and see which one gives the best results.

### Clustering
pbmc <- FindNeighbors(pbmc, dims = 1:10)

pbmc <- FindClusters(pbmc, resolution = 0.5) 
# Usually the resolution is between 0.4 and 1.2 (higher resolution means
# more and smaller clusters)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
head(pbmc[[]],5)
DimPlot(pbmc, reduction = "pca")

# But we know that for visualization and 2D plotting there are better 
# strategies, like t-SNE, done on the PC dimensions chosen for clustering:
pbmc <- RunTSNE(pbmc, dims=1:10) # Projection of 10 PCs used for clustering into 2D
DimPlot(pbmc, reduction = "tsne")

# If you haven't installed UMAP, you can do so via 
# reticulate::py_install(packages = 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)
# if you cannot install UMAP, t-SNE is anyway ok for your project!
DimPlot(pbmc, reduction = "umap")

# In this case we have 3 clouds of points: I can't say that one cloud is more
# close to another wrt to the other cloud!! I can just say they are separate.

# Distribution of number of reads per cell divided by cluster 
# (in cluster 8 I tend to have less reads)
VlnPlot(pbmc,features="nCount_RNA")
# It may be that cluster 8 is formed just by grouping cells with not too many reads
# which shared technical problems.
# So a solution is to say that cluster 8 is a problematic cluster.

VlnPlot(pbmc,features="nFeature_RNA")

VlnPlot(pbmc,features="percent.mt")

VlnPlot(pbmc,features="percent.rbp")
# In cluster 8 it's very low 

# After all these checks I identify cluster 8 as problematic.

library(ggplot2)
pbmc@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")

# Find all markers of cluster 2 versus all the others
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25, test.use = "wilcox")
head(cluster2.markers, n = 5)

# FindMarkers function: 
# ident.1 = 2: it means that I want to compare cluster 2 with all the others
# min.pct = 0.25: to be considered a marker gene of a cluster, it must be expressed
# in at least 25% of the cells of the cluster
# wilcox: it's the Mann-Whitney U Test

# I find DE genes for each cluster (each cluster is compared to all the others)
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
# Do I look for the best 5 in terms of log2FC or fdr?
# I try both and I pick the genes which best characterize the cell type,
# because for some clusters may be better to choose depending on log2FC,
# for others instead on adjusted p-value.

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", 
                               "FCGR3A", "LYZ", "PPBP", "CD8A"))

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Here we show the top10 genes (by log2FC) of each cluster.
# Yellow: expression goes up
# Purple: expression goes down
# ==> I can see for each cluster that the respective best 10 goes up in the
# respective cluster.
# Cluster 1,3,7 and 8 are ok, but the couples 0/2 and 4/6 have problems,
# because the selected genes going up in 4 goes up also in 6 (same for 0 and 2).

# So they seem to be characterized by the same DE genes.
# Are 0 and 2 the same cluster splitted in 2 by mistake or are they very similar
# cells but there is something more making the difference between cluster 0 and 2?
# Maybe hidden in this plot there is 1 gene going up in 0 but not in 2; 
# likewiese in 4 and 6
# In fact, looking at the UMAP plot, I can see that they are close in the space.

# Now let's assume that 0 and 2 are a unique cell type.
# I perform DE gene analysis considering 0+2 as one cluster against all the others,
# because if they are the same cell type I will find DE genes characterizing both.
# In this case, if I find a very clear marker gene, I consider them to be
# the same cell type.
# Now we do the second check: I do a comparison 0 vs 2 ==> I find DE genes between
# cluster 0 and 2: above we were considering them together to see if there were
# markers telling me they were the same cell type, while here I compare them so to
# see if there are genes making differences between one cluster and the other.
# If I find some genes here, I can for example say they are different subtype
# of the same cell type.

# Here Prof. skipped to the end but the final plot is the UMAP with the cell
# types names (see on Vignetta).

