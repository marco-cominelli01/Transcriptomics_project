library(dplyr)
library(Seurat)
library(patchwork)

load("./SRA779509_SRS3805265.sparse.RData")

for (i in 1:length(row.names(sm))) {
  gene <- row.names(sm)[i]
  if (is.character(gene)) {
    if (grepl("ENSG", gene)) {
      name <- strsplit(gene, "_ENSG")[[1]]
      row.names(sm)[i] <- name[1]
    } else if (grepl("ENSMU", gene)) {
      name <- strsplit(gene, "_ENSMU")[[1]]
      row.names(sm)[i] <- name[1]
    } else {
      row.names(sm)[i] <- gene
    }
  }
}

row.names(sm) <- make.unique(as.character(row.names(sm)))

rm(name)
rm(gene)
rm(i)

bmmc <- CreateSeuratObject(counts = sm, project = "bone_marrow", 
                           min.cells = 3, min.features = 200)


head(colnames(bmmc))

bmmc[["percent.mt"]] <- PercentageFeatureSet(bmmc, pattern = "^MT-")

bmmc[["percent.rbp"]] <- PercentageFeatureSet(bmmc, pattern = "^RP[LS]")

#VlnPlot(bmmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)

bmmc <- subset(bmmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


# Normalize counts
bmmc <- NormalizeData(bmmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Genes with highest mean expression across cells
apply(bmmc[["RNA"]]$data,1,mean) -> gene.expression
sort(gene.expression, decreasing = TRUE) -> gene.expression
head(gene.expression, n=50)

# VlnPlot(bmmc, features = c("ACTB", "MALAT1"), pt.size = 0)

CellCycleScoring(bmmc, s.features = cc.genes.updated.2019$s.genes, 
                 g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = TRUE) -> bmmc

bmmc[[]]

bmmc <- FindVariableFeatures(bmmc, selection.method = "vst", nfeatures = 2000)
# top10 <- head(VariableFeatures(bmmc), 10)

all.genes <- rownames(bmmc)
bmmc <- ScaleData(bmmc, features = all.genes)

# PCA (done on 2000 most variable genes)
bmmc <- RunPCA(bmmc, features = VariableFeatures(object = bmmc))

print(bmmc[["pca"]], dims = 1:5, nfeatures = 5)

VizDimLoadings(bmmc, dims = 1:2, reduction = "pca")

#DimPlot(bmmc, reduction = "pca")

#ElbowPlot(bmmc, ndims=20)

# Clustering

bmmc8 <- FindNeighbors(bmmc, dims = 1:8)
bmmc8_05 <- FindClusters(bmmc8, resolution = 0.5)
bmmc8_0.7 <- FindClusters(bmmc8, resolution = 0.7)
bmmc8_0.8 <- FindClusters(bmmc8, resolution = 0.8)
bmmc8_1.2 <- FindClusters(bmmc8, resolution = 1.2)

#bmmc10 <- FindNeighbors(bmmc, dims = 1:10)
#bmmc10_05 <- FindClusters(bmmc10, resolution = 0.5)
#bmmc10_0.8 <- FindClusters(bmmc10, resolution = 0.8)
#bmmc10_1.2 <- FindClusters(bmmc10, resolution = 1.2)


#x11()
#DimPlot(bmmc10_05, reduction = "pca")
#x11()
#DimPlot(bmmc8_05, reduction = "pca")

#bmmc_tsne8 <- RunTSNE(bmmc8_05, dims=1:8)
#DimPlot(bmmc_tsne8, reduction = "tsne", pt.size = 1)

#bmmc_tsne10 <- RunTSNE(bmmc10_05, dims=1:10)
#DimPlot(bmmc_tsne10, reduction = "tsne", pt.size = 1)

#bmmc_tsne10_1.2 <- RunTSNE(bmmc10_1.2, dims=1:10)
#DimPlot(bmmc_tsne10_1.2, reduction = "tsne", pt.size = 1)

## UMAP


#x11()
#bmmc_umap10 <- RunUMAP(bmmc10_05, dims = 1:10)
#DimPlot(bmmc_umap10, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap10_0.8 <- RunUMAP(bmmc10_0.8, dims = 1:10)
#DimPlot(bmmc_umap10_0.8, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap10_1.2 <- RunUMAP(bmmc10_1.2, dims = 1:10)
#DimPlot(bmmc_umap10_1.2, reduction = "umap", pt.size = 1)

#x11()
bmmc_umap8 <- RunUMAP(bmmc8_05, dims = 1:8)
#DimPlot(bmmc_umap8, reduction = "umap", pt.size = 1)

x11()
bmmc_umap8_0.8 <- RunUMAP(bmmc8_0.8, dims = 1:8)
DimPlot(bmmc_umap8_0.8, reduction = "umap", pt.size = 1)

#x11()
bmmc_umap8_1.2 <- RunUMAP(bmmc8_1.2, dims = 1:8)
DimPlot(bmmc_umap8_1.2, reduction = "umap", pt.size = 1)

## Quality controls on clustering
#x11()
#VlnPlot(bmmc8_05,features="nCount_RNA")

#x11()
#VlnPlot(bmmc8_05,features="nFeature_RNA")

#x11()
#VlnPlot(bmmc8_05,features="percent.mt")

#x11()
#VlnPlot(bmmc8_05,features="percent.rbp")

# This function calculates the overlap of cells between clusters from two different clustering results and visualizes the output as an heatmap.
# The brighter the color, the greater the number of cells in common between the two clusters.
# Each value is normalized by the total number of cells in the respective row cluster.
# This visualization can help us retrieve information about the original cluster from which the new clusters are derived.
# In addition, the heatmap can guide us in the cell type labelling.

library(viridisLite)
library(viridis)


overlap_clusters <- function(clustering1, clustering2){
  clusters1 <- as.numeric(levels(clustering1@meta.data$seurat_clusters))
  clusters2 <- as.numeric(levels(clustering2@meta.data$seurat_clusters))
  
  
  mat <- matrix(0, nrow = length(clusters2), ncol=length(clusters1), dimnames = list(levels(clustering2@meta.data$seurat_clusters), levels(clustering1@meta.data$seurat_clusters)))
  
  for (cluster1 in clusters1){
    for (cluster2 in clusters2){
      cells1 <- rownames(clustering1@meta.data[which(clustering1@meta.data$seurat_clusters==cluster1),])
      n1 <- length(cells1)
      cells2 <- rownames(clustering2@meta.data[which(clustering2@meta.data$seurat_clusters==cluster2),])
      n2 <- length(cells2)
      overlap <- intersect(cells1,cells2)
      mat[as.character(cluster2), as.character(cluster1)]<- length(overlap) #/min(n1,n2)
    }
  }
  
  
  df <- as.data.frame(as.table(mat))
  for (i in unique(df$Var1)){
    tot_col = sum(df[df$Var1==i,3])
    df[df$Var1==i,3] <- df[df$Var1==i,3]/tot_col
  }
  
  heatmap <- ggplot(df, aes(x = Var2, y = Var1, fill = Freq)) +
    geom_tile(color = "white", lwd = 0.55, linetype = 1) +
    scale_fill_viridis(option = "viridis", name = "Similarity", discrete = FALSE) +
    labs(title = "Overlap between clusters") +
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(color = "black")) +
    coord_fixed() +
    guides(fill = guide_colourbar(ticks = FALSE, title = NULL, ticks.linewidth = 0)) +
    scale_x_discrete(position = "bottom")
  return(heatmap)
}

library(ggplot2)
#overlap_clusters(bmmc8_05, bmmc10_05)


pbmc.markers <- FindAllMarkers(bmmc8_05, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

top10_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

View(top5_each_cluster_by_logfc)
View(top10_each_cluster_by_logfc)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

#DoHeatmap(bmmc8_05, features = top10$gene) + NoLegend()

# Top 10 for cluster 0 

x11()
FeaturePlot(bmmc_umap8, features = c("LDHBP2", "LEF1", "CCR7", "TCF7", "BCL11B", "CD27", 
                                     "TRAT1", "PIK3IP1", "IL7R", "AC107075.1"), reduction = 'umap')

# Top 10 for cluster 1

x11()
FeaturePlot(bmmc_umap8, features = c("KLRF1", "CCL3", "KLRC2", "PRF1", "SPON2", "GNLY", 
                                     "FGFBP2", "KLRD1", "CLIC3", "GZMB"), reduction = 'umap')


# Top 10 for cluster 2

x11()
FeaturePlot(bmmc_umap8, features = c("ARHGAP24", "BANK1", "IGLC3", "TNFRSF13C", "LINC00926", 
                                     "MS4A1", "ADAM28", "CD79A", "IGHD", "IGHG4"), reduction = 'umap')

# Top 10 for cluster 3

x11()
FeaturePlot(bmmc_umap8, features = c("GZMK", "A2M", "LINC01871", "KLRG1", "DUSP2", "LYAR", 
                                     "CCL5", "TRGC2", "IL32", "GZMA"), reduction = 'umap')

# Top 10 for cluster 4

x11()
FeaturePlot(bmmc_umap8, features = c("CLEC4E", "LYZ", "QPCT", "S100A12", "CSTA", 
                                     "S100A8", "CYP1B1", "VCAN", "S100A9", "MGST1"), reduction = 'umap')

# Top 10 for cluster 5

x11()
FeaturePlot(bmmc_umap8, features = c("FXYD5", "ACTG1", "ZFP36L2", "LTB", "SNU13", 
                                     "EIF3H", "EIF3D", "OCIAD2", "TRBC2", "MCUB"), reduction = 'umap')

# Top 10 for cluster 6

x11()
FeaturePlot(bmmc_umap8, features = c("ITGB1", "LIMS1", "AC058791.1", "TTC39C", "AC016831.5", 
                                     "LDHBP2", "C16orf54", "HNRNPH1", "TNRC6B", "ZC3HAV1"),reduction = 'umap')

# Top 10 for cluster 7

x11()
FeaturePlot(bmmc_umap8, features = c("CENPF", "TOP2A", "MYL4", "CDK1", "PCLAF", 
                                     "TYMS", "CA1", "GYPA", "SMIM1", "KLF1"), reduction = 'umap')
# Top 10 for cluster 8

x11()
FeaturePlot(bmmc_umap8, features = c("LINC00570", "IFIT1B", "TMCC2", "HBA1", "HBA2",
                                     "KRT1", "SNCA", "GMPR", "TRIM58", "OSBP2"), reduction = 'umap')

# Top 10 for cluster 9

x11()
FeaturePlot(bmmc_umap8, features = c("LYPD2", "CDKN1C", "VMO1", "AC064805.1", "HES4", 
                                     "TCF7L2", "AC245884.12", "CSF1R", "MS4A7", "SIGLEC10"), reduction = 'umap')

# Top 10 for cluster 10

x11()
FeaturePlot(bmmc_umap8, features = c("IGKC", "MZB1", "JCHAIN", "ITM2C", "SPNS1", "FXYD5", 
                                     "CHMP1B", "IGHM", "CD74", "CNOT4"), reduction = 'umap')

# Top 10 for cluster 11

x11()
FeaturePlot(bmmc_umap8, features = c("AL357143.1", "SCT", "AC097375.1", "PPP1R14BP2", 
                                     "CLEC4C", "LILRA4", "NRP1", "TNFRSF21", "PPP1R14BP3", "DNASE1L3")
            ,reduction = 'umap')

### Cluster 5 + 10 vs all
cluster5AND10.markers <- FindMarkers(bmmc8_05, ident.1 = c(5,10), min.pct = 0.25, test.use = "wilcox")
cluster5AND10.markers <- cluster5AND10.markers[order(-cluster5AND10.markers$avg_log2FC),]
head(cluster5AND10.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8, features = c("IGKC", "FXYD5", "ACTG1", "LDHB", "ZFP36L2", "LTB", 
                                     "CALM1", "TRAC", "EIF3H", "SNU13"), reduction = 'umap')


# cluster 5 vs 10
cluster5_10.markers <- FindMarkers(bmmc8_05, ident.1 = 5, ident.2 = 10, min.pct = 0.25, test.use = "wilcox")
cluster5_10.markers <- cluster5_10.markers[order(-cluster5_10.markers$avg_log2FC),]
head(cluster5_10.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8, features = c("TRAC", "IL32", "NKG7", "CCL5", "TRBC2", "CD3D", 
                                     "LDHB", "SARAF", "PTPRC", "GNLY"), reduction = 'umap')


# Cluster 5 vs cluster 5 (splitted in two by looking at bmmc8_0.8)
cluster5a_5b.markers <- FindMarkers(bmmc_umap8_0.8, ident.1 = 6, ident.2 = 12, min.pct = 0.25, test.use = "wilcox")
cluster5a_5b.markers <- cluster5a_5b.markers[order(-cluster5a_5b.markers$avg_log2FC),]
head(cluster5a_5b.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("EEF1B2", "LTB", "GAS5", "LRRC75A-AS1", "SKP1", "NPM1", 
                                     "RPL7A", "RPS8", "RPL8", "PABPC1"), reduction = 'umap')


# Cluster 6 vs all others (bmmc8_0.8)
cluster6_all.markers <- FindMarkers(bmmc_umap8_0.8, ident.1 = 6, min.pct = 0.25, test.use = "wilcox")
head(cluster6_all.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("AL136968.2", "AC016739.1", "PFN1P1"), reduction = 'umap')

# Cluster 12 vs all others (bmmc8_0.8)
cluster12_all.markers <- FindMarkers(bmmc_umap8_0.8, ident.1 = 12, min.pct = 0.25, test.use = "wilcox")
head(cluster12_all.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("B2M", "NKG7", "EEF1A1", "CCL5"), reduction = 'umap')


# Cluster 1 + 2 + 12 vs all (bmm8_0.8)
cluster1AND2_12.markers <- FindMarkers(bmmc8_0.8, ident.1 = c(1,2,12), min.pct = 0.25, test.use = "wilcox")
cluster1AND2_12.markers <- cluster1AND2_12.markers[order(-cluster1AND2_12.markers$avg_log2FC),]
head(cluster1AND2_12.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("KLRF1", "KLRD1", "TRGC1", "PRF1",
                                         "FGFBP2", "TRGC2", "CMC1","GZMH","GNLY","TRDC"), reduction = 'umap')


# Cluster 6 + 9 vs all (bmm8_0.8)
cluster6AND0.markers <- FindMarkers(bmmc8_0.8, ident.1 = c(6,0), min.pct = 0.25, test.use = "wilcox")
cluster6AND0.markers <- cluster6AND0.markers[order(-cluster6AND0.markers$avg_log2FC),]
head(cluster6AND0.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("CCR7", "LEF1", "TCF7", "TRAT1",
                                         "PIK3IP1", "BCL11B", "IL7R","LDHB","LDHBP2","NOSIP"), reduction = 'umap')
            
            
# Cluster 6 vs 9 (bmm8_0.8)
cluster6_9.markers <- FindMarkers(bmmc_umap8_0.8, ident.1 = 6, ident.2 = 9, min.pct = 0.25, test.use = "wilcox")
cluster6_9.markers <- cluster6_9.markers[order(-cluster6_9.markers$avg_log2FC),]
head(cluster6_9.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_0.8, features = c("IL32", "LDHB", "PTPRC", "CD3D", "MT-ND5", "SARAF", 
                                         "GAS5", "BTG1", "CALM1", "NPM1"), reduction = 'umap')


# Cluster 9 + 10 vs all (bmm8_1.2)
cluster9AND10.markers <- FindMarkers(bmmc8_1.2, ident.1 = c(9,10), min.pct = 0.25, test.use = "wilcox")
cluster9AND10.markers <- cluster9AND10.markers[order(-cluster9AND10.markers$avg_log2FC),]
head(cluster9AND10.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_1.2, features = c("HBA1", "HBA2", "IFIT1B", "GYPA",
                                         "HBB", "HEMGN", "HBD","HBM","SNCA","ALAS2"), reduction = 'umap')


# Cluster 11 + 12 vs all (bmm8_1.2)
cluster11AND12.markers <- FindMarkers(bmmc8_1.2, ident.1 = c(11,12), min.pct = 0.25, test.use = "wilcox")
cluster11AND12.markers <- cluster11AND12.markers[order(-cluster11AND12.markers$avg_log2FC),]
head(cluster11AND12.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_1.2, features = c("IFI27", "IGKC", "BLVRB", "AHSP",
                                         "MPC2", "CA1", "ATPIF1","TUBB","PRDX2","HMGB2", "MZB1"), reduction = 'umap')

# Cluster 11 vs 9 (bmm8_1.2)
cluster11_9.markers <- FindMarkers(bmmc_umap8_1.2, ident.1 = 11, ident.2 = 9, min.pct = 0.25, test.use = "wilcox")
cluster11_9.markers <- cluster11_9.markers[order(-cluster11_9.markers$avg_log2FC),]
head(cluster11_9.markers, n = 20)

x11()
FeaturePlot(bmmc_umap8_1.2, features = c("IFI27", "NDUFS5", "EEF1B2", "ATPIF1", "NCL", "NDUFB9", 
                                         "EIF3K", "TOMM7", "MZB1"), reduction = 'umap')


# Cluster 11 vs all (bmm8_1.2)
cluster11.markers <- FindMarkers(bmmc8_1.2, ident.1 = 11, min.pct = 0.25, test.use = "wilcox")
cluster11.markers <- cluster11.markers[order(-cluster11.markers$avg_log2FC),]
head(cluster11.markers, n = 10)

x11()
FeaturePlot(bmmc_umap8_1.2, features = c("PCLAF", "BLVRB", "IFI27", "AHSP", "TMEM14C", "MPC2", 
                                         "HIST1H4C", "NUCB2", "CA1", "ATPIF1"), reduction = 'umap')

