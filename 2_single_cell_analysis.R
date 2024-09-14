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

VlnPlot(bmmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rbp"), ncol = 4, pt.size=0)

bmmc <- subset(bmmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 8)


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

DimPlot(bmmc, reduction = "pca")

ElbowPlot(bmmc, ndims=20)

# Clustering

#bmmc8 <- FindNeighbors(bmmc, dims = 1:8)
#bmmc8_0.5 <- FindClusters(bmmc8, resolution = 0.5)
#bmmc8_0.7 <- FindClusters(bmmc8, resolution = 0.7)
#bmmc8_0.8 <- FindClusters(bmmc8, resolution = 0.8)
#bmmc8_1.2 <- FindClusters(bmmc8, resolution = 1.2)

bmmc10 <- FindNeighbors(bmmc, dims = 1:10)
bmmc10_0.5 <- FindClusters(bmmc10, resolution = 0.5)
#bmmc10_0.8 <- FindClusters(bmmc10, resolution = 0.8)
#bmmc10_1.2 <- FindClusters(bmmc10, resolution = 1.2)

bmmc12 <- FindNeighbors(bmmc, dims = 1:12)
bmmc12_0.5 <- FindClusters(bmmc12, resolution = 0.5)


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
#bmmc_umap8_0.5 <- RunUMAP(bmmc8_0.5, dims = 1:8)
#DimPlot(bmmc_umap8_0.5, reduction = "umap", pt.size = 1)

x11()
bmmc_umap10_0.5 <- RunUMAP(bmmc10_0.5, dims = 1:10)
DimPlot(bmmc_umap10_0.5, reduction = "umap", pt.size = 1)

bmmc_umap12_0.5 <- RunUMAP(bmmc12_0.5, dims = 1:12)
#DimPlot(bmmc_umap12_0.5, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap10_1.2 <- RunUMAP(bmmc10_1.2, dims = 1:10)
#DimPlot(bmmc_umap10_1.2, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap8 <- RunUMAP(bmmc8_05, dims = 1:8)
#DimPlot(bmmc_umap8, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap8_0.8 <- RunUMAP(bmmc8_0.8, dims = 1:8)
#DimPlot(bmmc_umap8_0.8, reduction = "umap", pt.size = 1)

#x11()
#bmmc_umap8_1.2 <- RunUMAP(bmmc8_1.2, dims = 1:8)
#DimPlot(bmmc_umap8_1.2, reduction = "umap", pt.size = 1)

## Quality controls on clustering
#x11()
VlnPlot(bmmc10_0.5,features="nCount_RNA")

#x11()
VlnPlot(bmmc10_0.5,features="nFeature_RNA")

#x11()
VlnPlot(bmmc10_0.5,features="percent.mt")

#x11()
VlnPlot(bmmc10_0.5,features="percent.rbp")

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
overlap_clusters(bmmc10_0.5, bmmc12_0.5)

overlap_clusters(bmmc12_0.5, bmmc10_0.5)


pbmc.markers <- FindAllMarkers(bmmc10_0.5, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

top10_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

top50_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 50, order_by = avg_log2FC)

View(top5_each_cluster_by_logfc)
View(top10_each_cluster_by_logfc)
View(top50_each_cluster_by_logfc)

pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(bmmc10_0.5, features = top10$gene) + NoLegend()

# Top 10 for cluster 0 

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("LDHBP2", "LEF1", "CD27", "BCL11B", "TCF7", "CCR7", 
                                          "TRAT1", "PIK3IP1", "IL7R", "AC010343.1"), reduction = 'umap')

# Top 10 for cluster 1

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("KLRC2", "TRGC2", "A2M", "KLRG1", "TRGC1", "GZMH", 
                                          "CCL5", "CMC1", "TRDC"), reduction = 'umap')


# Top 10 for cluster 2

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("TNFRSF13C", "BANK1", "IGLC3", "TNFRSF13C", "LINC00926", 
                                     "MS4A1", "ADAM28", "IGLC2", "IGHD", "IGHG4", "ARHGAP24", "CD79A"), reduction = 'umap')

# Top 10 for cluster 3

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("SPON2", "IGFBP7", "CLIC3", "AKR1C3", "CCL3", "KLRF1", 
                                     "PRF1", "PRSS23", "FGFBP2", "GZMB", "S100A4", "CCR7", "IL7R"), reduction = 'umap')

# Top 10 for cluster 4

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("LYZ", "CD14", "QPCT", "S100A12", "CSTA", 
                                     "S100A8", "CYP1B1", "VCAN", "S100A9", "MGST1"), reduction = 'umap')

# Top 10 for cluster 5

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("FXYD5", "ACTG1", "ZFP36L2", "LTB", "IL7R", 
                                     "LCK", "TMEM123", "MCUB", "EIF3H", "SNU13"), reduction = 'umap')

# Top 10 for cluster 6

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("PCLAF", "TYMS", "GYPA", "CA1", "SMIM1", 
                                     "CA2", "KLF1", "HEMGN", "HBM", "AHSP"),reduction = 'umap')

# Top 10 for cluster 7

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("MS4A6A", "MNDA", "LST1", "S100A9", "CST3", 
                                     "AIF1", "VSIG4", "S100A12", "FCN1", "PYCARD", "TREM2"), reduction = 'umap')
# Top 10 for cluster 8

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("IGKC", "JCHAIN", "IGHM", "CD74", "SPNS1",
                                     "HLA-DPA1", "HLA-DPB1", "HLA-DRA"), reduction = 'umap')

# Top 10 for cluster 9

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("IFIT1B", "TMCC2", "HBA1", "HBA2", "KRT1", 
                                     "SNCA", "PDZK1IP1", "TRIM58", "HBB"), reduction = 'umap')


# Top 10 for cluster 10

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("LYPD2", "CDKN1C", "VMO1", "TCF7L2", "CSF1R", "HES4", 
                                     "SIGLEC10", "MS4A7", "FCGR3A"), reduction = 'umap')

# Top 10 for cluster 11

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("SCT", "PPP1R14BP2", "GAS6", "FCER1A","CST3",
                                     "CLEC4C", "LILRA4", "NRP1", "TNFRSF21", "PPP1R14BP3", "DNASE1L3")
            ,reduction = 'umap')


# Top 10 for cluster 12

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("FCER1A", "CLEC10A", "PKIB", "FLT3", 
                                     "CACNA2D3", "CDC1C", "LGALS2", "NDRG2", "HLA-DQB2", "HLA-DQA1")
            ,reduction = 'umap')



# Cluster 8 vs 2

cluster8_2.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 8, ident.2 = 2, min.pct = 0.25, test.use = "wilcox")
cluster8_2.markers <- cluster8_2.markers[order(-cluster8_2.markers$avg_log2FC),]
head(cluster8_2.markers, n = 50)

# Cluster 1 vs all
cluster1_all.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 1, min.pct = 0.25, test.use = "wilcox")
cluster1_all.markers <- cluster1_all.markers[order(-cluster1_all.markers$avg_log2FC),]
head(cluster1_all.markers, n = 50)

# Cluster 1 vs 0

cluster1_0.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 1, ident.2 = 0, min.pct = 0.25, test.use = "wilcox")
cluster1_0.markers <- cluster1_0.markers[order(-cluster1_0.markers$avg_log2FC),]
head(cluster1_0.markers, n = 50)

# Cluster 3 vs all
cluster3_all.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 3, min.pct = 0.25, test.use = "wilcox")
cluster3_all.markers <- cluster3_all.markers[order(-cluster3_all.markers$avg_log2FC),]
head(cluster3_all.markers, n = 50)


# Cluster 8 vs all
cluster8_all.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 8, min.pct = 0.25, test.use = "wilcox")
cluster8_all.markers <- cluster8_all.markers[order(-cluster8_all.markers$avg_log2FC),]
head(cluster8_all.markers, n = 50)

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("IGKC", "JCHAIN", "IGHM", "IGHA1")
            ,reduction = 'umap')


# Cluster 5 vs all
cluster5_all.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 5, min.pct = 0.25, test.use = "wilcox")
cluster5_all.markers <- cluster5_all.markers[order(-cluster5_all.markers$avg_log2FC),]
head(cluster5_all.markers, n = 50)

x11()
FeaturePlot(bmmc_umap10_0.5, features = c("S100A4")
            ,reduction = 'umap')

# Cluster 5 vs 0

cluster5_0.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 5, ident.2 = 0, min.pct = 0.25, test.use = "wilcox")
cluster5_0.markers <- cluster5_0.markers[order(-cluster5_0.markers$avg_log2FC),]
head(cluster5_0.markers, n = 50)
x11()
FeaturePlot(bmmc_umap10_0.5, features = c("NKG7", "FXYD5", "ACTG1", "CCL5")
            ,reduction = 'umap')

# Cluster 6 vs 9
cluster6_9.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 6, ident.2 = 9, min.pct = 0.25, test.use = "wilcox")
cluster6_9.markers <- cluster6_9.markers[order(-cluster6_9.markers$avg_log2FC),]
head(cluster6_9.markers, n = 50)

# Cluster 9 vs 6
cluster9_6.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 9, ident.2 = 6, min.pct = 0.25, test.use = "wilcox")
cluster9_6.markers <- cluster9_6.markers[order(-cluster9_6.markers$avg_log2FC),]
head(cluster9_6.markers, n = 50)

# Cluster 11 vs 12
cluster11_12.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 11, ident.2 = 12, min.pct = 0.25, test.use = "wilcox")
cluster11_12.markers <- cluster11_12.markers[order(-cluster11_12.markers$avg_log2FC),]
head(cluster11_12.markers, n = 50)

# Cluster 12 vs 11
cluster12_11.markers <- FindMarkers(bmmc_umap10_0.5, ident.1 = 12, ident.2 = 11, min.pct = 0.25, test.use = "wilcox")
cluster12_11.markers <- cluster12_11.markers[order(-cluster12_11.markers$avg_log2FC),]
head(cluster12_11.markers, n = 50)

# Marker genes plot for presentation
FeaturePlot(bmmc_umap10_0.5, features = c("TSPAN13", "DERL3", "CLEC10A", "CD1C")
            ,reduction = 'umap')







# Final plot with names of cell types


new.cluster.ids <- c("Naive CD4+ T", "NK", "Naive B", "NK", "CD14+ Mono", "T-cells", "Erythroid",
                     "Macrophages", "Plasma", "Erythroid", "FCGR3A+ Mono", "Plasmacytoid DC", "Dendritic")

names(new.cluster.ids) <- levels(bmmc_umap10_0.5)
bmmc_umap10_0.5 <- RenameIdents(bmmc_umap10_0.5, new.cluster.ids)
DimPlot(bmmc_umap10_0.5, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()















