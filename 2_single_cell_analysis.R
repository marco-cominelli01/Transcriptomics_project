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

DimPlot(bmmc, reduction = "pca")

ElbowPlot(bmmc, ndims=20)

# Clustering

bmmc8 <- FindNeighbors(bmmc, dims = 1:8)
bmmc8_05 <- FindClusters(bmmc8, resolution = 0.5)
bmmc8_0.8 <- FindClusters(bmmc8, resolution = 0.8)
bmmc8_1.2 <- FindClusters(bmmc8, resolution = 1.2)

bmmc10 <- FindNeighbors(bmmc, dims = 1:10)
bmmc10_05 <- FindClusters(bmmc10, resolution = 0.5)
bmmc10_0.8 <- FindClusters(bmmc10, resolution = 0.8)
bmmc10_1.2 <- FindClusters(bmmc10, resolution = 1.2)


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


x11()
bmmc_umap10 <- RunUMAP(bmmc10_05, dims = 1:10)
DimPlot(bmmc_umap10, reduction = "umap", pt.size = 1)

x11()
bmmc_umap10_0.8 <- RunUMAP(bmmc10_0.8, dims = 1:10)
DimPlot(bmmc_umap10_0.8, reduction = "umap", pt.size = 1)

x11()
bmmc_umap10_1.2 <- RunUMAP(bmmc10_1.2, dims = 1:10)
DimPlot(bmmc_umap10_1.2, reduction = "umap", pt.size = 1)

x11()
bmmc_umap8 <- RunUMAP(bmmc8_05, dims = 1:8)
DimPlot(bmmc_umap8, reduction = "umap", pt.size = 1)

x11()
bmmc_umap8_0.8 <- RunUMAP(bmmc8_0.8, dims = 1:8)
DimPlot(bmmc_umap8_0.8, reduction = "umap", pt.size = 1)

x11()
bmmc_umap8_1.2 <- RunUMAP(bmmc8_1.2, dims = 1:8)
DimPlot(bmmc_umap8_1.2, reduction = "umap", pt.size = 1)

## Quality controls on clustering
x11()
VlnPlot(bmmc8_05,features="nCount_RNA")

x11()
VlnPlot(bmmc8_05,features="nFeature_RNA")

x11()
VlnPlot(bmmc8_05,features="percent.mt")

x11()
VlnPlot(bmmc8_05,features="percent.rbp")

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
overlap_clusters(bmmc8_05, bmmc10_05)


pbmc.markers <- FindAllMarkers(bmmc8_05, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top5_each_cluster_by_logfc <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

View(top5_each_cluster_by_logfc)


pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(bmmc8_05, features = top10$gene)

##########

library(Seurat)
library(dplyr)

# Supponiamo che 'pbmc.markers' sia già stato creato e contenga i marker differenziali
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

# Crea la heatmap con parametri aggiuntivi per migliorare la visualizzazione
DoHeatmap(bmmc8_05, features = top10$gene) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10))

# Se le colonne sono ancora schiacciate, prova a regolare il layout della heatmap
DoHeatmap(bmmc8_05, features = top10$gene, slot = "scale.data") +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10)) +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )





