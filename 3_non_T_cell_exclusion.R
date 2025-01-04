#Load the necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(knitr)
library(Rgb)
library(reshape2)
library(tidyverse)


#Read in the Seurat objects
dg0d_dm <- readRDS(file = "out/seurat_objects/dg0d_adt_norm.rds")
dg0h_dm <- readRDS(file = "out/seurat_objects/dg0h_adt_norm.rds")
v113_dm <- readRDS(file = "out/seurat_objects/v113_adt_norm.rds")
v015_dm <- readRDS(file = "out/seurat_objects/v015_adt_norm.rds")

###############################################
##                                           ##
## Cluster the cells based on their RNA data ##
##                                           ##
## Overlay ADT plots with ADT/RNA Expression ## 
##                                           ##
## Overlay RNA plots with ADt/RNA Expression ##
##                                           ##
###############################################

#DG0D
DefaultAssay(dg0d_dm) <- "RNA"
dg0d_dm <- NormalizeData(dg0d_dm, normalization.method = "LogNormalize", scale.factor = 10000)
dg0d_dm <- FindVariableFeatures(dg0d_dm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0d_dm), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0d_dm)
dg0d_dm <- ScaleData(dg0d_dm, features = all.genes)

dg0d_dm <- RunPCA(dg0d_dm, features = VariableFeatures(object = dg0d_dm))
print(dg0d_dm[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0d_dm, dims = 1:5, reduction = "pca")
DimPlot(dg0d_dm, reduction = "pca")
ElbowPlot(dg0d_dm)

dg0d_dm <- FindNeighbors(dg0d_dm, dims = 1:17)
dg0d_dm <- FindClusters(dg0d_dm, resolution = 0.5)
dg0d_dm <- RunUMAP(dg0d_dm, dims = 1:17)
DimPlot(dg0d_dm, reduction = "umap", label = TRUE)
dg0d_markers <- FindAllMarkers(dg0d_dm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0d_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0d_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5


avgexp <- AverageExpression(dg0d_dm, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))
ggsave(file = "out/images/dg0d_rna_hm.png", width = 8, height = 8)

DefaultAssay(dg0d_dm) <- "adt"
FeaturePlot(dg0d_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/dg0d_nonT.png", width = 15, height = 10)


DefaultAssay(dg0d_dm) <- "RNA"
FeaturePlot(dg0d_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/dg0d_nonT_RNA.png", width = 15, height = 10)



#DG0H
DefaultAssay(dg0h_dm) <- "RNA"
dg0h_dm <- NormalizeData(dg0h_dm, normalization.method = "LogNormalize", scale.factor = 10000)
dg0h_dm <- FindVariableFeatures(dg0h_dm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0h_dm), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0h_dm)
dg0h_dm <- ScaleData(dg0h_dm, features = all.genes)

dg0h_dm <- RunPCA(dg0h_dm, features = VariableFeatures(object = dg0h_dm))
print(dg0h_dm[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0h_dm, dims = 1:5, reduction = "pca")
DimPlot(dg0h_dm, reduction = "pca")
ElbowPlot(dg0h_dm)

dg0h_dm <- FindNeighbors(dg0h_dm, dims = 1:17)
dg0h_dm <- FindClusters(dg0h_dm, resolution = 0.5)
dg0h_dm <- RunUMAP(dg0h_dm, dims = 1:17)
DimPlot(dg0h_dm, reduction = "umap", label = TRUE)
dg0h_markers <- FindAllMarkers(dg0h_dm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0h_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0h_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5


avgexp <- AverageExpression(dg0h_dm, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))
ggsave(file = "out/images/dg0h_rna_hm.png", width = 8, height = 8)


DefaultAssay(dg0h_dm) <- "adt"
FeaturePlot(dg0h_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9")) 
ggsave(file = "out/images/dg0h_nonT.png", width = 15, height = 10)

DefaultAssay(dg0h_dm) <- "RNA"
FeaturePlot(dg0h_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9")) 
ggsave(file = "out/images/dg0h_nonT_RNA.png", width = 15, height = 10)


#V113
DefaultAssay(v113_dm) <- "RNA"
v113_dm <- NormalizeData(v113_dm, normalization.method = "LogNormalize", scale.factor = 10000)
v113_dm <- FindVariableFeatures(v113_dm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v113_dm), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v113_dm)
v113_dm <- ScaleData(v113_dm, features = all.genes)

v113_dm <- RunPCA(v113_dm, features = VariableFeatures(object = v113_dm))
print(v113_dm[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v113_dm, dims = 1:5, reduction = "pca")
DimPlot(v113_dm, reduction = "pca")
ElbowPlot(v113_dm)

v113_dm <- FindNeighbors(v113_dm, dims = 1:17)
v113_dm <- FindClusters(v113_dm, resolution = 0.5)
v113_dm <- RunUMAP(v113_dm, dims = 1:17)
DimPlot(v113_dm, reduction = "umap", label = TRUE)
v113_markers <- FindAllMarkers(v113_dm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v113_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v113_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
DoHeatmap(v113_dm, features = top5$gene) + NoLegend()


avgexp <- AverageExpression(v113_dm, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))
ggsave(file = "out/images/v113_rna_hm.png", width = 8, height = 8)


DefaultAssay(v113_dm) <- "adt"
FeaturePlot(v113_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/v113_nonT.png", width = 15, height = 10)


DefaultAssay(v113_dm) <- "RNA"
FeaturePlot(v113_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/v113_nonT_RNA.png", width = 15, height = 10)



#V015
DefaultAssay(v015_dm) <- "RNA"
v015_dm <- NormalizeData(v015_dm, normalization.method = "LogNormalize", scale.factor = 10000)
v015_dm <- FindVariableFeatures(v015_dm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v015_dm), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v015_dm)
v015_dm <- ScaleData(v015_dm, features = all.genes)

v015_dm <- RunPCA(v015_dm, features = VariableFeatures(object = v015_dm))
print(v015_dm[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v015_dm, dims = 1:5, reduction = "pca")
DimPlot(v015_dm, reduction = "pca")
ElbowPlot(v015_dm)

v015_dm <- FindNeighbors(v015_dm, dims = 1:17)
v015_dm <- FindClusters(v015_dm, resolution = 0.5)
v015_dm <- RunUMAP(v015_dm, dims = 1:17)
DimPlot(v015_dm, reduction = "umap", label = TRUE)
v015_markers <- FindAllMarkers(v015_dm, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v015_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v015_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
DoHeatmap(v015_dm, features = top5$gene) + NoLegend()


avgexp <- AverageExpression(v015_dm, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))
ggsave(file = "out/images/v015_rna_hm.png", width = 8, height = 8)


DefaultAssay(v015_dm) <- "adt"
FeaturePlot(v015_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/v015_nonT.png", width = 15, height = 10)


DefaultAssay(v015_dm) <- "RNA"
FeaturePlot(v015_dm, features = c("rna_MS4A1", "adt_CD20", "rna_CD14", "adt_CD14", "rna_CD3E", "rna_TRAC", "rna_TRDV1", "rna_TRDV3", "rna_TRDV2", "rna_TRGV9"))
ggsave(file = "out/images/v015_nonT_RNA.png", width = 15, height = 10)



## Save the DE genes and Seurat Objects ##
saveRDS(dg0d_dm, file = "out/seurat_objects/dg0d_clust.rds")
saveRDS(dg0h_dm, file = "out/seurat_objects/dg0h_clust.rds")
saveRDS(v113_dm, file = "out/seurat_objects/v113_clust.rds")
saveRDS(v015_dm, file = "out/seurat_objects/v015_clust.rds")


###########################################
## Exclude the non-T cells from the data ##
###########################################

#Read in the Seurat objects and DE genes
dg0d_dm <- readRDS(file = "out/seurat_objects/dg0d_clust.rds")
dg0h_dm <- readRDS(file = "out/seurat_objects/dg0h_clust.rds")
v113_dm <- readRDS(file = "out/seurat_objects/v113_clust.rds")
v015_dm <- readRDS(file = "out/seurat_objects/v015_clust.rds")

#Exclude the non-T cell clusters in each PTID
DefaultAssay(dg0d_dm) <- "adt"
Idents(dg0d_dm) <- "adt_snn_res.0.4"
dg0d_T <- dg0d_dm %>%
  subset(idents = c(0, 1, 2, 3, 4))
unique(dg0d_T$adt_snn_res.0.4)

DefaultAssay(dg0h_dm) <- "RNA"
Idents(dg0h_dm) <- "RNA_snn_res.0.5"
dg0h_T <- dg0h_dm %>%
  subset(idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
unique(dg0h_T$RNA_snn_res.0.5)

DefaultAssay(v113_dm) <- "RNA"
Idents(v113_dm) <- "RNA_snn_res.0.5"
v113_T <- v113_dm %>%
  subset(idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12))
unique(v113_T$RNA_snn_res.0.5)

DefaultAssay(v015_dm) <- "RNA"
Idents(v015_dm) <- "RNA_snn_res.0.5"
v015_T <- v015_dm %>%
  subset(idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13))
unique(v015_T$RNA_snn_res.0.5)

#Save the new T-only objects
saveRDS(dg0d_T, file = "out/seurat_objects/dg0d_T.rds")
saveRDS(dg0h_T, file = "out/seurat_objects/dg0h_T.rds")
saveRDS(v113_T, file = "out/seurat_objects/v113_T.rds")
saveRDS(v015_T, file = "out/seurat_objects/v015_T.rds")

dg0d <- readRDS("out/seurat_objects/dg0d_T.rds")
dg0h <- readRDS("out/seurat_objects/dg0h_T.rds")
v015 <- readRDS("out/seurat_objects/v015_T.rds")
v113 <- readRDS("out/seurat_objects/v113_T.rds")

qc_table <- read.csv('qc_table.csv')
qc_table$tcells <- c(ncol(dg0d), ncol(dg0h), ncol(v015), ncol(v113))

write.csv(qc_table, 'qc_table.csv')
