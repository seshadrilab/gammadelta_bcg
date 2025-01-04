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
dg0d_dm <- readRDS(file = "out/seurat_objects/dg0d_demux.rds")
dg0h_dm <- readRDS(file = "out/seurat_objects/dg0h_demux.rds")
v113_dm <- readRDS(file = "out/seurat_objects/v113_demux.rds")
v015_dm <- readRDS(file = "out/seurat_objects/v015_demux.rds")

#Clean up the ADT names
DefaultAssay(dg0d_dm) <- "adt"
DefaultAssay(dg0h_dm) <- "adt"
DefaultAssay(v113_dm) <- "adt"
DefaultAssay(v015_dm) <- "adt"

dg0d_adt <- as.data.frame(dg0d_dm@assays$adt@counts)
rownames(dg0d_adt) <- rownames(dg0d_adt) %>% str_replace(fixed(".1"), "")
dg0d_dm[["adt"]] <- CreateAssayObject(counts = dg0d_adt)

dg0h_adt <- as.data.frame(dg0h_dm@assays$adt@counts)
rownames(dg0h_adt) <- rownames(dg0h_adt) %>% str_replace(fixed(".1"), "")
dg0h_dm[["adt"]] <- CreateAssayObject(counts = dg0h_adt)

v113_adt <- as.data.frame(v113_dm@assays$adt@counts)
rownames(v113_adt) <- rownames(v113_adt) %>% str_replace(fixed(".1"), "")
v113_dm[["adt"]] <- CreateAssayObject(counts = v113_adt)

v015_adt <- as.data.frame(v015_dm@assays$adt@counts)
rownames(v015_adt) <- rownames(v015_adt) %>% str_replace(fixed(".1"), "")
v015_dm[["adt"]] <- CreateAssayObject(counts = v015_adt)


#Normalize the data and visualize the signal from each ADT
dg0d_dm <- NormalizeData(dg0d_dm, normalization.method = "CLR", margin = 2, assay = "adt")
dg0h_dm <- NormalizeData(dg0h_dm, normalization.method = "CLR", margin = 2, assay = "adt")
v113_dm <- NormalizeData(v113_dm, normalization.method = "CLR", margin = 2, assay = "adt")
v015_dm <- NormalizeData(v015_dm, normalization.method = "CLR", margin = 2, assay = "adt")

RidgePlot(dg0d_dm, assay = "adt", features = rownames(dg0d_dm), ncol = 4)
ggsave(file = "out/images/dg0d_adt_ridge.png", width = 15, height = 15)

RidgePlot(dg0h_dm, assay = "adt", features = rownames(dg0h_dm), ncol = 4)
ggsave(file = "out/images/dg0h_adt_ridge.png", width = 15, height = 15)

RidgePlot(v113_dm, assay = "adt", features = rownames(v113_dm), ncol = 4)
ggsave(file = "out/images/v113_adt_ridge.png", width = 15, height = 15)

RidgePlot(v015_dm, assay = "adt", features = rownames(v015_dm), ncol = 4)
ggsave(file = "out/images/v015_adt_ridge.png", width = 15, height = 15)


##################
## ADT Analysis ##
##################

## Cluster the cells based on their ADT data ##

# DGOD
DefaultAssay(dg0d_dm) <- "adt"
VariableFeatures(dg0d_dm) <- rownames(dg0d_dm[["adt"]])
dg0d_dm <- ScaleData(dg0d_dm)

dg0d_dm_adt <- GetAssayData(dg0d_dm, slot = "data", assay = "adt")
dg0d_adt_dist <- dist(t(dg0d_dm_adt))
dg0d_dm[["umap_adt"]] <- RunUMAP(dg0d_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
dg0d_dm[["adt_snn"]] <- FindNeighbors(dg0d_adt_dist)$snn
dg0d_dm <- FindClusters(dg0d_dm, resolution = 0.4, graph.name = "adt_snn")
DimPlot(dg0d_dm, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(dg0d_dm, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# DGOH
DefaultAssay(dg0h_dm) <- "adt"
VariableFeatures(dg0h_dm) <- rownames(dg0h_dm[["adt"]])
dg0h_dm <- ScaleData(dg0h_dm)

dg0h_dm_adt <- GetAssayData(dg0h_dm, slot = "data", assay = "adt")
dg0h_adt_dist <- dist(t(dg0h_dm_adt))
dg0h_dm[["umap_adt"]] <- RunUMAP(dg0h_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
dg0h_dm[["adt_snn"]] <- FindNeighbors(dg0h_adt_dist)$snn
dg0h_dm <- FindClusters(dg0h_dm, resolution = 0.4, graph.name = "adt_snn")
DimPlot(dg0h_dm, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(dg0h_dm, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# V113
DefaultAssay(v113_dm) <- "adt"
VariableFeatures(v113_dm) <- rownames(v113_dm[["adt"]])
v113_dm <- ScaleData(v113_dm)

v113_dm_adt <- GetAssayData(v113_dm, slot = "data", assay = "adt")
v113_adt_dist <- dist(t(v113_dm_adt))
v113_dm[["umap_adt"]] <- RunUMAP(v113_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
v113_dm[["adt_snn"]] <- FindNeighbors(v113_adt_dist)$snn
v113_dm <- FindClusters(v113_dm, resolution = 0.4, graph.name = "adt_snn")
DimPlot(v113_dm, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(v113_dm, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# V015
DefaultAssay(v015_dm) <- "adt"
VariableFeatures(v015_dm) <- rownames(v015_dm[["adt"]])
v015_dm <- ScaleData(v015_dm)

v015_dm_adt <- GetAssayData(v015_dm, slot = "data", assay = "adt")
v015_adt_dist <- dist(t(v015_dm_adt))
v015_dm[["umap_adt"]] <- RunUMAP(v015_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
v015_dm[["adt_snn"]] <- FindNeighbors(v015_adt_dist)$snn
v015_dm <- FindClusters(v015_dm, resolution = 0.4, graph.name = "adt_snn")
DimPlot(v015_dm, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(v015_dm, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# Save the updated Seurat objects
saveRDS(dg0d_dm, file = "out/seurat_objects/dg0d_adt_norm.rds")
saveRDS(dg0h_dm, file = "out/seurat_objects/dg0h_adt_norm.rds")
saveRDS(v113_dm, file = "out/seurat_objects/v113_adt_norm.rds")
saveRDS(v015_dm, file = "out/seurat_objects/v015_adt_norm.rds")

