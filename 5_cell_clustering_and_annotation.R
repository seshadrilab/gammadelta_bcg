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
dg0d <- readRDS(file = "out/seurat_objects/dg0d_T.rds")
dg0h <- readRDS(file = "out/seurat_objects/dg0h_T.rds")
v113 <- readRDS(file = "out/seurat_objects/v113_T.rds")
v015 <- readRDS(file = "out/seurat_objects/v015_T.rds")

#Confirm that the non-T cells were removed
DefaultAssay(dg0d) <- "adt"
DimPlot(dg0d, reduction = "umap_adt", label = TRUE)

DefaultAssay(dg0h) <- "RNA"
DimPlot(dg0h, reduction = "umap", label = TRUE)

DefaultAssay(v113) <- "RNA"
DimPlot(v113, reduction = "umap", label = TRUE)

DefaultAssay(v015) <- "RNA"
DimPlot(v015, reduction = "umap", label = TRUE)

#########################################
## Re-cluster the new T-cell-only data ##
#########################################

## DG0D ##

# ADT
DefaultAssay(dg0d) <- "adt"
VariableFeatures(dg0d) <- rownames(dg0d[["adt"]])
dg0d <- ScaleData(dg0d)

dg0d_adt <- GetAssayData(dg0d, slot = "data", assay = "adt")
dg0d_adt_dist <- dist(t(dg0d_adt))
dg0d[["umap_adt"]] <- RunUMAP(dg0d_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
dg0d[["adt_snn"]] <- FindNeighbors(dg0d_adt_dist)$snn
dg0d <- FindClusters(dg0d, resolution = 0.4, graph.name = "adt_snn")
DimPlot(dg0d, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(dg0d, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# GEX
DefaultAssay(dg0d) <- "RNA"
dg0d <- NormalizeData(dg0d, normalization.method = "LogNormalize", scale.factor = 10000)
dg0d <- FindVariableFeatures(dg0d, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0d), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0d)
dg0d <- ScaleData(dg0d, features = all.genes)

dg0d <- RunPCA(dg0d, features = VariableFeatures(object = dg0d))
print(dg0d[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0d, dims = 1:5, reduction = "pca")
DimPlot(dg0d, reduction = "pca")
ElbowPlot(dg0d)

dg0d <- FindNeighbors(dg0d, dims = 1:17)
dg0d <- FindClusters(dg0d, resolution = 0.5)
dg0d <- RunUMAP(dg0d, dims = 1:17)
DimPlot(dg0d, reduction = "umap", label = TRUE)
dg0d_markers <- FindAllMarkers(dg0d, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0d_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0d_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(dg0d, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))

Idents(dg0d) <- "Tissue"
DimPlot(dg0d, reduction = "umap", label = TRUE)

Idents(dg0d) <- "Timepoint"
DimPlot(dg0d, reduction = "umap", label = TRUE)



## DG0H ##

# ADT
DefaultAssay(dg0h) <- "adt"
VariableFeatures(dg0h) <- rownames(dg0h[["adt"]])
dg0h <- ScaleData(dg0h)

dg0h_adt <- GetAssayData(dg0h, slot = "data", assay = "adt")
dg0h_adt_dist <- dist(t(dg0h_adt))
dg0h[["umap_adt"]] <- RunUMAP(dg0h_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
dg0h[["adt_snn"]] <- FindNeighbors(dg0h_adt_dist)$snn
dg0h <- FindClusters(dg0h, resolution = 0.4, graph.name = "adt_snn")
DimPlot(dg0h, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(dg0h, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# GEX
DefaultAssay(dg0h) <- "RNA"
dg0h <- NormalizeData(dg0h, normalization.method = "LogNormalize", scale.factor = 10000)
dg0h <- FindVariableFeatures(dg0h, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0h), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0h)
dg0h <- ScaleData(dg0h, features = all.genes)

dg0h <- RunPCA(dg0h, features = VariableFeatures(object = dg0h))
print(dg0h[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0h, dims = 1:5, reduction = "pca")
DimPlot(dg0h, reduction = "pca")
ElbowPlot(dg0h)

dg0h <- FindNeighbors(dg0h, dims = 1:17)
dg0h <- FindClusters(dg0h, resolution = 0.5)
dg0h <- RunUMAP(dg0h, dims = 1:17)
DimPlot(dg0h, reduction = "umap", label = TRUE)
dg0h_markers <- FindAllMarkers(dg0h, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0h_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0h_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(dg0h, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))

Idents(dg0h) <- "Tissue"
DimPlot(dg0h, reduction = "umap", label = TRUE)

Idents(dg0h) <- "Timepoint"
DimPlot(dg0h, reduction = "umap", label = TRUE)



## V113 ##

# ADT
DefaultAssay(v113) <- "adt"
VariableFeatures(v113) <- rownames(v113[["adt"]])
v113 <- ScaleData(v113)

v113_adt <- GetAssayData(v113, slot = "data", assay = "adt")
v113_adt_dist <- dist(t(v113_adt))
v113[["umap_adt"]] <- RunUMAP(v113_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
v113[["adt_snn"]] <- FindNeighbors(v113_adt_dist)$snn
v113 <- FindClusters(v113, resolution = 0.4, graph.name = "adt_snn")
DimPlot(v113, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(v113, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# GEX
DefaultAssay(v113) <- "RNA"
v113 <- NormalizeData(v113, normalization.method = "LogNormalize", scale.factor = 10000)
v113 <- FindVariableFeatures(v113, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v113), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v113)
v113 <- ScaleData(v113, features = all.genes)

v113 <- RunPCA(v113, features = VariableFeatures(object = v113))
print(v113[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v113, dims = 1:5, reduction = "pca")
DimPlot(v113, reduction = "pca")
ElbowPlot(v113)

v113 <- FindNeighbors(v113, dims = 1:17)
v113 <- FindClusters(v113, resolution = 0.5)
v113 <- RunUMAP(v113, dims = 1:17)
DimPlot(v113, reduction = "umap", label = TRUE)
v113_markers <- FindAllMarkers(v113, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v113_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v113_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(v113, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))

Idents(v113) <- "Tissue"
DimPlot(v113, reduction = "umap", label = TRUE)

Idents(v113) <- "Timepoint"
DimPlot(v113, reduction = "umap", label = TRUE)



## V015 ##

# ADT
DefaultAssay(v015) <- "adt"
VariableFeatures(v015) <- rownames(v015[["adt"]])
v015 <- ScaleData(v015)

v015_adt <- GetAssayData(v015, slot = "data", assay = "adt")
v015_adt_dist <- dist(t(v015_adt))
v015[["umap_adt"]] <- RunUMAP(v015_adt_dist, assay = "adt", reduction.key = "adtUMAP_")
v015[["adt_snn"]] <- FindNeighbors(v015_adt_dist)$snn
v015 <- FindClusters(v015, resolution = 0.4, graph.name = "adt_snn")
DimPlot(v015, reduction = "umap_adt", label = TRUE)

avgexp <- AverageExpression(v015, assays = "adt", return.seurat = T)
hm_order <- c("CD20", "CD14", "CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


# GEX
DefaultAssay(v015) <- "RNA"
v015 <- NormalizeData(v015, normalization.method = "LogNormalize", scale.factor = 10000)
v015 <- FindVariableFeatures(v015, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v015), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v015)
v015 <- ScaleData(v015, features = all.genes)

v015 <- RunPCA(v015, features = VariableFeatures(object = v015))
print(v015[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v015, dims = 1:5, reduction = "pca")
DimPlot(v015, reduction = "pca")
ElbowPlot(v015)

v015 <- FindNeighbors(v015, dims = 1:17)
v015 <- FindClusters(v015, resolution = 0.5)
v015 <- RunUMAP(v015, dims = 1:17)
DimPlot(v015, reduction = "umap", label = TRUE)
v015_markers <- FindAllMarkers(v015, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v015_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v015_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(v015, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TRAC"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))

Idents(v015) <- "Tissue"
DimPlot(v015, reduction = "umap", label = TRUE)

Idents(v015) <- "Timepoint"
DimPlot(v015, reduction = "umap", label = TRUE)


## Save the Outputs ##

saveRDS(dg0d, file = "out/dg0d_analysis.rds")
saveRDS(dg0h, file = "out/dg0h_analysis.rds")
saveRDS(v113, file = "out/v113_analysis.rds")
saveRDS(v015, file = "out/v015_analysis.rds")

write.csv(dg0d_markers, file = "out/dg0d_markers_Tcell.csv")
write.csv(dg0h_markers, file = "out/dg0h_markers_Tcell.csv")
write.csv(v113_markers, file = "out/v113_markers_Tcell.csv")
write.csv(v015_markers, file = "out/v015_markers_Tcell.csv")




##################################################
## Add a gamma-delta annotation to the metadata ##
##################################################

# Read in the files
dg0d <- readRDS(file = "out/seurat_objects/dg0d_T.rds")
dg0h <- readRDS(file = "out/seurat_objects/dg0h_T.rds")
v113 <- readRDS(file = "out/seurat_objects/v113_T.rds")
v015 <- readRDS(file = "out/seurat_objects/v015_T.rds")

dg0d_gd_bc <- read.csv("out/tcr/dg0d_gd_bc.csv")
dg0h_gd_bc <- read.csv("out/tcr/dg0h_gd_bc.csv")
v113_gd_bc <- read.csv("out/tcr/v113_gd_bc.csv")
v015_gd_bc <- read.csv("out/tcr/v015_gd_bc.csv")

# Annotate the Seurat objects
dg0d_gd_bc <- dg0d_gd_bc$barcode
dg0d_gd_bc <- paste(dg0d_gd_bc, collapse = "|")
dg0d_gd_bc <- str_detect(colnames(dg0d), dg0d_gd_bc)
dg0d@meta.data$Celltype[dg0d_gd_bc] <- "gd"
dg0d@meta.data$Celltype[!dg0d_gd_bc] <- "other"

dg0h_gd_bc <- dg0h_gd_bc$barcode
dg0h_gd_bc <- paste(dg0h_gd_bc, collapse = "|")
dg0h_gd_bc <- str_detect(colnames(dg0h), dg0h_gd_bc)
dg0h@meta.data$Celltype[dg0h_gd_bc] <- "gd"
dg0h@meta.data$Celltype[!dg0h_gd_bc] <- "other"

v113_gd_bc <- v113_gd_bc$barcode
v113_gd_bc <- paste(v113_gd_bc, collapse = "|")
v113_gd_bc <- str_detect(colnames(v113), v113_gd_bc)
v113@meta.data$Celltype[v113_gd_bc] <- "gd"
v113@meta.data$Celltype[!v113_gd_bc] <- "other"

v015_gd_bc <- v015_gd_bc$barcode
v015_gd_bc <- paste(v015_gd_bc, collapse = "|")
v015_gd_bc <- str_detect(colnames(v015), v015_gd_bc)
v015@meta.data$Celltype[v015_gd_bc] <- "gd"
v015@meta.data$Celltype[!v015_gd_bc] <- "other"

#Save the updated RDS files
saveRDS(dg0d, file = "out/seurat_objects/dg0d_analysis.rds")
saveRDS(dg0h, file = "out/seurat_objects/dg0h_analysis.rds")
saveRDS(v113, file = "out/seurat_objects/v113_analysis.rds")
saveRDS(v015, file = "out/seurat_objects/v015_analysis.rds")

Idents(dg0d) <- "Celltype"
dg0d_gd <- dg0d %>%
  subset(idents = c("gd"))

Idents(dg0h) <- "Celltype"
dg0h_gd <- dg0h %>%
  subset(idents = c("gd"))

Idents(v015) <- "Celltype"
v015_gd <- v015 %>%
  subset(idents = c("gd"))

Idents(v113) <- "Celltype"
v113_gd <- v113 %>%
  subset(idents = c("gd"))


qc_table <- read.csv('qc_table.csv')
qc_table$gd_tcells <- c(ncol(dg0d_gd), ncol(dg0h_gd), ncol(v015_gd), ncol(v113_gd))

## Visualize the final cell counts at each QC step in stacked bar plot format
qc_table <- qc_table[,3:9]
qc_table <- qc_table %>% pivot_longer(cols = 2:7, names_to = 'step', values_to = 'count')

qc_table %>%
  ggplot(aes(x = step, y = count, fill = ptids)) +
  geom_col(position = 'stack') +
  scale_x_discrete(limits = c('initial_count', 'hto_minimum', 'gex_filters', 'hto_singlets', 'tcells', 'gd_tcells')) +
  scale_fill_manual(values = c("#463480FF", "#2C738EFF", "#1F998AFF", 'lightgrey')) + theme_bw()
ggsave('ivbcg_qc.pdf')

####################################################
## Re-cluster the gamma delta T cells only -- GEX ##
####################################################

## DG0D ##
Idents(dg0d) <- "Celltype"
dg0d_gd <- dg0d %>%
  subset(idents = c("gd"))

DefaultAssay(dg0d_gd) <- "RNA"
dg0d_gd <- NormalizeData(dg0d_gd, normalization.method = "LogNormalize", scale.factor = 10000)
dg0d_gd <- FindVariableFeatures(dg0d_gd, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0d_gd), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0d_gd)
dg0d_gd <- ScaleData(dg0d_gd, features = all.genes)

dg0d_gd <- RunPCA(dg0d_gd, features = VariableFeatures(object = dg0d_gd))
print(dg0d_gd[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0d_gd, dims = 1:5, reduction = "pca")
DimPlot(dg0d_gd, reduction = "pca")
ElbowPlot(dg0d_gd)

dg0d_gd <- FindNeighbors(dg0d_gd, dims = 1:17)
dg0d_gd <- FindClusters(dg0d_gd, resolution = 0.5)
dg0d_gd <- RunUMAP(dg0d_gd, dims = 1:17)
DimPlot(dg0d_gd, reduction = "umap", label = TRUE)
dg0d_gd_markers <- FindAllMarkers(dg0d_gd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0d_gd_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0d_gd_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(dg0d_gd, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))



## DG0H ##
Idents(dg0h) <- "Celltype"
dg0h_gd <- dg0h %>%
  subset(idents = c("gd"))

DefaultAssay(dg0h_gd) <- "RNA"
dg0h_gd <- NormalizeData(dg0h_gd, normalization.method = "LogNormalize", scale.factor = 10000)
dg0h_gd <- FindVariableFeatures(dg0h_gd, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(dg0h_gd), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(dg0h_gd)
dg0h_gd <- ScaleData(dg0h_gd, features = all.genes)

dg0h_gd <- RunPCA(dg0h_gd, features = VariableFeatures(object = dg0h_gd))
print(dg0h_gd[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(dg0h_gd, dims = 1:5, reduction = "pca")
DimPlot(dg0h_gd, reduction = "pca")
ElbowPlot(dg0h_gd)

dg0h_gd <- FindNeighbors(dg0h_gd, dims = 1:17)
dg0h_gd <- FindClusters(dg0h_gd, resolution = 0.5)
dg0h_gd <- RunUMAP(dg0h_gd, dims = 1:17)
DimPlot(dg0h_gd, reduction = "umap", label = TRUE)
dg0h_gd_markers <- FindAllMarkers(dg0h_gd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
dg0h_gd_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
dg0h_gd_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(dg0h_gd, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))


## V113 ##
Idents(v113) <- "Celltype"
v113_gd <- v113 %>%
  subset(idents = c("gd"))

DefaultAssay(v113_gd) <- "RNA"
v113_gd <- NormalizeData(v113_gd, normalization.method = "LogNormalize", scale.factor = 10000)
v113_gd <- FindVariableFeatures(v113_gd, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v113_gd), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v113_gd)
v113_gd <- ScaleData(v113_gd, features = all.genes)

v113_gd <- RunPCA(v113_gd, features = VariableFeatures(object = v113_gd))
print(v113_gd[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v113_gd, dims = 1:5, reduction = "pca")
DimPlot(v113_gd, reduction = "pca")
ElbowPlot(v113_gd)

v113_gd <- FindNeighbors(v113_gd, dims = 1:17)
v113_gd <- FindClusters(v113_gd, resolution = 0.5)
v113_gd <- RunUMAP(v113_gd, dims = 1:17)
DimPlot(v113_gd, reduction = "umap", label = TRUE)
v113_gd_markers <- FindAllMarkers(v113_gd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v113_gd_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v113_gd_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(v113_gd, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))


## V015 ##
Idents(v015) <- "Celltype"
v015_gd <- v015 %>%
  subset(idents = c("gd"))

DefaultAssay(v015_gd) <- "RNA"
v015_gd <- NormalizeData(v015_gd, normalization.method = "LogNormalize", scale.factor = 10000)
v015_gd <- FindVariableFeatures(v015_gd, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(v015_gd), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(v015_gd)
v015_gd <- ScaleData(v015_gd, features = all.genes)

v015_gd <- RunPCA(v015_gd, features = VariableFeatures(object = v015_gd))
print(v015_gd[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(v015_gd, dims = 1:5, reduction = "pca")
DimPlot(v015_gd, reduction = "pca")
ElbowPlot(v015_gd)

v015_gd <- FindNeighbors(v015_gd, dims = 1:17)
v015_gd <- FindClusters(v015_gd, resolution = 0.5)
v015_gd <- RunUMAP(v015_gd, dims = 1:17)
DimPlot(v015_gd, reduction = "umap", label = TRUE)
v015_gd_markers <- FindAllMarkers(v015_gd, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #This step can take several minutes.
v015_gd_markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
v015_gd_markers %>%
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(v015_gd, assays = "RNA", return.seurat = T)
avgexp <- ScaleData(avgexp, features = all.genes)
DoHeatmap(avgexp, features = c(top5$gene, "TRDV1", "TRDV3", "TRDV2", "TRGV9"), raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))



################################################################################################################################################
# TEST ANALYSIS 
################################################################################################################################################

############################################
## Reorder DE genes according to keywords ##
############################################

## V113 ##
gene_annot <- read.csv(file = "D:/Shared drives/Seshadri Lab/Lab Members/Maerz_Megan/Aim_2/Gene_annotations.csv") # Import the gene annotation csv.
gene_annot <- gene_annot %>% rename(gene = gene_symbol)
gene_annot_kw <- gene_annot %>% subset(select = c("gene", "keywords")) # Subset to only list the corresponding keyword for each gene.

top10 <- v113_markers %>% # Generate a list of the top 10 DE genes in each cluster.
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC)

top10.2 <- merge(top10, gene_annot_kw, by = "gene") # Annotate and sort the list of top 10 DE genes.
top10.2 <- top10.2 %>% arrange(keywords)

Idents(v113) <- "RNA_snn_res.0.5"
avgexp <- AverageExpression(v113, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
DoHeatmap(avgexp, features = top10.2$gene, raster = FALSE, draw.lines = F) + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Visualize the heat map.
ggsave(device = NULL, "D:/Shared drives/Seshadri Lab/Lab Members/Maerz_Megan/Aim_2/R_files/GEX/Out/heatmap_top10_sorted.png", width = 5, height = 8, units = c("in")) # Save an image of the heat map.

##################################################################
## Visualize distribution changes across tissues and timepoints ##
##################################################################

dg0d <- readRDS(file = "out/dg0d_analysis.rds")
dg0h <- readRDS(file = "out/dg0h_analysis.rds")
v113 <- readRDS(file = "out/v113_analysis.rds")
v015 <- readRDS(file = "out/v015_analysis.rds")

dg0d_markers <- read.csv(file = "out/dg0d_markers_Tcell.csv")
dg0h_markers <- read.csv(file = "out/dg0h_markers_Tcell.csv")
v113_markers <- read.csv(file = "out/v113_markers_Tcell.csv")
v015_markers <- read.csv(file = "out/v015_markers_Tcell.csv")

## DG0D ##
Idents(dg0d) <- "Tissue"
dg0d_p <- dg0d %>% subset(idents = c("PBMC"))
Idents(dg0d_p) <- "Timepoint"
DimPlot(dg0d_p, reduction = "umap", label = TRUE)


## V113 ##
Idents(v113) <- "Tissue"
v113_p <- v113 %>% subset(idents = c("PBMC"))
Idents(v113_p) <- "Timepoint"
DimPlot(v113_p, reduction = "umap", label = TRUE)

Idents(v113) <- "Timepoint"
v113_4 <- v113 %>% subset(idents = c("Week_4"))
Idents(v113_4) <- "Tissue"
DimPlot(v113_4, reduction = "umap", label = TRUE)


###############################################
## DE markers across Tissues and Time Points ##
###############################################

dg0d <- readRDS(file = "out/dg0d_analysis.rds")
dg0h <- readRDS(file = "out/dg0h_analysis.rds")
v113 <- readRDS(file = "out/v113_analysis.rds")
v015 <- readRDS(file = "out/v015_analysis.rds")

## PBMC Time Points, RNA ##

#DG0D
Idents(dg0d) <- "Tissue"
dg0d_pbmc <- dg0d %>% subset(idents = c("PBMC"))
Idents(dg0d_pbmc) <- "Timepoint"

markers <- FindAllMarkers(dg0d_pbmc)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20

avgexp <- AverageExpression(dg0d_pbmc, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


#DG0H
Idents(dg0h) <- "Tissue"
dg0h_pbmc <- dg0h %>% subset(idents = c("PBMC"))
Idents(dg0h_pbmc) <- "Timepoint"

markers <- FindAllMarkers(dg0h_pbmc)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20

avgexp <- AverageExpression(dg0h_pbmc, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


#V113
Idents(v113) <- "Tissue"
v113_pbmc <- v113 %>% subset(idents = c("PBMC"))
Idents(v113_pbmc) <- "Timepoint"
v113_pbmc$Timepoint <- v113_pbmc$Timepoint %>% factor(levels = c("Week_0, Week_4, Week_8"))

markers <- FindAllMarkers(v113_pbmc)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20

avgexp <- AverageExpression(v113_pbmc, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


#V015
Idents(v015) <- "Tissue"
v015_pbmc <- v015 %>% subset(idents = c("PBMC"))
Idents(v015_pbmc) <- "Timepoint"
v015_pbmc$Timepoint <- v015_pbmc$Timepoint %>% factor(levels = c("Week_0, Week_4, Week_8"))

markers <- FindAllMarkers(v015_pbmc)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20

avgexp <- AverageExpression(v015_pbmc, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


## PBMC Time Points, ADT ##

#DG0D
Idents(dg0d_pbmc) <- "Timepoint"
avgexp <- AverageExpression(dg0d_pbmc, assays = "adt", return.seurat = T)
avgexp <- ScaleData(avgexp, features = rownames(dg0d_pbmc$nCount_adt))
hm_order <- c("CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 14)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap

#DG0H
Idents(dg0h_pbmc) <- "Timepoint"
avgexp <- AverageExpression(dg0h_pbmc, assays = "adt", return.seurat = T)
avgexp <- ScaleData(avgexp, features = rownames(dg0h_pbmc$nCount_adt))
hm_order <- c("CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 14)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap

#V113
Idents(v113_pbmc) <- "Timepoint"
avgexp <- AverageExpression(v113_pbmc, assays = "adt", return.seurat = T)
avgexp <- ScaleData(avgexp, features = rownames(v113_pbmc$nCount_adt))
hm_order <- c("CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 14)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap

#V015
Idents(v015_pbmc) <- "Timepoint"
avgexp <- AverageExpression(v015_pbmc, assays = "adt", return.seurat = T)
avgexp <- ScaleData(avgexp, features = rownames(v015_pbmc$nCount_adt))
hm_order <- c("CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "CD11b", "CD163", "CD11c", "CD161", "CD16", "CD86", "CCR6", "CD95", "CXCR3", "PD-L1")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 14)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


###############################################################
## List the top 10 genes that haven't already been annotated ##
###############################################################

top10.1 <- dg0d_markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC)

top10.2 <- dg0h_markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC)

top10.3 <- v113_markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC)

top10.4 <- v015_markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt=avg_log2FC)

top10.all <- do.call("rbind", list(top10.1, top10.2, top10.3, top10.4))
mamu_genes <- top10.all$gene %>% as.data.frame() %>% distinct() %>% rename("gene_symbol" = ".")
mamu_genes$source <- "new"

gene_annot <- read.csv(file = "D:/Shared drives/Seshadri Lab/Lab Members/Maerz_Megan/Aim_2/Gene_annotations.csv") # Import the gene annotation csv.
gene_annot <- gene_annot %>% rename(gene = gene_symbol)
gene_annot_kw <- gene_annot %>% subset(select = c("gene")) %>% rename("gene_symbol" = "gene")
gene_annot_kw$source <- "anno"

new <- mamu_genes %>% anti_join(gene_annot_kw, by = "gene_symbol")
write.csv(new, file = "D:/Shared drives/Seshadri Lab/Lab Members/Maerz_Megan/Aim_2/genes_needing_anno.csv")



#######################################################
## Test Differentiation of AB T cells vs. GD T cells ##
#######################################################

## DG0D ##
dg0d_Tsubset <- rownames(dg0d@assays$RNA@data)
TCR_genes <- grepl("^TRAV|^TRAD|^TRAJ|^TRAC|^TRBV|^TRBJ|^TRBC|^TRDV|^TRDD|^TRDJ|^TRDC|^TRGV|^TRGJ|^TRGC", dg0d_Tsubset)
dg0d_Tsubset <- dg0d_Tsubset[TCR_genes]
dg0d_Tgenes <- dg0d@assays$RNA@data[dg0d_Tsubset,]
dg0d[["Tgex"]] <- CreateAssayObject(data = dg0d_Tgenes)

DefaultAssay(dg0d) <- "Tgex"
VariableFeatures(dg0d) <- rownames(dg0d[["Tgex"]])
dg0d <- ScaleData(dg0d)
dg0d_Tgex <- GetAssayData(dg0d, slot = "data", assay = "Tgex") %>% as.matrix()
dist <- dist(t(dg0d_Tgex))
dg0d[["umap_tcr"]] <- RunUMAP(dist, assay = "Tgex", reduction.key = "tcrUMAP_")
dg0d[["tcr_snn"]] <- FindNeighbors(dist)$snn
dg0d <- FindClusters(dg0d, resolution = 0.4, graph.name = "tcr_snn")
DimPlot(dg0d, reduction = "umap_tcr", label = TRUE)
markers <- FindAllMarkers(dg0d)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20
avgexp <- AverageExpression(dg0d, assays = "Tgex", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20

## DG0H ##
dg0h_Tsubset <- rownames(dg0h@assays$RNA@data)
TCR_genes <- grepl("^TRAV|^TRAD|^TRAJ|^TRAC|^TRBV|^TRBJ|^TRBC|^TRDV|^TRDD|^TRDJ|^TRDC|^TRGV|^TRGJ|^TRGC", dg0h_Tsubset)
dg0h_Tsubset <- dg0h_Tsubset[TCR_genes]
dg0h_Tgenes <- dg0h@assays$RNA@data[dg0h_Tsubset,]
dg0h[["Tgex"]] <- CreateAssayObject(data = dg0h_Tgenes)

DefaultAssay(dg0h) <- "Tgex"
VariableFeatures(dg0h) <- rownames(dg0h[["Tgex"]])
dg0h <- ScaleData(dg0h)
dg0h_Tgex <- GetAssayData(dg0h, slot = "data", assay = "Tgex") %>% as.matrix()
dist <- dist(t(dg0h_Tgex))
dg0h[["umap_tcr"]] <- RunUMAP(dist, assay = "Tgex", reduction.key = "tcrUMAP_")
dg0h[["tcr_snn"]] <- FindNeighbors(dist)$snn
dg0h <- FindClusters(dg0h, resolution = 0.4, graph.name = "tcr_snn")
DimPlot(dg0h, reduction = "umap_tcr", label = TRUE)
markers <- FindAllMarkers(dg0h)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20
avgexp <- AverageExpression(dg0h, assays = "Tgex", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


## V113 ##
v113_Tsubset <- rownames(v113@assays$RNA@data)
TCR_genes <- grepl("^TRAV|^TRAD|^TRAJ|^TRAC|^TRBV|^TRBJ|^TRBC|^TRDV|^TRDD|^TRDJ|^TRDC|^TRGV|^TRGJ|^TRGC", v113_Tsubset)
v113_Tsubset <- v113_Tsubset[TCR_genes]
v113_Tgenes <- v113@assays$RNA@data[v113_Tsubset,]
v113[["Tgex"]] <- CreateAssayObject(data = v113_Tgenes)

DefaultAssay(v113) <- "Tgex"
VariableFeatures(v113) <- rownames(v113[["Tgex"]])
v113 <- ScaleData(v113)
v113_Tgex <- GetAssayData(v113, slot = "data", assay = "Tgex") %>% as.matrix()
dist <- dist(t(v113_Tgex))
v113[["umap_tcr"]] <- RunUMAP(dist, assay = "Tgex", reduction.key = "tcrUMAP_")
v113[["tcr_snn"]] <- FindNeighbors(dist)$snn
v113 <- FindClusters(v113, resolution = 0.4, graph.name = "tcr_snn")
DimPlot(v113, reduction = "umap_tcr", label = TRUE)
markers <- FindAllMarkers(v113)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20
avgexp <- AverageExpression(v113, assays = "Tgex", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20


## V015 ##
v015_Tsubset <- rownames(v015@assays$RNA@data)
TCR_genes <- grepl("^TRAV|^TRAD|^TRAJ|^TRAC|^TRBV|^TRBJ|^TRBC|^TRDV|^TRDD|^TRDJ|^TRDC|^TRGV|^TRGJ|^TRGC", v015_Tsubset)
v015_Tsubset <- v015_Tsubset[TCR_genes]
v015_Tgenes <- v015@assays$RNA@data[v015_Tsubset,]
v015[["Tgex"]] <- CreateAssayObject(data = v015_Tgenes)

DefaultAssay(v015) <- "Tgex"
VariableFeatures(v015) <- rownames(v015[["Tgex"]])
v015 <- ScaleData(v015)
v015_Tgex <- GetAssayData(v015, slot = "data", assay = "Tgex") %>% as.matrix()
dist <- dist(t(v015_Tgex))
v015[["umap_tcr"]] <- RunUMAP(dist, assay = "Tgex", reduction.key = "tcrUMAP_")
v015[["tcr_snn"]] <- FindNeighbors(dist)$snn
v015 <- FindClusters(v015, resolution = 0.4, graph.name = "tcr_snn")
DimPlot(v015, reduction = "umap_tcr", label = TRUE)
markers <- FindAllMarkers(v015)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20
avgexp <- AverageExpression(v015, assays = "Tgex", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top20 <- DoHeatmap(avgexp, features = top20$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 12)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top20

