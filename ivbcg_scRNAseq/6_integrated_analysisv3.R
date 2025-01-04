#Load the necessary libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(stringr)
library(ggplot2)
library(cowplot)
library(hdf5r)
library(Rgb)
library(reshape2)
library(tidyverse)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(viridis)

## Visualize the final cell counts at each QC step in stacked bar plot format
qc <- read.csv('qc_cell_counts.csv')
qc <- qc %>% pivot_longer(cols = 2:7, names_to = 'step', values_to = 'count')

qc %>%
  ggplot(aes(x = step, y = count, fill = ptid)) +
  geom_col(position = 'stack') +
  scale_x_discrete(limits = c('initial_count', 'hto_minimum', 'gex_filters', 'hto_singlets', 'tcells', 'gd_tcells')) +
  scale_fill_manual(values = c("#463480FF", "#2C738EFF", "#1F998AFF", 'lightgrey')) + theme_bw()

#Read in the Seurat objects
dg0d <- readRDS(file = "out/seurat_objects/dg0d_analysis.rds")
dg0h <- readRDS(file = "out/seurat_objects/dg0h_analysis.rds")
v113 <- readRDS(file = "out/seurat_objects/v113_analysis.rds")
v015 <- readRDS(file = "out/seurat_objects/v015_analysis.rds")

#Merge the data into a single Seurat object
ivbcg <- merge(dg0d, y = c(dg0h, v015, v113), add.cell.ids = c("dg0d", "dg0h", "v015", "v113"), merge.data = FALSE, merge.dr = NULL, project = "ivbcg")

#Generate a new merged object including only the gamma-delta T cells
Idents(ivbcg) <- "Celltype"
ivbcg_gd <- ivbcg %>%
  subset(idents = "gd")

#Cluster the new merged data using an RNA-only approach
DefaultAssay(ivbcg_gd) <- "RNA"
ivbcg_gd <- NormalizeData(ivbcg_gd, normalization.method = "LogNormalize", scale.factor = 10000)
ivbcg_gd <- FindVariableFeatures(ivbcg_gd, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(ivbcg_gd), 10) # Identify the 10 most highly variable genes.
all.genes <- rownames(ivbcg_gd)
ivbcg_gd <- ScaleData(ivbcg_gd, features = all.genes)

ivbcg_gd <- RunPCA(ivbcg_gd, features = VariableFeatures(object = ivbcg_gd))
print(ivbcg_gd[["pca"]], dims = 1:5, nfeatures = 3)
VizDimLoadings(ivbcg_gd, dims = 1:5, reduction = "pca")
DimPlot(ivbcg_gd, reduction = "pca")
ElbowPlot(ivbcg_gd)

ivbcg_gd <- FindNeighbors(ivbcg_gd, dims = 1:17)
ivbcg_gd <- FindClusters(ivbcg_gd, resolution = 0.5)
ivbcg_gd <- RunUMAP(ivbcg_gd, dims = 1:17)
DimPlot(ivbcg_gd, reduction = "umap", label = TRUE)

ivbcg_gd <- FindNeighbors(ivbcg_gd, dims = 1:14)
ivbcg_gd <- FindClusters(ivbcg_gd, resolution = 0.5)
ivbcg_gd <- RunUMAP(ivbcg_gd, dims = 1:14)
DimPlot(ivbcg_gd, reduction = "umap", label = TRUE)

ivbcg_gd <- FindNeighbors(ivbcg_gd, dims = 1:12)
ivbcg_gd <- FindClusters(ivbcg_gd, resolution = 0.5)
ivbcg_gd <- RunUMAP(ivbcg_gd, dims = 1:12)
DimPlot(ivbcg_gd, reduction = "umap", label = TRUE)

ivbcg_gd <- FindNeighbors(ivbcg_gd, dims = 1:9)
ivbcg_gd <- FindClusters(ivbcg_gd, resolution = 0.5)
ivbcg_gd <- RunUMAP(ivbcg_gd, dims = 1:9)
DimPlot(ivbcg_gd, reduction = "umap", label = TRUE)

ivbcg_gd <- FindNeighbors(ivbcg_gd, dims = 1:7)
ivbcg_gd <- FindClusters(ivbcg_gd, resolution = 0.5)
ivbcg_gd <- RunUMAP(ivbcg_gd, dims = 1:7)
DimPlot(ivbcg_gd, reduction = "umap", label = TRUE)


DimPlot(ivbcg_gd, reduction = "umap", group.by = "PTID")
DimPlot(ivbcg_gd, reduction = "umap", group.by = "Tissue")
DimPlot(ivbcg_gd, reduction = "umap", group.by = "Timepoint")

FeaturePlot(ivbcg_gd, pt.size = 0.5, features = c("TRDV1", 
                                                "TRDV2", 
                                                "TRDV3",
                                                "TRGV9"))


Idents(ivbcg_gd) <- "Tissue"
ivbcg_test <- ivbcg_gd %>%
  subset(idents = "PBMC")

Idents(ivbcg_test) <- "Timepoint"
ivbcg_test1 <- ivbcg_test %>%
  subset(idents = "Week_0")
ivbcg_test2 <- ivbcg_test %>%
  subset(idents = "Week_4")

Idents(ivbcg_test1) <- "PTID"
markers <- FindAllMarkers(ivbcg_test1)
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
avgexp <- AverageExpression(ivbcg_test1, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top5 <- DoHeatmap(avgexp, features = top5$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top5

Idents(ivbcg_test2) <- "PTID"
markers <- FindAllMarkers(ivbcg_test2)
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
avgexp <- AverageExpression(ivbcg_test2, assays = "RNA", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top5 <- DoHeatmap(avgexp, features = top5$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top5



#########################################
## Try the Seurat Integration Pipeline ##
#########################################

#Read in the Seurat objects
dg0d <- readRDS(file = "out/seurat_objects/dg0d_analysis.rds")
dg0h <- readRDS(file = "out/seurat_objects/dg0h_analysis.rds")
v113 <- readRDS(file = "out/seurat_objects/v113_analysis.rds") 
v015 <- readRDS(file = "out/seurat_objects/v015_analysis.rds")

DefaultAssay(dg0d) <- "RNA"
DefaultAssay(dg0h) <- "RNA"
DefaultAssay(v113) <- "RNA"
DefaultAssay(v015) <- "RNA"

Idents(dg0d) <- "Celltype"
Idents(dg0h) <- "Celltype"
Idents(v113) <- "Celltype"
Idents(v015) <- "Celltype"

# Add cell names
dg0d <- RenameCells(dg0d, add.cell.id = "dg0d")
dg0h <- RenameCells(dg0h, add.cell.id = "dg0h")
v113 <- RenameCells(v113, add.cell.id = "v113")
v015 <- RenameCells(v015, add.cell.id = "v015")

dg0d <- dg0d %>%
  subset(idents = "gd")
dg0h <- dg0h %>%
  subset(idents = "gd")
v113 <- v113 %>%
  subset(idents = "gd")
v015 <- v015 %>%
  subset(idents = "gd")


# Select features that are repeatedly variable across datasets for integration
list <- list(dg0d, dg0h, v113, v015)
features <- SelectIntegrationFeatures(object.list = list)


# Perform integration
PrepSCTIntegration(
  list,
  assay = NULL,
  anchor.features = 2000,
  sct.clip.range = NULL,
  verbose = TRUE
)

immune.anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)

ivbcg <- IntegrateData(anchorset = immune.anchors)


#Identify the frequency of each paired_tcr group in each cluster
ggplot(ivbcg@meta.data, aes(x=seurat_clusters, fill=Timepoint)) + geom_bar(position = "fill")

# Cluster the integrated data
DefaultAssay(ivbcg) <- "integrated"
ivbcg <- ScaleData(ivbcg, verbose = FALSE)
ivbcg <- RunPCA(ivbcg, npcs = 30, verbose = FALSE)
ivbcg <- RunUMAP(ivbcg, reduction = "pca", dims = 1:12)
ivbcg <- FindNeighbors(ivbcg, reduction = "pca", dims = 1:12)
ivbcg <- FindClusters(ivbcg, resolution = 0.4)

DimPlot(ivbcg, reduction = "umap")
DimPlot(ivbcg, reduction = "umap", group.by = "PTID")
DimPlot(ivbcg, reduction = "umap", group.by = "Tissue")
DimPlot(ivbcg, reduction = "umap", group.by = "Timepoint")

FeaturePlot(ivbcg, pt.size = 0.5, order = T, min.cutoff = 1.5, ncol = 3, 
            features = c("TRDV1", 
                         "TRDV2", 
                         "TRDV3")) &
  scale_color_viridis_c() &
  theme(text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


# Annotate the GD subsets
tcr_grouped <- read.csv(file = "out/tcr/tcr_grouped.csv")

group_bc <- tcr_grouped %>% select(c(barcode, tcr_group, paired_tcr))
group_bc2 <- group_bc[,-1]
rownames(group_bc2) <- group_bc$barcode
ivbcg <- AddMetaData(ivbcg, metadata = group_bc2)

# Visualize the GD subsets on the UMAP
DimPlot(ivbcg, reduction = "umap", group.by = "tcr_group")
Idents(ivbcg) <- "tcr_group"
ivbcg %>% subset(idents = c("Vd1")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = c("Vd3")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = c("Vd4")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = c("Vg9Vd2")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = c("Other_Vd2")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = c("Gamma_only")) %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")

# Make a new Vd1/3/4-only object
Idents(ivbcg) <- "tcr_group"
ivbcg_vd1 <- ivbcg %>%
  subset(idents = c("Vd1", "Vd3", "Vd4"))

DimPlot(ivbcg_vd1, reduction = "umap")

# Save the Seurat Objects
saveRDS(ivbcg, file = "out/seurat_objects/ivbcg_integrated_2024.rds")
saveRDS(ivbcg_vd1, file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")


#####################################################################################
##                                                                                 ##
##                              TEST ANALYSIS                                      ##
##                                                                                 ##
#####################################################################################

# Read in the Seurat Objects
ivbcg <- readRDS(file = "out/seurat_objects/ivbcg_integrated_2024.rds")
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds") 

# Quantify the number of cell barcodes in each TCR-defined group
Idents(ivbcg) <- "tcr_group"

ivbcg_vd1 <- ivbcg %>%
  subset(idents = "Vd1")

ivbcg_vd3 <- ivbcg %>%
  subset(idents = "Vd1")

ivbcg_vd1 <- ivbcg %>%
  subset(idents = "Vd1")


# Set color palettes
my_cols12 <- viridis(12, begin = 0.15)
my_cols6 <- viridis(6, end = 0.9)
my_cols4 <- viridis(4, end = 0.9)
my_cols3 <- viridis(3, end = 0.9)
my_cols2 <- viridis(2, begin = 0.3, end = 0.9)

# Visualize the UMAPs
Idents(ivbcg) <- "seurat_clusters"
DimPlot(ivbcg, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.1,
        label.size = 4,
        cols = my_cols12)
ggsave('out/images/fig2b.pdf', width = 4, height = 3, units = 'in')

Idents(ivbcg) <- "Tissue"
DimPlot(ivbcg, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.1,
        label.size = 5,
        cols = my_cols2)
ggsave('out/images/fig2c.pdf', width = 4.3, height = 3, units = 'in')

Idents(ivbcg) <- "PTID"
DimPlot(ivbcg, reduction = "umap", cols = my_cols4, shuffle = T)
ggsave('out/images/supp_ptid_locations.pdf')

Idents(ivbcg) <- "Tissue"
ivbcg_bal <- ivbcg %>% subset(idents = c('BAL'))
ivbcg_pbmc <- ivbcg %>% subset(idents = c('PBMC'))


Idents(ivbcg) <- "tcr_group"
ivbcg %>% subset(idents = c("Vd1", "Vg9Vd2", "Vd3")) %>%
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "tcr_group",
          pt.size = 0.01)
ggsave('out/images/fig2d.pdf', width = 9, height = 3, units = 'in')

ivbcg_pbmc %>% 
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "Timepoint",
          pt.size = 0.01)
ggsave('out/images/fig2e.pdf', width = 9, height = 3, units = 'in')

ivbcg_bal %>% 
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "Timepoint",
          pt.size = 0.01)
ggsave('out/images/fig2f.pdf', width = 6.5, height = 3, units = 'in')



# Visualize the DE genes between PTIDs

Idents(ivbcg) <- "Tissue"
ivbcg_test <- ivbcg %>%
  subset(idents = "PBMC")

Idents(ivbcg_test) <- "Timepoint"
ivbcg_test1 <- ivbcg_test %>%
  subset(idents = "Week_0")
ivbcg_test2 <- ivbcg_test %>%
  subset(idents = "Week_4")

Idents(ivbcg_test1) <- "PTID"
markers <- FindAllMarkers(ivbcg_test1)
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
avgexp <- AverageExpression(ivbcg_test1, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top5 <- DoHeatmap(avgexp, features = top5$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top5

Idents(ivbcg_test2) <- "PTID"
markers <- FindAllMarkers(ivbcg_test2)
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5
avgexp <- AverageExpression(ivbcg_test2, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top5 <- DoHeatmap(avgexp, features = top5$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top5


# Visualize the distribution of each PTID on the UMAP
Idents(ivbcg) <- "PTID"
ivbcg %>% subset(idents = "DG0D") %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = "DG0H") %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = "V015") %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")
ivbcg %>% subset(idents = "V113") %>% DimPlot(reduction = "umap", group.by = "seurat_clusters")

# Visualize some of the key markers on the UMAP
FeaturePlot(ivbcg, pt.size = 0.75, features = c("TRDV1", 
                                "TRDV2", 
                                "TRDV3", 
                                "adt_CD4", 
                                "adt_CD8a", 
                                "adt_CCR7", 
                                "adt_CD28", 
                                "SELL", 
                                "TCF7", 
                                "IFNG", 
                                "adt_CXCR3", 
                                "RORC", 
                                "RORA", 
                                "MAMU-DRB1", 
                                "GZMB", 
                                "NKG7")) &
  scale_color_viridis_c() &
  theme(text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


# Abbreviated Feature Plots
FeaturePlot(ivbcg, pt.size = 0.75, order = T, min.cutoff = 1.5, ncol = 4, 
                                   features = c("adt_CD28", 
                                                "MAMU-DRB1", 
                                                "NKG7",
                                                "GZMA",
                                                "GZMK",
                                                "PDCD1",
                                                "STAT1",
                                                "ITGA1")) &
  scale_color_viridis_c() &
  theme(text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


# Find all markers and view a heat map
Idents(ivbcg) <- "seurat_clusters"

markers <- FindAllMarkers(ivbcg)
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC) -> top10
markers %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=5, wt = avg_log2FC) -> top5

avgexp <- AverageExpression(ivbcg, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.

heatmap_top10 <- DoHeatmap(avgexp, features = top10$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 6)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top10
ggsave("out/images/integrated_hm.png", width = 5, height = 10, units = "in")

heatmap_top5 <- DoHeatmap(avgexp, features = top5$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 7)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top5
ggsave("out/images/integrated_hm2.png", width = 5, height = 5, units = "in")



#Generate a heat map of T cell genes of interest in each cluster
marker_list <- c("CD4", "CD8A", "CD8B", "TRDV1", "TRDV3", "TRDV2", "TRGV9", "TCF7", "CCR7", "SELL", "CD28", "IL7R", "CD69", "IL2RA", "MAMU-DRA", "MAMU-DRB1", "LAT", "PTPRC", "LCK", "ZAP70", "TIGIT", "KLRG1", "BTLA", "CTLA4", "HAVCR2", "LAG3", "PDCD1", "TOX", "TBX21", "IFNG", "CXCR3", "TNF", "CCR5", "STAT1", "IFNGR1", "GATA3", "STAT6", "IL4R", "RORA", "RORC", "CCR6", "STAT3", "BCL6", "CXCR5", "ICOS", "CD40LG", "TGFB1", "IL10", "FOXP3", "FGFBP2", "GNLY", "GZMB", "GZMH", "GZMA", "GZMK", "SRGN", "PRF1", "GZMM", "KLRD1", "NKG7", "ITGB1", "ZNF683","CST7", "KLRB1", "NCR3", "EOMES", "NCAM1", "IL10RA", "IL16", "IL18", "IL18R1", "IL21R", "IL6ST", "IL7", "CCL4", "CCL5", "CXCR3", "CXCL10", "CCL4L1", "CCR9", "CXCR4", "CXCR6")

avgexp <- AverageExpression(ivbcg, assays = "integrated", return.seurat = T)
avgexp <- ScaleData(avgexp, features = marker_list)
DoHeatmap(avgexp, features = marker_list, raster = FALSE, draw.lines = F) + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10))
ggsave("out/images/ivbcgvd1_int_goi_timepoint.png", width = 4, height = 8)


# Generate a heat map of the surface proteins in each cluster
Idents(ivbcg) <- "seurat_clusters"
DefaultAssay(ivbcg) <- "adt"
VariableFeatures(ivbcg) <- rownames(ivbcg[["adt"]])
ivbcg <- ScaleData(ivbcg)
avgexp <- AverageExpression(ivbcg, assays = "adt", return.seurat = T)
hm_order <- c("CD4", "CD8a", "CD8b", "CCR7", "CD28", "CD127-IL-7R", "CD25", "CD69", "HLA-DR", "PD-L1", "CD11b", "CD86", "CD11c", "CD161", "CD16",  "CCR6", "CD95", "CXCR3")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap




####################################
## Final Heat Map and Annotations ##
####################################

#Generate a heat map of genes of interest in each cluster
marker_list <- c("CCR7", "TCF7", "ID3",
                 "TNFRSF9", "NR4A1", "NR4A2", "NR4A3",
                 "ICOS", "ICAM1", "CD69", "MAMU-DRB1", "TIGIT", "ITGAX", "CD40LG", "TNFRSF4",
                 "TOX", "PDCD1", "CD96", "LAG3",
                 "XCL1", "MX1", "STAT1", "IFNG", "TNF",
                 "RORA", "RORC", "IL21", "CCR6",
                 "ITGB1", "NCAM1", "TNFSF8","TNFRSF10B", "GZMK", "GZMM", "NKG7", "CST7", "GZMB", "GZMA", "PRF1", "TNFRSF18",
                 "ITGA1", "ITGB7", "ITGAE"
                 )
avgexp <- AverageExpression(ivbcg, assays = "integrated", return.seurat = T)
avgexp <- ScaleData(avgexp, features = marker_list)
DoHeatmap(avgexp, features = marker_list, 
          raster = FALSE, 
          draw.lines = F,
          group.colors = my_cols12,
          group.bar.height = 0.03,
          label = F) + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) +
  guides(color = F)
ggsave("out/images/hm_240522.pdf", units = c("in"), width = 4, height = 3.5)


# Generate a heat map in dot plot form
Idents(ivbcg) <- "seurat_clusters"

DotPlot(
  ivbcg,
  features = marker_list,
  assay = 'integrated',
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0.0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA) + 
  RotatedAxis() +
  scale_colour_gradient2(low = 'steelblue', mid = 'white', high = '#E54500')

ggsave('out/images/fig1g.pdf', units = c("in"), width = 11.5, height = 4)



# Generate a heat map of surface proteins in each cluster
Idents(ivbcg) <- "seurat_clusters"
DefaultAssay(ivbcg) <- "adt"
VariableFeatures(ivbcg) <- rownames(ivbcg[["adt"]])
ivbcg <- ScaleData(ivbcg)
avgexp <- AverageExpression(ivbcg, assays = "adt", return.seurat = T)
hm_order <- c("CD4", "CD8a", "CD28", "CD69", "HLA-DR", "CD11b", "CD161", "CXCR3", "CD95")
DoHeatmap(avgexp, features = hm_order, 
          raster = FALSE, 
          draw.lines = F,
          group.colors = my_cols12,
          group.bar.height = 0.05,
          label = F) + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) +
  guides(color = F)
ggsave("out/images/hm_adt_240522.pdf", units = c("in"), width = 4, height = 1.5)

# Generate a dot plot of surface proteins in each cluster
DotPlot(
  ivbcg,
  features = hm_order,
  assay = 'adt',
  col.min = -2.5,
  col.max = 2.5,
  dot.min = 0,
  dot.scale = 6,
  idents = NULL,
  group.by = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.min = NA,
  scale.max = NA) + 
  RotatedAxis() +
  scale_colour_gradient2(low = 'steelblue', mid = 'white', high = '#E54500') +
  scale_size(breaks = c(25, 50, 75, 100))
ggsave("out/images/hm_adt_240522.pdf", units = c("in"), width = 4.8, height = 3.8)



#############################################
## Differential Abundance Analysis -- PBMC ##
#############################################

#Identify the frequency of each paired_tcr group in each cluster
ivbcg@meta.data %>% filter(tcr_group != "Gamma_only") %>% ggplot(aes(x=seurat_clusters, fill=tcr_group)) + geom_bar(position = "fill") + scale_fill_viridis_d()

# Export the data from the Seurat object and identify the proportion of cells in each cluster over time
ivbcg_md <- ivbcg@meta.data

summary <- ivbcg_md %>% group_by(Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))

summary %>% filter(Tissue == 'PBMC') %>%
  ggplot(aes(x = Timepoint, y = count, fill = seurat_clusters)) +
  geom_bar(position = 'fill', stat = 'identity')

summary$seurat_clusters <- as.character(summary$seurat_clusters)

summary <- summary %>% mutate(clust2 = case_when(
  seurat_clusters %in% c('6', '7', '8', '9', '10', '11') ~ '6+',
  .default = seurat_clusters)
  )
                               
summary %>% filter(Tissue == 'PBMC') %>%
  ggplot(aes(x = Timepoint, y = count, fill = clust2)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual(values = c("#463480FF", "#3E4C8AFF", "black",  "#2C738EFF", "#24868EFF", "#1F998AFF", 'lightgrey')) + theme_bw()
ggsave('out/images/fig2h.pdf', units = c("in"), width = 3.5, height = 3)

  
sum_p <- summary %>% filter(Tissue == 'PBMC')
table(sum_p$Timepoint, sum_p$clust2)



# Plots of cluster Freq over time

summary <- ivbcg_md %>% group_by(PTID, Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))

a <- summary %>% filter(seurat_clusters == 0 & Tissue == 'PBMC') %>% arrange(Timepoint) %>%
    ggpaired(x = "Timepoint", y = "freq", id = "PTID",
           line.color = "black",
           line.size = 1, 
           point.size = 1,
           short.panel.labs = FALSE) +
  ggtitle("Cluster 0, 'Resting'")+
  ylab("Prop. of GD T cells") +
  xlab("Timepoint")

a$layers <- a$layers[-1]

 
b <- summary %>% filter(seurat_clusters == 3 & Tissue == 'PBMC') %>% arrange(Timepoint) %>%
  ggpaired(x = "Timepoint", y = "freq", id = "PTID",
           line.color = "black",
           line.size = 1, 
           point.size = 1,
           short.panel.labs = FALSE) +
  ggtitle("Cluster 3, 'GZMK+/GZMM+'")+
  ylab("Prop. of GD T cells") +
  xlab("Timepoint")
b$layers <- b$layers[-1]


c <- summary %>% filter(seurat_clusters == 4 & Tissue == 'PBMC') %>% arrange(Timepoint) %>%
  ggpaired(x = "Timepoint", y = "freq", id = "PTID",
           line.color = "black",
           line.size = 1, 
           point.size = 1,
           short.panel.labs = FALSE) +
  ggtitle("Cluster 4, 'GZMA+/PRF1+'")+
  ylab("Prop. of GD T cells") +
  xlab("Timepoint")
c$layers <- c$layers[-1]

ggarrange(a, b, c, ncol = 3, legend = F)
ggsave('out/images/fig2i.pdf', width = 9, height = 3, units = 'in')



############################################
## Differential Abundance Analysis -- BAL ##
############################################


# Export the data from the Seurat object and identify the proportion of cells in each cluster over time
ivbcg_md <- ivbcg@meta.data

summary <- ivbcg_md %>% group_by(Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))

summary %>% filter(Tissue == 'BAL') %>%
  ggplot(aes(x = Timepoint, y = count, fill = seurat_clusters)) +
  geom_bar(position = 'fill', stat = 'identity')

summary$seurat_clusters <- as.character(summary$seurat_clusters)

summary <- summary %>% mutate(clust2 = case_when(
  seurat_clusters %in% c('3', '4', '5', '6', '7', '8', '9', '10', '11') ~ '3+',
  .default = seurat_clusters)
)

summary %>% filter(Tissue == 'BAL') %>%
  ggplot(aes(x = Timepoint, y = count, fill = clust2)) +
  geom_bar(position = 'fill', stat = 'identity') +
  scale_fill_manual(values = c("#463480FF", "black", "#24868EFF", 'lightgrey')) + theme_bw()

sum_p <- summary %>% filter(Tissue == 'PBMC')
table(sum_p$Timepoint, sum_p$clust2)
ggsave('out/images/supp_bal_clust_abundance.pdf', width = 3, height = 3, units = 'in')


# Box plots of cluster Freq over time

summary <- ivbcg_md %>% group_by(PTID, Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))

a <- summary %>% filter(seurat_clusters == 0 & Tissue == 'BAL') %>% arrange(Timepoint) %>%
  ggpaired(x = "Timepoint", y = "freq", id = "PTID",
           line.color = "black",
           line.size = 1, 
           point.size = 1,
           short.panel.labs = FALSE) +
  ggtitle("Cluster 0")+
  ylab("Prop. of GD T cells") +
  xlab("Timepoint") +
  stat_compare_means(method = 't.test', 
                     paired = T,
                     ref.group = 'Week_4',
                     label = 'p.signif')

a$layers <- a$layers[-1]


b <- summary %>% filter(seurat_clusters == 2 & Tissue == 'BAL') %>% arrange(Timepoint) %>%
  ggpaired(x = "Timepoint", y = "freq", id = "PTID",
           line.color = "black",
           line.size = 1, 
           point.size = 1,
           short.panel.labs = FALSE) +
  ggtitle("Cluster 2")+
  ylab("Prop. of GD T cells") +
  xlab("Timepoint") +
  stat_compare_means(method = 't.test', 
                     paired = T,
                     ref.group = 'Week_4',
                     label = 'p.signif')
b$layers <- b$layers[-1]

ggarrange(a, b, ncol = 2, legend = F)
ggsave('out/images/supp_bal_clust_lineplots.pdf', width = 5, height = 3, units = 'in')



################################
## DE GENES ACROSS CONDITIONS ##
################################

#Identify the DE genes across Week 0, Week 4, and Week 8 PBMC
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("PBMC"))
Idents(ivbcg_vd1) <- "Timepoint"

markers <- FindAllMarkers(ivbcg_vd1)
markers %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30

top30 <- top30 %>% filter(!grepl('ENSMMUG', gene))

avgexp <- AverageExpression(ivbcg_vd1, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30$gene, raster = FALSE, draw.lines = F) + NoLegend() + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30



## MAKE A VOLCANO PLOT OF THE GENES UP AND DOWN AFTER VACCINATION ##
## PRE VS POST ##
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

# Omit the BAL data
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("PBMC"))

# Annotate the cell barcodes as "pre" or "post"
Idents(ivbcg_vd1) <- "Timepoint"
pre <- ivbcg_vd1 %>% subset(idents = c("Week_0"))
post <- ivbcg_vd1 %>% subset(idents = c("Week_4", "Week_8"))

pre <- pre@meta.data
post <- post@meta.data

pre$Vax <- "Pre"
post$Vax <- "Post" 

anno <- rbind(pre, post) %>% subset(select = c("Vax"))

ivbcg_vd1 <- AddMetaData(ivbcg_vd1, metadata = anno)


# Find DE genes between Week 0 and Week 4 PBMC
Idents(ivbcg_vd1) <- "Timepoint"
markers <- FindAllMarkers(ivbcg_vd1)
markers <- markers %>% filter(!grepl('ENSMMUG', gene)) %>% filter(cluster == "Week_4")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC> 0.6 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -0.6 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  scale_color_manual(values = c("blue", "gray", "red")) +
   geom_text_repel(label = markers$delabel)

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "azure4") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)

ggsave("out/images/volcano_prepost.png", width = 5, height = 3, bg = "white")


# Using the same list of markers, make a heat map

markers <- markers %>% filter(p_val_adj < 0.05 & (avg_log2FC < -0.6 | avg_log2FC > 0.6))
write.csv(markers, file = 'out/Week4_vs_Week0_PBMC_genes.csv')






## MAKE A VOLCANO PLOT OF THE GENES UP AND DOWN AFTER VACCINATION ##
## PRE VS WEEK 4 ONLY ##
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")

# Omit the BAL data and Week 8 data
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("PBMC"))

Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("Week_0", "Week_4"))

# Find DE genes between Pre and Post
markers <- FindAllMarkers(ivbcg_vd1)
markers <- markers %>% filter(!grepl('ENSMMUG', gene)) %>% filter(cluster == "Week_4")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC> 0.6 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -0.6 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)



## MAKE A VOLCANO PLOT OF THE GENES UP AND DOWN IN BAL VS PBMC ##

## WEEK 4 ##
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("Week_4"))

Idents(ivbcg_vd1) <- "Tissue"

markers <- FindAllMarkers(ivbcg_vd1)
markers <- markers %>% filter(!grepl('ENSMMUG', gene)) %>% filter(cluster == "BAL")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC > 0.6 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -0.6 & markers$p_val_adj < 0.05] <- "DOWN"
markers$p_val_adj[markers$p_val_adj< 0.05]

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)


## WEEK 8 ##
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")

Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("Week_8"))

Idents(ivbcg_vd1) <- "Tissue"

markers <- FindAllMarkers(ivbcg_vd1)
markers <- markers %>% filter(!grepl('ENSMMUG', gene)) %>% filter(cluster == "BAL")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC > 0.6 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -0.6 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)


## MAKE A VOLCANO PLOT OF THE GENES UP AND DOWN IN WEEK 8 BAL VS WEEK 4 BAL ##

ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")

Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("BAL"))

Idents(ivbcg_vd1) <- "Timepoint"

markers <- FindAllMarkers(ivbcg_vd1)
markers <- markers %>% filter(!grepl('ENSMMUG', gene)) %>% filter(cluster == "Week_8")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC > 0.6 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -0.6 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "black") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)




#######################################
## HEATMAPS INSTEAD OF VOLCANO PLOTS ##
#######################################

ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")


## PBMC, compare timepoints
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1_p <- ivbcg_vd1 %>% subset(idents = c('PBMC'))
Idents(ivbcg_vd1_p) <- "Timepoint"
markers_pbmc <- FindAllMarkers(ivbcg_vd1_p)

markers_pbmc %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_pbmc
top30_pbmc <- top30_pbmc %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

avgexp <- AverageExpression(ivbcg_vd1_p, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_pbmc$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_pbmc, 'out/markers_vd1_pbmc.csv')
write.csv(top30_pbmc, 'out/top30_vd1_pbmc.csv')


## BAL, compare timepoints
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1_b <- ivbcg_vd1 %>% subset(idents = c('BAL'))
Idents(ivbcg_vd1_b) <- "Timepoint"
markers_bal <- FindAllMarkers(ivbcg_vd1_b)

markers_bal %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_bal
top30_bal <- top30_bal %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

avgexp <- AverageExpression(ivbcg_vd1_b, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_bal$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_bal, 'out/markers_vd1_bal.csv')
write.csv(top30_bal, 'out/top30_vd1_bal.csv')


## PBMC v BAL, Week 4
Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1_wk4 <- ivbcg_vd1 %>% subset(idents = c('Week_4'))
Idents(ivbcg_vd1_wk4) <- "Tissue"
markers_wk4 <- FindAllMarkers(ivbcg_vd1_wk4)

markers_wk4 %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_wk4
top30_wk4 <- top30_wk4 %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

avgexp <- AverageExpression(ivbcg_vd1_wk4, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_wk4$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_wk4, 'out/markers_vd1_wk4.csv')
write.csv(top30_wk4, 'out/top30_vd1_wk4.csv')


## PBMC v BAL, Week 8
Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1_wk8 <- ivbcg_vd1 %>% subset(idents = c('Week_8'))
Idents(ivbcg_vd1_wk8) <- "Tissue"
markers_wk8 <- FindAllMarkers(ivbcg_vd1_wk8)

markers_wk8 %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_wk8
top30_wk8 <- top30_wk8 %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

avgexp <- AverageExpression(ivbcg_vd1_wk8, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_wk8$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_wk8, 'out/markers_vd1_wk8.csv')
write.csv(top30_wk8, 'out/top30_vd1_wk8.csv')


## ALL TISSUES AND TIMEPOINTS

ivbcg_vd1_md <- ivbcg_vd1@meta.data
ivbcg_vd1_md$SampleType <- paste(ivbcg_vd1_md$Tissue, ivbcg_vd1_md$Timepoint, sep = "_")
ST <- ivbcg_vd1_md %>% select(c('SampleType'))
ivbcg_vd1 <- AddMetaData(ivbcg_vd1, metadata = ST)

Idents(ivbcg_vd1) <- "SampleType"
markers_all <- FindAllMarkers(ivbcg_vd1)

markers_all %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_all
top30_all <- top30_all %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

avgexp <- AverageExpression(ivbcg_vd1, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_all$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_all, 'out/markers_vd1_all.csv')
write.csv(top30_all, 'out/top30_vd1_all.csv')

filtered <- read.csv('out/filtered_samptype_markers2.csv')
filtered <- filtered$gene

test <- markers_all[markers_all$gene %in% filtered, ]
test <- test %>% arrange(cluster, desc(avg_log2FC))

heatmap_filtered <- DoHeatmap(avgexp, features = test$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_filtered


# Filter to top FC genes
test %>% # Generate a list of the top 30 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_test
heatmap_filtered <- DoHeatmap(avgexp, features = top30_test$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 8)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_filtered

test %>% # Generate a list of the top 20 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=20, wt = avg_log2FC) -> top20_test
heatmap_filtered <- DoHeatmap(avgexp, features = top20_test$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_filtered

test %>% # Generate a list of the top 10 genes in each cluster.
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC) -> top10_test
heatmap_filtered <- DoHeatmap(avgexp, features = top10_test$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_filtered


## FINALIZED HEAT MAP

filtered <- read.csv('out/filtered_samptype_markers2.csv')
filtered <- filtered$gene

heatmap_filtered <- DoHeatmap(avgexp, features = filtered, raster = FALSE, draw.lines = F) + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 11))

heatmap_filtered
ggsave('out/images/fig2j.pdf', width = 4, height = 6.5, units = 'in')






#######################################
## DIFFERENTIALLY EXPRESSED PROTEINS ##
#######################################

## DE PROTEINS BETWEEN TIMEPOINTS IN PBMC

# Generate a heat map of surface proteins in each cluster
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")

Idents(ivbcg_vd1) <- "Tissue"
temp <- ivbcg_vd1 %>% subset(idents = c("PBMC"))

Idents(temp) <- "Timepoint"

DefaultAssay(temp) <- "adt"
VariableFeatures(temp) <- rownames(ivbcg[["adt"]])
temp <- ScaleData(temp)
avgexp <- AverageExpression(temp, assays = "adt", return.seurat = T)
hm_order <- c("CD4", "CD8a", "CD8b", "SELL", "CCR7", "CD28", "CD69", "CD25", "HLA-DR", "CD11b", "CCR6", "CD16", "CD11c", "CD161")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap


## DE PROTEINS BETWEEN TIMEPOINTS IN BAL

# Generate a heat map of surface proteins in each cluster
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")

Idents(ivbcg_vd1) <- "Tissue"
temp <- ivbcg_vd1 %>% subset(idents = c("BAL"))

Idents(temp) <- "Timepoint"

DefaultAssay(temp) <- "adt"
VariableFeatures(temp) <- rownames(ivbcg[["adt"]])
temp <- ScaleData(temp)
avgexp <- AverageExpression(temp, assays = "adt", return.seurat = T)
hm_order <- c("CD4", "CD8a", "CD8b", "SELL", "CCR7", "CD28", "CD69", "CD25", "HLA-DR", "CD11b", "CCR6", "CD16", "CD11c", "CD161")
heatmap <- DoHeatmap(avgexp, features = hm_order, raster = FALSE, draw.lines = F) + scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap




##########
## GSEA ##
##########

######################
## Hypothesis-based ##
######################

library(fgsea)
library(presto)
library(msigdbr)
library(scales)
library(dplyr)
library(Seurat)
library(tidyverse)

ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

# Exclude the BAL cells for now
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = "PBMC")

# Generate a reference table containing the human gene sets to which human genes are assigned. Genes appear more than once depending on the number of gene sets they are associated with. We are only looking in the MSigDB C7 category (immunologic signatures).
msig_<- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")

sets <- msig_%>% filter(
  str_detect(gs_name, ("_T_CELL"))|str_detect(gs_name, ("CYTOTOXIC"))
)

sets <- sets %>% split(x = .$gene_symbol, f = .$gs_name)


# GSEA -- ALL COMBINED
Idents(ivbcg_vd1) <- "Timepoint"

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes <- wilcoxauc(ivbcg_vd1, 'Timepoint')


## T cell  gene sets ##
# Week 0
wk0_genes <- fgsea.genes %>%
  dplyr::filter(group == "Week_0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk0.ranks <- deframe(wk0_genes) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk0 <- fgsea(sets, stats = wk0.ranks, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk0 <- gsea.wk0 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 4
wk4_genes <- fgsea.genes %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk4.ranks <- deframe(wk4_genes) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk4 <- fgsea(sets, stats = wk4.ranks, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk4 <- gsea.wk4 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


# Week 8
wk8_genes <- fgsea.genes %>%
  dplyr::filter(group == "Week_8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk8.ranks <- deframe(wk8_genes) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk8 <- fgsea(sets, stats = wk8.ranks, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk8 <- gsea.wk8 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))



# Try making a dot plot
gsea.wk0$timepoint <- "Week_0"
gsea.wk4$timepoint <- "Week_4"
gsea.wk8$timepoint <- "Week_8"
gsea.all <- bind_rows(gsea.wk0, gsea.wk4, gsea.wk8)




gsea.all %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.01)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500")



gsea.all %>% filter(pathway == "GOBP_T_CELL_ACTIVATION"|pathway == "GOBP_T_CELL_PROLIFERATION"|pathway == "GOBP_NK_T_CELL_DIFFERENTIATION"|pathway == "GOBP_T_CELL_MIGRATION"|pathway == "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY") %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500")




# Break into a new object for each PTID #
## Just compare Week 0 vs. Week 4, since those timepoints have the most differences ##

Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = c("Week_0", "Week_4"))

Idents(ivbcg_vd1) <- "PTID"
dg0d <- ivbcg_vd1 %>% subset(idents = "DG0D")
dg0h <- ivbcg_vd1 %>% subset(idents = "DG0H")
v113 <- ivbcg_vd1 %>% subset(idents = "V113")
v015 <- ivbcg_vd1 %>% subset(idents = "V015")

## DG0D ##
Idents(dg0d) <- "Timepoint"
DimPlot(dg0d, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0d <- wilcoxauc(dg0d, 'Timepoint')

## T cell  gene sets ##
# Week 0
wk0_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "Week_0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk0.ranks.dg0d <- deframe(wk0_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk0.dg0d <- fgsea(sets, stats = wk0.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk0.dg0d <- gsea.wk0.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 4
wk4_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk4.ranks.dg0d <- deframe(wk4_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk4.dg0d <- fgsea(sets, stats = wk4.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk4.dg0d <- gsea.wk4.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Dotplot
gsea.wk0.dg0d$timepoint <- "Week_0"
gsea.wk4.dg0d$timepoint <- "Week_4"
gsea.all.dg0d <- bind_rows(gsea.wk0.dg0d, gsea.wk4.dg0d)

gsea.all.dg0d %>% filter(padj < 0.05) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.01)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500")



## DG0H ##

Idents(dg0h) <- "Timepoint"
DimPlot(dg0h, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0h <- wilcoxauc(dg0h, 'Timepoint')


## T cell  gene sets ##
# Week 0
wk0_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "Week_0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk0.ranks.dg0h <- deframe(wk0_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk0.dg0h <- fgsea(sets, stats = wk0.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk0.dg0h <- gsea.wk0.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 4
wk4_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk4.ranks.dg0h <- deframe(wk4_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk4.dg0h <- fgsea(sets, stats = wk4.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk4.dg0h <- gsea.wk4.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))



## v113 ##

Idents(v113) <- "Timepoint"
DimPlot(v113, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v113 <- wilcoxauc(v113, 'Timepoint')


## T cell  gene sets ##
# Week 0
wk0_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "Week_0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk0.ranks.v113 <- deframe(wk0_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk0.v113 <- fgsea(sets, stats = wk0.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk0.v113 <- gsea.wk0.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 4
wk4_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk4.ranks.v113 <- deframe(wk4_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk4.v113 <- fgsea(sets, stats = wk4.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk4.v113 <- gsea.wk4.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))



## v015 ##

Idents(v015) <- "Timepoint"
DimPlot(v015, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v015 <- wilcoxauc(v015, 'Timepoint')


## T cell  gene sets ##
# Week 0
wk0_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "Week_0") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk0.ranks.v015 <- deframe(wk0_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk0.v015 <- fgsea(sets, stats = wk0.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk0.v015 <- gsea.wk0.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 4
wk4_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)

wk4.ranks.v015 <- deframe(wk4_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.

gsea.wk4.v015 <- fgsea(sets, stats = wk4.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster

gsea.wk4.v015 <- gsea.wk4.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DOT PLOT SEPARATED BY PTID ##
gsea.wk0.dg0d$timepoint <- "Week_0"
gsea.wk0.dg0d$ptid <- "DG0D"
gsea.wk4.dg0d$timepoint <- "Week_4"
gsea.wk4.dg0d$ptid <- "DG0D"


gsea.wk0.dg0h$timepoint <- "Week_0"
gsea.wk0.dg0h$ptid <- "DG0H"
gsea.wk4.dg0h$timepoint <- "Week_4"
gsea.wk4.dg0h$ptid <- "DG0H"


gsea.wk0.v113$timepoint <- "Week_0"
gsea.wk0.v113$ptid <- "V113"
gsea.wk4.v113$timepoint <- "Week_4"
gsea.wk4.v113$ptid <- "V113"


gsea.wk0.v015$timepoint <- "Week_0"
gsea.wk0.v015$ptid <- "V015"
gsea.wk4.v015$timepoint <- "Week_4"
gsea.wk4.v015$ptid <- "V015"


gsea.all <- bind_rows(gsea.wk0.dg0d, 
                      gsea.wk4.dg0d, 
                      gsea.wk0.dg0h, 
                      gsea.wk4.dg0h, 
                      gsea.wk0.v113, 
                      gsea.wk4.v113, 
                      gsea.wk0.v015, 
                      gsea.wk4.v015)

gsea.all$id <- paste(gsea.all$ptid, gsea.all$timepoint, sep = "_")


gsea.all %>% filter(padj < 0.05) %>%
   ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05, facet_by = "ptid")) + 
  geom_point(size = 5) + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)


ggsave("out/images/gsea_all_faceted.png", width = 16, height = 10, units = c("in"))

gsea.all$pathway <- gsub("GOBP_", "", as.character(gsea.all$pathway))
gsea.all$pathway <- gsub("_", " ", as.character(gsea.all$pathway))
gsea.all$timepoint <- gsub("Week_0", "Pre", as.character(gsea.all$timepoint))
gsea.all$timepoint <- gsub("Week_4", "Wk.4", as.character(gsea.all$timepoint))

gsea.all %>%
  filter(pathway == "ALPHA BETA T CELL ACTIVATION"|
           pathway == "T CELL PROLIFERATION"|
           pathway == "NK T CELL DIFFERENTIATION"|
           pathway == "ALPHA BETA T CELL DIFFERENTIATION"|
           pathway == "LEUKOCYTE MEDIATED CYTOTOXICITY"|
           pathway == "T CELL MEDIATED IMMUNITY"|
           pathway == "POSITIVE THYMIC T CELL SELECTION"
       ) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05, facet_by = "ptid")) + 
  scale_y_discrete(labels = label_wrap(30)) +
  geom_point() + 
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid) +
  theme_bw()




# Visualize the NES plots for the pathways that are consistently up after IV-BCG

a <- plotEnrichment(sets$GOBP_T_CELL_PROLIFERATION, wk4.ranks.dg0d)
b <- plotEnrichment(sets$GOBP_T_CELL_PROLIFERATION, wk4.ranks.dg0h)
c <- plotEnrichment(sets$GOBP_T_CELL_PROLIFERATION, wk4.ranks.v113)
d <- plotEnrichment(sets$GOBP_T_CELL_PROLIFERATION, wk4.ranks.v015)
ggarrange(a, b, c, d, nrow = 2, ncol = 2, labels = c("DG0D", "DG0H", "V113", "V015"))

a <- plotEnrichment(sets$GOBP_POSITIVE_THYMIC_T_CELL_SELECTION, wk4.ranks.dg0d)
b <- plotEnrichment(sets$GOBP_POSITIVE_THYMIC_T_CELL_SELECTION, wk4.ranks.dg0h)
c <- plotEnrichment(sets$GOBP_POSITIVE_THYMIC_T_CELL_SELECTION, wk4.ranks.v113)
d <- plotEnrichment(sets$GOBP_POSITIVE_THYMIC_T_CELL_SELECTION, wk4.ranks.v015)
ggarrange(a, b, c, d, nrow = 2, ncol = 2, labels = c("DG0D", "DG0H", "V113", "V015"))

a <- plotEnrichment(sets$GOBP_NK_T_CELL_DIFFERENTIATION, wk4.ranks.dg0d)
b <- plotEnrichment(sets$GOBP_NK_T_CELL_DIFFERENTIATION, wk4.ranks.dg0h)
c <- plotEnrichment(sets$GOBP_NK_T_CELL_DIFFERENTIATION, wk4.ranks.v113)
d <- plotEnrichment(sets$GOBP_NK_T_CELL_DIFFERENTIATION, wk4.ranks.v015)
ggarrange(a, b, c, d, nrow = 2, ncol = 2, labels = c("DG0D", "DG0H", "V113", "V015"))

a <- plotEnrichment(sets$GOBP_ALPHA_BETA_T_CELL_ACTIVATION, wk4.ranks.dg0d)
b <- plotEnrichment(sets$GOBP_ALPHA_BETA_T_CELL_ACTIVATION, wk4.ranks.dg0h)
c <- plotEnrichment(sets$GOBP_ALPHA_BETA_T_CELL_ACTIVATION, wk4.ranks.v113)
d <- plotEnrichment(sets$GOBP_ALPHA_BETA_T_CELL_ACTIVATION, wk4.ranks.v015)
ggarrange(a, b, c, d, nrow = 2, ncol = 2, labels = c("DG0D", "DG0H", "V113", "V015"))


# Review the leading edge genes

gsea.wk4.dg0d <- gsea.wk4.dg0d %>% filter(pathway == "GOBP_T_CELL_PROLIFERATION"|
                                            pathway == "GOBP_POSITIVE_THYMIC_T_CELL_SELECTION"|
                                            pathway == "GOBP_NK_T_CELL_DIFFERENTIATION" |
                                            pathway == "GOBP_ALPHA_BETA_T_CELL_ACTIVATION")
gsea.wk4.dg0d$leadingEdge

gsea.wk4.dg0h <- gsea.wk4.dg0h %>% filter(pathway == "GOBP_T_CELL_PROLIFERATION"|
                                            pathway == "GOBP_POSITIVE_THYMIC_T_CELL_SELECTION"|
                                            pathway == "GOBP_NK_T_CELL_DIFFERENTIATION" |
                                            pathway == "GOBP_ALPHA_BETA_T_CELL_ACTIVATION")
gsea.wk4.dg0h$leadingEdge

gsea.wk4.v113 <- gsea.wk4.v113 %>% filter(pathway == "GOBP_T_CELL_PROLIFERATION"|
                                            pathway == "GOBP_POSITIVE_THYMIC_T_CELL_SELECTION"|
                                            pathway == "GOBP_NK_T_CELL_DIFFERENTIATION" |
                                            pathway == "GOBP_ALPHA_BETA_T_CELL_ACTIVATION")
gsea.wk4.v113$leadingEdge

gsea.wk4.v015 <- gsea.wk4.v015 %>% filter(pathway == "GOBP_T_CELL_PROLIFERATION"|
                                            pathway == "GOBP_POSITIVE_THYMIC_T_CELL_SELECTION"|
                                            pathway == "GOBP_NK_T_CELL_DIFFERENTIATION" |
                                            pathway == "GOBP_ALPHA_BETA_T_CELL_ACTIVATION")
gsea.wk4.v015$leadingEdge


# Try making a dot plot
gsea.wk0$timepoint <- "Week_0"
gsea.wk4$timepoint <- "Week_4"
gsea.wk8$timepoint <- "Week_8"
gsea.all <- bind_rows(gsea.wk0, gsea.wk4, gsea.wk8)




gsea.all %>% filter(padj < 0.05) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.01)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500")



gsea.all %>% filter(pathway == "T CELL ACTIVATION"|pathway == "T CELL PROLIFERATION"|pathway == "NK T CELL DIFFERENTIATION"|pathway == "T CELL MIGRATION"|pathway == "T CELL MEDIATED CYTOTOXICITY"|pathway == "LEUKOCYTE MEDIATED CYTOTOXICITY") %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500")



##########
## GSEA ##
##########

ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

# BAL WEEK 8 VS PBMC WEEK 8
Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = "Week_8")

# Generate a reference table containing the human gene sets to which human genes are assigned. Genes appear more than once depending on the number of gene sets they are associated with. We are only looking in the MSigDB C7 category (immunologic signatures).
msig_<- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
sets <- msig_%>% filter(
  str_detect(gs_name, ("_T_CELL"))|str_detect(gs_name, ("CYTOTOXIC")))
sets <- sets %>% split(x = .$gene_symbol, f = .$gs_name)

# Break into a new object for each PTID #

Idents(ivbcg_vd1) <- "PTID"
dg0d <- ivbcg_vd1 %>% subset(idents = "DG0D")
dg0h <- ivbcg_vd1 %>% subset(idents = "DG0H")
v113 <- ivbcg_vd1 %>% subset(idents = "V113")
v015 <- ivbcg_vd1 %>% subset(idents = "V015")


## DG0D ##
Idents(dg0d) <- "Tissue"
DimPlot(dg0d, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0d <- wilcoxauc(dg0d, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.dg0d <- deframe(bal_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.dg0d <- fgsea(sets, stats = bal.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.dg0d <- gsea.bal.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.dg0d <- deframe(pbmc_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.dg0d <- fgsea(sets, stats = pbmc.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.dg0d <- gsea.pbmc.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DG0H ##
Idents(dg0h) <- "Tissue"
DimPlot(dg0h, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0h <- wilcoxauc(dg0h, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.dg0h <- deframe(bal_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.dg0h <- fgsea(sets, stats = bal.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.dg0h <- gsea.bal.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.dg0h <- deframe(pbmc_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.dg0h <- fgsea(sets, stats = pbmc.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.dg0h <- gsea.pbmc.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V113 ##
Idents(v113) <- "Tissue"
DimPlot(v113, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v113 <- wilcoxauc(v113, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.v113 <- deframe(bal_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.v113 <- fgsea(sets, stats = bal.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.v113 <- gsea.bal.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.v113 <- deframe(pbmc_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.v113 <- fgsea(sets, stats = pbmc.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.v113 <- gsea.pbmc.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V015 ##
Idents(v015) <- "Tissue"
DimPlot(v015, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v015 <- wilcoxauc(v015, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.v015 <- deframe(bal_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.v015 <- fgsea(sets, stats = bal.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.v015 <- gsea.bal.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.v015 <- deframe(pbmc_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.v015 <- fgsea(sets, stats = pbmc.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.v015 <- gsea.pbmc.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DOT PLOT SEPARATED BY PTID ##
gsea.bal.dg0d$tissue <- "BAL"
gsea.bal.dg0d$ptid <- "DG0D"
gsea.pbmc.dg0d$tissue <- "PBMC"
gsea.pbmc.dg0d$ptid <- "DG0D"


gsea.bal.dg0h$tissue <- "BAL"
gsea.bal.dg0h$ptid <- "DG0H"
gsea.pbmc.dg0h$tissue <- "PBMC"
gsea.pbmc.dg0h$ptid <- "DG0H"


gsea.bal.v113$tissue <- "BAL"
gsea.bal.v113$ptid <- "V113"
gsea.pbmc.v113$tissue <- "PBMC"
gsea.pbmc.v113$ptid <- "V113"


gsea.bal.v015$tissue <- "BAL"
gsea.bal.v015$ptid <- "V015"
gsea.pbmc.v015$tissue <- "PBMC"
gsea.pbmc.v015$ptid <- "V015"


gsea.all <- bind_rows(gsea.bal.dg0d, 
                      gsea.pbmc.dg0d, 
                      gsea.bal.dg0h, 
                      gsea.pbmc.dg0h, 
                      gsea.bal.v113, 
                      gsea.pbmc.v113, 
                      gsea.bal.v015, 
                      gsea.pbmc.v015)

gsea.all$id <- paste(gsea.all$ptid, gsea.all$tissue, sep = "_")


gsea.all %>% filter(padj < 0.05) %>%
  ggplot(aes(x = tissue, y = pathway, color = NES, facet_by = "ptid")) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)


ggsave("out/images/gsea_all_faceted.png", width = 16, height = 10, units = c("in"))

gsea.all$pathway <- gsub("GOBP_", "", as.character(gsea.all$pathway))
gsea.all$pathway <- gsub("_", " ", as.character(gsea.all$pathway))

gsea.all %>%
  filter(pathway == "ALPHA BETA T CELL ACTIVATION"|
           pathway == "T CELL PROLIFERATION"|
           pathway == "NK T CELL DIFFERENTIATION"|
           pathway == "ALPHA BETA T CELL DIFFERENTIATION"|
           pathway == "LEUKOCYTE MEDIATED CYTOTOXICITY"|
           pathway == "T CELL MEDIATED IMMUNITY"|
           pathway == "POSITIVE THYMIC T CELL SELECTION"
  ) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05, facet_by = "ptid")) + 
  scale_y_discrete(labels = label_wrap(30)) +
  geom_point() + 
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)








# BAL WEEK 4 VS PBMC WEEK 4
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")
Idents(ivbcg_vd1) <- "Timepoint"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = "Week_4")

# Generate a reference table containing the human gene sets to which human genes are assigned. Genes appear more than once depending on the number of gene sets they are associated with. We are only looking in the MSigDB C7 category (immunologic signatures).
msig_<- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
sets <- msig_%>% filter(
  str_detect(gs_name, ("_T_CELL"))|str_detect(gs_name, ("CYTOTOXIC")))
sets <- sets %>% split(x = .$gene_symbol, f = .$gs_name)

# Break into a new object for each PTID #

Idents(ivbcg_vd1) <- "PTID"
dg0d <- ivbcg_vd1 %>% subset(idents = "DG0D")
dg0h <- ivbcg_vd1 %>% subset(idents = "DG0H")
v113 <- ivbcg_vd1 %>% subset(idents = "V113")
v015 <- ivbcg_vd1 %>% subset(idents = "V015")


## DG0D ##
Idents(dg0d) <- "Tissue"
DimPlot(dg0d, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0d <- wilcoxauc(dg0d, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.dg0d <- deframe(bal_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.dg0d <- fgsea(sets, stats = bal.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.dg0d <- gsea.bal.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.dg0d <- deframe(pbmc_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.dg0d <- fgsea(sets, stats = pbmc.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.dg0d <- gsea.pbmc.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DG0H ##
Idents(dg0h) <- "Tissue"
DimPlot(dg0h, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0h <- wilcoxauc(dg0h, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.dg0h <- deframe(bal_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.dg0h <- fgsea(sets, stats = bal.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.dg0h <- gsea.bal.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.dg0h <- deframe(pbmc_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.dg0h <- fgsea(sets, stats = pbmc.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.dg0h <- gsea.pbmc.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V113 ##
Idents(v113) <- "Tissue"
DimPlot(v113, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v113 <- wilcoxauc(v113, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.v113 <- deframe(bal_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.v113 <- fgsea(sets, stats = bal.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.v113 <- gsea.bal.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.v113 <- deframe(pbmc_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.v113 <- fgsea(sets, stats = pbmc.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.v113 <- gsea.pbmc.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V015 ##
Idents(v015) <- "Tissue"
DimPlot(v015, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v015 <- wilcoxauc(v015, 'Tissue')

## T cell  gene sets ##
# BAL
bal_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "BAL") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
bal.ranks.v015 <- deframe(bal_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.bal.v015 <- fgsea(sets, stats = bal.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.bal.v015 <- gsea.bal.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# PBMC
pbmc_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "PBMC") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
pbmc.ranks.v015 <- deframe(pbmc_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.pbmc.v015 <- fgsea(sets, stats = pbmc.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.pbmc.v015 <- gsea.pbmc.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DOT PLOT SEPARATED BY PTID ##
gsea.bal.dg0d$tissue <- "BAL"
gsea.bal.dg0d$ptid <- "DG0D"
gsea.pbmc.dg0d$tissue <- "PBMC"
gsea.pbmc.dg0d$ptid <- "DG0D"


gsea.bal.dg0h$tissue <- "BAL"
gsea.bal.dg0h$ptid <- "DG0H"
gsea.pbmc.dg0h$tissue <- "PBMC"
gsea.pbmc.dg0h$ptid <- "DG0H"


gsea.bal.v113$tissue <- "BAL"
gsea.bal.v113$ptid <- "V113"
gsea.pbmc.v113$tissue <- "PBMC"
gsea.pbmc.v113$ptid <- "V113"


gsea.bal.v015$tissue <- "BAL"
gsea.bal.v015$ptid <- "V015"
gsea.pbmc.v015$tissue <- "PBMC"
gsea.pbmc.v015$ptid <- "V015"


gsea.all <- bind_rows(gsea.bal.dg0d, 
                      gsea.pbmc.dg0d, 
                      gsea.bal.dg0h, 
                      gsea.pbmc.dg0h, 
                      gsea.bal.v113, 
                      gsea.pbmc.v113, 
                      gsea.bal.v015, 
                      gsea.pbmc.v015)

gsea.all$id <- paste(gsea.all$ptid, gsea.all$tissue, sep = "_")


gsea.all %>%
  ggplot(aes(x = tissue, y = pathway, color = NES, facet_by = "ptid")) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)

gsea.all %>%
  ggplot(aes(x = tissue, y = pathway, color = NES, size = padj < 0.01)) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5))+
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)


ggsave("out/images/gsea_all_faceted.png", width = 16, height = 10, units = c("in"))

gsea.all$pathway <- gsub("GOBP_", "", as.character(gsea.all$pathway))
gsea.all$pathway <- gsub("_", " ", as.character(gsea.all$pathway))

gsea.all %>%
  filter(pathway == "ALPHA BETA T CELL ACTIVATION"|
           pathway == "T CELL PROLIFERATION"|
           pathway == "NK T CELL DIFFERENTIATION"|
           pathway == "ALPHA BETA T CELL DIFFERENTIATION"|
           pathway == "LEUKOCYTE MEDIATED CYTOTOXICITY"|
           pathway == "T CELL MEDIATED IMMUNITY"|
           pathway == "POSITIVE THYMIC T CELL SELECTION"
  ) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05, facet_by = "ptid")) + 
  scale_y_discrete(labels = label_wrap(30)) +
  geom_point() + 
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)




##########
## GSEA ##
##########

ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

# BAL WEEK 4 VS BAL WEEK 8
Idents(ivbcg_vd1) <- "Tissue"
ivbcg_vd1 <- ivbcg_vd1 %>% subset(idents = "BAL")

# Generate a reference table containing the human gene sets to which human genes are assigned. Genes appear more than once depending on the number of gene sets they are associated with. We are only looking in the MSigDB C7 category (immunologic signatures).
msig_<- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
sets <- msig_%>% filter(
  str_detect(gs_name, ("_T_CELL"))|str_detect(gs_name, ("CYTOTOXIC")))
sets <- sets %>% split(x = .$gene_symbol, f = .$gs_name)

# Break into a new object for each PTID #

Idents(ivbcg_vd1) <- "PTID"
dg0d <- ivbcg_vd1 %>% subset(idents = "DG0D")
dg0h <- ivbcg_vd1 %>% subset(idents = "DG0H")
v113 <- ivbcg_vd1 %>% subset(idents = "V113")
v015 <- ivbcg_vd1 %>% subset(idents = "V015")


## DG0D ##
Idents(dg0d) <- "Timepoint"
DimPlot(dg0d, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0d <- wilcoxauc(dg0d, 'Timepoint')

## T cell  gene sets ##
# Week 4 BAL
wk4_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk4.ranks.dg0d <- deframe(wk4_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk4.dg0d <- fgsea(sets, stats = wk4.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk4.dg0d <- gsea.wk4.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 8 BAL
wk8_genes_dg0d <- fgsea.genes.dg0d %>%
  dplyr::filter(group == "Week_8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk8.ranks.dg0d <- deframe(wk8_genes_dg0d) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk8.dg0d <- fgsea(sets, stats = wk8.ranks.dg0d, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk8.dg0d <- gsea.wk8.dg0d %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## DG0H ##
Idents(dg0h) <- "Timepoint"
DimPlot(dg0h, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.dg0h <- wilcoxauc(dg0h, 'Timepoint')

## T cell  gene sets ##
# Week 4 BAL
wk4_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk4.ranks.dg0h <- deframe(wk4_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk4.dg0h <- fgsea(sets, stats = wk4.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk4.dg0h <- gsea.wk4.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 8 BAL
wk8_genes_dg0h <- fgsea.genes.dg0h %>%
  dplyr::filter(group == "Week_8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk8.ranks.dg0h <- deframe(wk8_genes_dg0h) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk8.dg0h <- fgsea(sets, stats = wk8.ranks.dg0h, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk8.dg0h <- gsea.wk8.dg0h %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V113 ##
Idents(v113) <- "Timepoint"
DimPlot(v113, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v113 <- wilcoxauc(v113, 'Timepoint')

## T cell  gene sets ##
# Week 4 BAL
wk4_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk4.ranks.v113 <- deframe(wk4_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk4.v113 <- fgsea(sets, stats = wk4.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk4.v113 <- gsea.wk4.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 8 BAL
wk8_genes_v113 <- fgsea.genes.v113 %>%
  dplyr::filter(group == "Week_8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk8.ranks.v113 <- deframe(wk8_genes_v113) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk8.v113 <- fgsea(sets, stats = wk8.ranks.v113, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk8.v113 <- gsea.wk8.v113 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))


## V015 ##
Idents(v015) <- "Timepoint"
DimPlot(v015, reduction = "umap")

# Generate a table containing the DE genes, with p-values, for each group (wilcoxon rank sum test)
fgsea.genes.v015 <- wilcoxauc(v015, 'Timepoint')

## T cell  gene sets ##
# Week 4 BAL
wk4_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "Week_4") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk4.ranks.v015 <- deframe(wk4_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk4.v015 <- fgsea(sets, stats = wk4.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk4.v015 <- gsea.wk4.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))

# Week 8 BAL
wk8_genes_v015 <- fgsea.genes.v015 %>%
  dplyr::filter(group == "Week_8") %>%
  arrange(desc(auc)) %>% 
  dplyr::select(feature, auc)
wk8.ranks.v015 <- deframe(wk8_genes_v015) ## Convert the table into a vector. Fgsea requires this format to make comparisons.
gsea.wk8.v015 <- fgsea(sets, stats = wk8.ranks.v015, nproc=1) ## Calculate p-values and enrichment scores for the gene sets associated with each cluster
gsea.wk8.v015 <- gsea.wk8.v015 %>% ## Tidy up the table
  as_tibble() %>%
  arrange(desc(NES))





## DOT PLOT SEPARATED BY PTID ##
gsea.wk4.dg0d$timepoint <- "wk4"
gsea.wk4.dg0d$ptid <- "DG0D"
gsea.wk8.dg0d$timepoint <- "wk8"
gsea.wk8.dg0d$ptid <- "DG0D"


gsea.wk4.dg0h$timepoint <- "wk4"
gsea.wk4.dg0h$ptid <- "DG0H"
gsea.wk8.dg0h$timepoint <- "wk8"
gsea.wk8.dg0h$ptid <- "DG0H"


gsea.wk4.v113$timepoint <- "wk4"
gsea.wk4.v113$ptid <- "V113"
gsea.wk8.v113$timepoint <- "wk8"
gsea.wk8.v113$ptid <- "V113"


gsea.wk4.v015$timepoint <- "wk4"
gsea.wk4.v015$ptid <- "V015"
gsea.wk8.v015$timepoint <- "wk8"
gsea.wk8.v015$ptid <- "V015"


gsea.all <- bind_rows(gsea.wk4.dg0d, 
                      gsea.wk8.dg0d, 
                      gsea.wk4.dg0h, 
                      gsea.wk8.dg0h, 
                      gsea.wk4.v113, 
                      gsea.wk8.v113, 
                      gsea.wk4.v015, 
                      gsea.wk8.v015)

gsea.all$id <- paste(gsea.all$ptid, gsea.all$timepoint, sep = "_")


gsea.all %>% filter(padj < 0.05) %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, facet_by = "ptid")) + 
  geom_point() + 
  scale_y_discrete(labels = label_wrap(15)) +
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid)

ggsave("out/images/gsea_BAL_wk4_v_wk8.png", width = 16, height = 10, units = c("in"))

gsea.all$pathway <- gsub("GOBP_", "", as.character(gsea.all$pathway))
gsea.all$pathway <- gsub("_", " ", as.character(gsea.all$pathway))
gsea.all$timepoint <- gsub("wk4", "Wk.4", as.character(gsea.all$timepoint))
gsea.all$timepoint <- gsub("wk8", "Wk.8", as.character(gsea.all$timepoint))

gsea.all %>%
  filter(pathway == "T CELL SELECTION"|
           pathway == "ALPHA BETA T CELL DIFFERENTIATION"|
           pathway == "LEUKOCYTE MEDIATED CYTOTOXICITY") %>%
  ggplot(aes(x = timepoint, y = pathway, color = NES, size = padj < 0.05, facet_by = "ptid")) + 
  scale_y_discrete(labels = label_wrap(30)) +
  geom_point() + 
  scale_size_manual(values = c(2, 5)) +
  scale_color_gradient2(low = "#00A2DA", mid = "#FEFFB2", high = "#E54500") +
  facet_grid(. ~ ptid) +
  theme_bw()

