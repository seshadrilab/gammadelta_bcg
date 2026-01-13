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


# Read in the Seurat objects
dg0d <- readRDS(file = "out/seurat_objects/dg0d_analysis.rds")
dg0h <- readRDS(file = "out/seurat_objects/dg0h_analysis.rds")
v113 <- readRDS(file = "out/seurat_objects/v113_analysis.rds") 
v015 <- readRDS(file = "out/seurat_objects/v015_analysis.rds")

# Set the default assay to RNA and identity to CellType
DefaultAssay(dg0d) <- "RNA"
DefaultAssay(dg0h) <- "RNA"
DefaultAssay(v113) <- "RNA"
DefaultAssay(v015) <- "RNA"

Idents(dg0d) <- "Celltype"
Idents(dg0h) <- "Celltype"
Idents(v113) <- "Celltype"
Idents(v015) <- "Celltype"

# Add cell IDs for each donor
dg0d <- RenameCells(dg0d, add.cell.id = "dg0d")
dg0h <- RenameCells(dg0h, add.cell.id = "dg0h")
v113 <- RenameCells(v113, add.cell.id = "v113")
v015 <- RenameCells(v015, add.cell.id = "v015")

# Filter each Seurat object to GD T cells only
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

# Cluster the integrated data
DefaultAssay(ivbcg) <- "integrated"
ivbcg <- ScaleData(ivbcg, verbose = FALSE)
ivbcg <- RunPCA(ivbcg, npcs = 30, verbose = FALSE)
ivbcg <- RunUMAP(ivbcg, reduction = "pca", dims = 1:12)
ivbcg <- FindNeighbors(ivbcg, reduction = "pca", dims = 1:12)
ivbcg <- FindClusters(ivbcg, resolution = 0.4)

# Visualise the clusters labeled by Cluster ID, donor (ptid), tissue, and timepoint
DimPlot(ivbcg, reduction = "umap")
DimPlot(ivbcg, reduction = "umap", group.by = "PTID")
DimPlot(ivbcg, reduction = "umap", group.by = "Tissue")
DimPlot(ivbcg, reduction = "umap", group.by = "Timepoint")

# Visualize the GD T cell subsets (V-gene usage) across each cluster
FeaturePlot(ivbcg, pt.size = 0.5, order = T, min.cutoff = 1.5, ncol = 3, 
            features = c("TRDV1", 
                         "TRDV2", 
                         "TRDV3")) &
  scale_color_viridis_c() &
  theme(text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())

# Annotate the GD subsets according to their TCR gene usage
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


# Save the Seurat Objects
saveRDS(ivbcg, file = "out/seurat_objects/ivbcg_integrated_2024.rds")
saveRDS(ivbcg_vd1, file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")



###############
## Heat maps ##
###############

# Read in the Seurat Objects
ivbcg <- readRDS(file = "out/seurat_objects/ivbcg_integrated_2024.rds")
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds") 


# Set color palettes
my_cols12 <- viridis(12, begin = 0.15)
my_cols6 <- viridis(6, end = 0.9)
my_cols4 <- viridis(4, end = 0.9)
my_cols3 <- viridis(3, end = 0.9)
my_cols2 <- viridis(2, begin = 0.3, end = 0.9)

## Visualize the UMAPs according to:
# Cluster ID
Idents(ivbcg) <- "seurat_clusters"
DimPlot(ivbcg, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.1,
        label.size = 4,
        cols = my_cols12)
ggsave('out/images/fig2b.pdf', width = 4, height = 3, units = 'in')

# Tissue
Idents(ivbcg) <- "Tissue"
DimPlot(ivbcg, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.1,
        label.size = 5,
        cols = my_cols2)
ggsave('out/images/fig2c.pdf', width = 4.3, height = 3, units = 'in')


# TCR group
Idents(ivbcg) <- "tcr_group"
ivbcg %>% subset(idents = c("Vd1", "Vg9Vd2", "Vd3")) %>%
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "tcr_group",
          pt.size = 0.01)
ggsave('out/images/fig2d.pdf', width = 9, height = 3, units = 'in')

# Generate seurat objects containing only BAL or PBMC cells
Idents(ivbcg) <- "Tissue"
ivbcg_bal <- ivbcg %>% subset(idents = c('BAL'))
ivbcg_pbmc <- ivbcg %>% subset(idents = c('PBMC'))

# Visualize the PBMC-derived clusters at each time point
ivbcg_pbmc %>% 
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "Timepoint",
          pt.size = 0.01)
ggsave('out/images/fig2e.pdf', width = 9, height = 3, units = 'in')

# Visualize the BAL-derived cells at each time point
ivbcg_bal %>% 
  DimPlot(reduction = "umap", 
          cols = c('black', 'black', 'black'), 
          split.by = "Timepoint",
          pt.size = 0.01)
ggsave('out/images/fig2f.pdf', width = 6.5, height = 3, units = 'in')



# Generate a heat map of genes of interest in each cluster
# Genes are grouped by their associated phenotype/function
marker_list <- c("CCR7", "TCF7", "ID3", # Naive-like
                 "TNFRSF9", "NR4A1", "NR4A2", "NR4A3", "ICOS", "ICAM1", "CD69", "MAMU-DRB1", "TIGIT", "ITGAX", "CD40LG", "TNFRSF4", # Activation
                 "TOX", "PDCD1", "CD96", "LAG3", # Exhaustion
                 "XCL1", "MX1", "STAT1", "IFNG", "TNF", # Type 1 function
                 "RORA", "RORC", "IL21", "CCR6", # Type 3 function
                 "ITGB1", "NCAM1", "TNFSF8","TNFRSF10B", "GZMK", "GZMM", "NKG7", "CST7", "GZMB", "GZMA", "PRF1", "TNFRSF18", # Cytotoxicity
                 "ITGA1", "ITGB7", "ITGAE" # Tissue residency
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


# Generate the heat map in dot plot form
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

ggsave('out/images/fig2g.pdf', units = c("in"), width = 11.5, height = 4)



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

# Identify the frequency of each tcr group in each cluster
ivbcg@meta.data %>% filter(tcr_group != "Gamma_only") %>% ggplot(aes(x=seurat_clusters, fill=tcr_group)) + geom_bar(position = "fill") + scale_fill_viridis_d()

# Export the metadata from the Seurat object and summarize the proportion of cells in each cluster over time
ivbcg_md <- ivbcg@meta.data

summary <- ivbcg_md %>% group_by(Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))

# Plot the proportion of PBMC cells in each cluster over time
# Clusters 6-11 can be grouped together as they are less abundant than clusters 0-5
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


# Generate plots showing the frequency of clusters 0, 3, and 4 over time
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


# Export the metadata from the Seurat object and summarize the proportion of cells in each cluster over time
ivbcg_md <- ivbcg@meta.data

summary <- ivbcg_md %>% group_by(Tissue, Timepoint, seurat_clusters) %>%
  summarize(count = n()) %>%
  mutate(freq = count/sum(count))


# Plot the proportion of BAL cells in each cluster over time
# Clusters 3-11 can be grouped together as they are less abundant than clusters 0-2
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


# Generate plots showing the frequency of clusters 0 and 2 over time
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



############################################
## DE GENES ACROSS TISSUES AND TIMEPOINTS ##
############################################

# Generate a new metadata column containing SampleType: Tissue_Timepoint
ivbcg_vd1 <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
ivbcg_vd1_md <- ivbcg_vd1@meta.data
ivbcg_vd1_md$SampleType <- paste(ivbcg_vd1_md$Tissue, ivbcg_vd1_md$Timepoint, sep = "_")
ST <- ivbcg_vd1_md %>% select(c('SampleType'))
ivbcg_vd1 <- AddMetaData(ivbcg_vd1, metadata = ST)

# Find the differentially expressed genes across all sample types
Idents(ivbcg_vd1) <- "SampleType"
markers_all <- FindAllMarkers(ivbcg_vd1)

# Generate a list of the top 30 genes in each cluster, excluding genes with no symbol and alpha-beta TCR genes
markers_all %>% 
  group_by(cluster) %>%
  top_n(n=30, wt = avg_log2FC) -> top30_all
top30_all <- top30_all %>% filter(!grepl('ENSMMUG', gene)) %>%
  filter(!grepl('TRAV', gene)) %>%
  filter(!grepl('TRBV', gene))

# Visualize the top 30 genes in a heat map
avgexp <- AverageExpression(ivbcg_vd1, assays = "integrated", return.seurat = T) # Average the expression across each cluster to make the heat map easier to interpret.
heatmap_top30 <- DoHeatmap(avgexp, features = top30_all$gene, raster = FALSE, draw.lines = F) + 
  NoLegend() + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) # Generate a heatmap of the averaged expression, using custom colors.
heatmap_top30
write.csv(markers_all, 'out/markers_vd1_all.csv')
write.csv(top30_all, 'out/top30_vd1_all.csv')


# Visualize a heat map containing the top genes of interest
filtered <- read.csv('out/filtered_samptype_markers2.csv')
filtered <- filtered$gene

heatmap_filtered <- DoHeatmap(avgexp, features = filtered, raster = FALSE, draw.lines = F) + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 11))

heatmap_filtered
ggsave('out/images/fig2j.pdf', width = 4, height = 6.5, units = 'in')
