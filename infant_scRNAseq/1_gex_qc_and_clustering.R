
## Load the necessary packages.
library(dplyr)
library(Seurat)
library(patchwork)
library(presto)
library(Rcpp)
library(msigdbr)
library(devtools)
library(fgsea)
library(ggplot2)
library(tibble)
library(stringr)
library(viridis)
library(ggpubr)


## Read in the dataset.
# Confirm that the directory is set correctly and that the barcodes, genes, and matrix files are returned.
data_dir <- getwd()
list.files(data_dir)
# Convert the cellranger output to a UMI matrix using the function Read10x(). The UMI matrix describes the number of molecules for each gene that are detected for each cell.
pbmc.data <- Read10X(data.dir = data_dir)

## Initialize the Seurat Object. 
# min.cells = Include features detected in at least this many cells. 3 is the Satija lab recommended value.
# min.features = Include cells where at least this many features are detected. 200 is the Satija lab recommended value.
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc.gd", min.cells = 3, min.features = 200)
# Confirm that the object was initialized correctly. "pbmc" should return "An object of class Seurat" with a brief description of the object.
pbmc

initial_count <- ncol(pbmc)

## Begin the pre-processing workflow.
# Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern= "^MT-")
# The Seurat Object also contains the automatically calculated metrics nCount_RNA and nFeature_RNA. These describe the total RNA count and the total number of unique features (genes) for each cell.
# Visualize all three of these QC metrics as a violin plot.
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Use the violin plots to define the aberrent values. Exclude the aberrant values.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10)

gex_filters <- ncol(pbmc)

## After removing the unwanted cells, normalize the data. The following is the recommended default method.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Select the most variable features in the dataset. By default, return 2000 features.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes.
top10 <- head(VariableFeatures(pbmc), 10)


## Scale the data from each gene so that the mean expression is 0 and the variance is 1. This gives equal weight to each gene in the downstream analysis and prevents highly-expressed genes from dominating. Scaling is also necessary if you would like to generate a heat map.
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Perform linear dimensional reduction using the 2000 most variable features defined previously.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

## Visualize the PCA in a few different ways.
# Generate a list of the genes at the extreme ends of the first 5 principal component vectors.
print(pbmc[["pca"]], dims = 1:5, nfeatures = 3)
# Visualize the top genes associated with principal components 1 through 5.
VizDimLoadings(pbmc, dims = 1:5, reduction = "pca")
# Graph the output of dimensional reduction on a 2D scatter plot.
DimPlot(pbmc, reduction = "pca")
# Draw a heatmap focusing on each principal component, with cells along the x axis and genes along the y axis. Set the number of cells at a low number (like 500) to visualize the extreme ends of each component and speed up the plotting.
DimHeatmap(pbmc, dims = 1:20, cells = 500, balanced = TRUE)

## Determine the dimensionality of the dataset (Elbow plot).
ElbowPlot(pbmc)

## Cluster the cells.
# Construct a K nearest-neighbor graph.
pbmc <- FindNeighbors(pbmc, dims = 1:10)
# Find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

## Run non-linear dimensional reduction to map the cells in two-dimensional space.
pbmc <- RunUMAP(pbmc, dims = 1:10)
# Visualize the clusters
DimPlot(pbmc, reduction = "umap")

## Identify what the clusters represent.
# Find the markers that distinguish each cluster from the remaining clusters.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)
# Generate a heatmap showing the differentially expressed genes.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
# Overlay the clusters with expression levels for a particular gene or genes.
FeaturePlot(pbmc, features = c("TRDV1", "TRDV2", "TRDV3", "SELL", "IL2", "IFNG", "TNF", "RORC"))
FeaturePlot(pbmc, features = c("TRDV1", "TRDV2", "LAG3", "PDCD1", "HAVCR2", "CTLA4", "CD244", "CD160"))
# Visualize the expression level of a particular gene(s) across all clusters.
VlnPlot(pbmc, features = c("TRDV1", "TRDV2", "MS4A1"))

## Using these data, we can classify each cluster.
new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

## Exclude all cells assigned to cluster 14. These are B cells and not of interest to us here.
pbmc <- pbmc %>%
  subset(idents = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,15))
unique(pbmc$seurat_clusters)

tcells <- ncol(pbmc)


## After excluding the B cells, re-cluster the data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")

# Save this object so it can be easily loaded back in without needing to re-run the scripts.
saveRDS(pbmc, file = "GEX/Out/pbmc_tcells_2024.rds")
pbmc <- readRDS(file = "GEX/Out/pbmc_tcells_2024.rds")










## Next, identify what the clusters represent.
# Find the markers that distinguish each cluster from the remaining clusters.
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n=2, order_by = avg_log2FC)

# Generate a heatmap showing the differentially expressed genes.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n=10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
write.csv(pbmc.markers, file = "GEX/Out/pbmc.markers_2024.csv")
pbmc.markers <- read.csv(file = "GEX/Out/pbmc.markers_2024.csv")








##################################################################
## TO START WITH THE RDS FILE -- GD T cells only, 10 dimensions ##
##################################################################
pbmc_sct <- readRDS(file = "GEX/Out/pbmc_sct_dim10.rds")
pbmc.markers <- read.csv(file = "GEX/Out/pbmc.markers_sct_dim10.csv")

# Visualize the expression level of a particular gene(s) across all clusters.
VlnPlot(pbmc_sct, features = c("TRDV1", "TRDV2", 'TRDV3'))
FeaturePlot(pbmc_sct, features = c('TRDV1', 'TRDV2', 'TRDV3'))
VlnPlot(pbmc_sct, features = c("IKZF2", "CD4", "CD8A", "CD8B", "PDCD1"))

FeaturePlot(pbmc_sct, pt.size = 0.75, order = T, min.cutoff = 0.5, ncol = 4, 
            features = c("TRDV2",
                         "SELL", 
                         "HLA-DRB1", 
                         "NKG7", 
                         "GZMA", 
                         "GZMK", 
                         "PDCD1",
                         "ITGA1")) &
  scale_color_viridis_c() &
  theme(text = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank())


## Assign a TCR group to each cell barcode
tcr_grouped <- read.csv("TCR/Out/tcr_grouped.csv")
group_bc <- tcr_grouped %>% select(c(barcode, tcr_group, paired_tcr))
group_bc2 <- group_bc[,-1]
rownames(group_bc2) <- group_bc$barcode
pbmc_sct <- AddMetaData(pbmc_sct, metadata = group_bc2)



## Volcano plot of DE genes between Vd1-C5 and other Vd1 clusters

# Annotate the cell barcodes as "pre" or "post"
Idents(inf_vd1) <- "seurat_clusters"
c5 <- inf_vd1 %>% subset(idents = c("5"))
not_c5 <- inf_vd1 %>% subset(idents = (c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11", "12", "13")))
c5 <- c5@meta.data
not_c5 <- not_c5@meta.data
c5$gsea <- "C5"
not_c5$gsea <- "Not_C5" 
anno <- rbind(c5, not_c5) %>% subset(select = c("gsea"))
inf_vd1 <- AddMetaData(inf_vd1, metadata = anno)

# Find DE genes between Pre and Post
Idents(inf_vd1) <- "gsea"
markers <- FindAllMarkers(inf_vd1)

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC> 1 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -1 & markers$p_val_adj < 0.05] <- "DOWN"
markers$p_val_adj[markers$p_val_adj<6.968133e-300] <- as.numeric(2.0e-300)

markers$delabel <- NA
markers$delabel[markers$diffex != "NO"] <- markers$gene[markers$diffexp != "NO"]
markers <- markers %>% filter(cluster == "C5")

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline(yintercept = -log10(0.05), col = "red") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel)+
  ylim(-5, 320)











###############
## HEAT MAPS ##
###############

Idents(pbmc_sct) <- "seurat_clusters"
goi <- pbmc.markers %>% filter(gene == 'SELL'|
                                 gene == 'CCR7' |
                                 gene == 'TCF7' |
                                 gene == 'HLA-DPA1' |
                                 gene == 'HLA-DPB1' |
                                 gene == 'HLA-DRB1' |
                                 gene == 'TIGIT' |
                                 gene == 'CD69' |
                                 gene == 'TOX' |
                                 gene == 'PDCD1' |
                                 gene == 'STAT1' |
                                 gene == 'IFNG' |
                                 gene == 'TNF' |
                                 gene == 'TBX21' |
                                 gene == 'XCL1' |
                                 gene == 'MX1' |
                                 gene == 'RORA' |
                                 gene == 'RORC' |
                                 gene == 'CCR6' |
                                 
                                 gene == 'CST7' |
                                 gene == 'GZMA' |
                                 gene == 'GZMK' |
                                 gene == 'NKG7' |
                                 gene == 'PRF1' |
                                 gene == 'GNLY' |
                                 
                                 gene == 'CD8A' |
                                 gene == 'GZMM' |
                                 gene == 'GZMB' |
                                 gene == 'GZMH' |
                                 gene == 'NCAM1')
goi_order <- c('SELL',
               'CCR7',
               'TCF7',
               'CD69',
               
               'HLA-DPA1',
               'HLA-DPB1',
               'HLA-DRB1',
               'TIGIT',
               'TOX',
               'PDCD1',
               'TBX21',
               'STAT1',
               'IFNG',
               'TNF',
               'XCL1',
               'MX1',
               
               'RORA',
               'RORC',
               'CCR6',
               
               'CST7',
               'GZMA',
               'GZMK',
               'NKG7',
               'PRF1',
               'GNLY',
               
               'CD8A',
               'GZMM',
               'GZMB',
               'GZMH',
               'NCAM1')

goi_order <- c('SELL',
               'CCR7',
               'TCF7',
               'CD69',
               
               'HLA-DPA1',
               'HLA-DPB1',
               'HLA-DRB1',
               'TIGIT',
               'TBX21',
               'STAT1',
               'MX1',
               'RORA',
               
               'GZMK',
               'GZMA',
               'CST7',
               'PRF1',
               'GZMM',
               'GNLY',
               
               'CD8A',
               'GZMB',
               'GZMH',
               'NCAM1')
goi <- goi %>% arrange(match(gene, goi_order))
avgexp <- AverageExpression(pbmc_sct, assays = "SCT", return.seurat = T)
new.cluster.ids <- c('Naive Vd1', 'Naive Vd1', 'Memory Vd1', 'CD69-hi Vd2', 'Mixed cytotoxic',
                     'Activated Vd1', 'Naive Vd1', 'Mixed CD69-hi', 'CD69-lo Vd2',
                     'Naive Vd1', 'Memory Vd1', 'STAT1-hi Vd1', 'STAT1-hi Vd2', 'RORA-hi Vd2')
names(new.cluster.ids) <- levels(pbmc_sct)
pbmc_sct <- RenameIdents(pbmc_sct, new.cluster.ids)

# Set color palettes
my_cols10 <- viridis(10, begin = 0.15)
my_cols3 <- viridis(3, end = 0.9)


# Visualize the UMAPs
DimPlot(pbmc_sct, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.75,
        label.size = 4,
        cols = my_cols10)

new.cluster.ids <- c('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
names(new.cluster.ids) <- levels(pbmc_sct)
pbmc_sct <- RenameIdents(pbmc_sct, new.cluster.ids)

a <- DimPlot(pbmc_sct, reduction = "umap", 
        label = TRUE,
        label.color = 'black',
        pt.size = 0.75,
        label.size = 5,
        cols = my_cols10)

Idents(pbmc_sct) <- "tcr_group"
b <- pbmc_sct %>% subset(idents = c("Vd1", "Vg9Vd2", "Vd3")) %>%
  DimPlot(reduction = "umap", cols = my_cols3, split.by = "tcr_group")


ggarrange(a,b, ncol = 2, nrow = 1, widths = c(0.5, 1))

# Visualize the heat map
DoHeatmap(avgexp, features = goi$gene, 
          raster = FALSE, 
          draw.lines = F,
          group.colors = my_cols10,
          group.bar.height = 0.06,
          label = F) + 
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) + 
  theme(text = element_text(size = 10)) +
  guides(color = F)



###########################
## COMBINE WITH TCR DATA ##
###########################

tcr_infant <- read.csv("TCR/Out/tcr_combined_annotated.csv")

# Remove the unnecessary columns and rename for clarity
tcr_md <- tcr_infant %>% subset(select = c('barcode', 'cdr3_aa1', 'cdr3_aa2', 'TRD', 'TRG', 'paired_tcr', 'Freq')) %>%
  rename(cdr3_delta = cdr3_aa1,
         cdr3_gamma = cdr3_aa2,
         trd = TRD,
         trg = TRG,
         clono_count = Freq) %>%
  mutate(count = ifelse(clono_count == 1, 'singlet', 'expanded'))


# Add to seurat object metadata
rownames(tcr_md) <- tcr_md$barcode
tcr_md <- tcr_md %>% subset(select = -c(barcode))
pbmc_sct <- AddMetaData(pbmc_sct, metadata = tcr_md)

Idents(pbmc_sct) <- 'count'

pbmc_sct %>% subset(idents = c('singlet', 'expanded')) %>% 
  DimPlot(reduction = "umap", cols = c('lightgrey', '#440154FF'),
          order = 'expanded', split.by = "tcr_group")

