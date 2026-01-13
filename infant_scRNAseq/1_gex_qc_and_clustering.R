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
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc.gd", min.cells = 3, min.features = 200)
# Confirm that the object was initialized correctly. "pbmc" should return "An object of class Seurat" with a brief description of the object.
pbmc

initial_count <- ncol(pbmc) # record initial cell count

## Begin the pre-processing workflow.
# Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize RNA features, RNA counts, and % mitochondrial DNA as a violin plot.
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Use the violin plots to define the aberrent values. Exclude the aberrant values.
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10)

gex_filters <- ncol(pbmc) # record cell count after QC

## Normalize the data. The following is the recommended default method.
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

## Select the 2000 most variable features in the dataset.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes.
top10 <- head(VariableFeatures(pbmc), 10)


## Scale the data from each gene so that the mean expression is 0 and the variance is 1.
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

## Perform linear dimensional reduction using the 2000 most variable features defined previously.
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Graph the output of dimensional reduction on a 2D scatter plot.
DimPlot(pbmc, reduction = "pca")

## Determine the dimensionality of the dataset using an elbow plot.
ElbowPlot(pbmc) # 10 dimensions

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
  slice_max(n = 2, order_by = avg_log2FC)

# Generate a heatmap showing the differentially expressed genes.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# Overlay the clusters with expression levels for genes of interest.
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
  subset(idents = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15))
unique(pbmc$seurat_clusters)

tcells <- ncol(pbmc) # record the number of cells identified as T cells


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

# Save the object
saveRDS(pbmc, file = "GEX/Out/pbmc_tcells_2024.rds")

