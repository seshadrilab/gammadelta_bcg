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


###############################
## FILTER TO GD T CELLS ONLY ##
###############################

# Read in the Seurat object and quality-checked TCRs
gex <- readRDS(file = "GEX/Out/pbmc_tcells_2024.rds")
tcr <- read.csv(file = "TCR/out/tcr_QCd_productive.csv")

# Here we are defining gamma-delta T cells as cells expressing at least one productive, quality-checked gamma or delta chain
# Cells not meeting these criteria are filtered out
gd_bc <- tcr$barcode
gd_bc <- paste(gd_bc, collapse = "|")
gd_bc <- str_detect(colnames(gex), gd_bc)
gex@meta.data$Celltype[gd_bc] <- "gd"
gex@meta.data$Celltype[!gd_bc] <- "other"
Idents(gex) <- "Celltype"
gex <- gex %>%
  subset(idents = c("gd"))

# After filtering to GD T cells only, re-cluster the data
gex <- NormalizeData(gex, normalization.method = "LogNormalize", scale.factor = 10000)
gex <- FindVariableFeatures(gex, selection.method = "vst", nfeatures = 2000)
gex <- ScaleData(gex)
gex <- RunPCA(gex, features = VariableFeatures(object = gex))
ElbowPlot(gex)
gex <- FindNeighbors(gex, dims = 1:10)
gex <- FindClusters(gex, resolution = 0.5)
gex <- RunUMAP(gex, dims = 1:10)
DimPlot(gex, reduction = "umap")

# Save this object so it can be easily loaded back in without needing to re-run the scripts.
saveRDS(gex, file = "GEX/Out/pbmc_2024_gd.rds")

# Generate the plot of cell counts as quality checks are applied
# Initial_count, gex_filters_ and tcells are generated in script 1 GEX_QC_and_clustering
gd_tcells <- ncol(gex)
qc_table <- data.frame(initial_count, gex_filters, tcells, gd_tcells)
qc_table <- qc_table %>% pivot_longer(cols = 1:4, names_to = "step", values_to = "count")
order <- c("initial_count", "gex_filters", "tcells", "gd_tcells")
ggplot(qc_table, aes(x = factor(step, order), y = count)) +
  geom_bar(stat = "identity")
ggsave("SuppFig1A.pdf")


####################################
## CELL CLUSTERING AND ANNOTATION ##
####################################

# Read in the TCR groups so we can use them for annotation
tcr_grouped <- read.csv("TCR/Out/tcr_grouped.csv")

# Enumerate the number of samples from each ptid (sample donor)
md <- gex@meta.data
temp <- table(md$ptid) %>% as.data.frame()
mean(temp$Freq)

# Find the differentially-expressed genes across each cluster
pbmc.markers <- FindAllMarkers(gex, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Generate a heatmap showing the differentially expressed genes.
pbmc.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(gex, features = top10$gene) + NoLegend()

# Save the cluster markers
write.csv(pbmc.markers, file = "GEX/Out/pbmc.markers_2024.csv")


## Identify the GD TCR subsets in each cluster
# Assign a TCR group to each cell barcode
group_bc <- tcr_grouped %>% select(c(barcode, tcr_group, paired_tcr))
group_bc2 <- group_bc[, -1]
rownames(group_bc2) <- group_bc$barcode
gex <- AddMetaData(gex, metadata = group_bc2)

# Visualize the TCR groups on the UMAP
DimPlot(gex, split.by = "tcr_group")

# Enumerate the cell barcodes in each TCR group
md <- gex@meta.data
tcr_group_table <- table(md$tcr_group) %>% as.data.frame()
total <- sum(tcr_group_table$Freq)
tcr_group_table <- tcr_group_table %>% mutate(pct = Freq / total * 100)

## Generate a heat map of key genes
# Genes have been separated into groups for convenience
Idents(gex) <- "seurat_clusters"
goi <- pbmc.markers %>% filter(
  gene == "SELL" | # High in naive-like cells
  gene == "CCR7" |
  gene == "TCF7" |
    
  gene == "HLA-DPA1" | # Activation
  gene == "HLA-DPB1" |
  gene == "HLA-DRB1" |
  gene == "TIGIT" |
  gene == "CD69" |
  gene == "PDCD1" |
    
  gene == "STAT1" | # Pro-inflammatory or Th17-like
  gene == "IFNG" |
  gene == "TNF" |
  gene == "TBX21" |
  gene == "XCL1" |
  gene == "MX1" |
  gene == "RORA" |
  gene == "RORC" |
  gene == "CCR6" |

  gene == "CST7" | # Cytotoxicity
  gene == "GZMA" |
  gene == "GZMK" |
  gene == "NKG7" |
  gene == "PRF1" |
  gene == "GNLY" |
  gene == "CD8A" |
  gene == "GZMM" |
  gene == "GZMB" |
  gene == "GZMH" |
  gene == "NCAM1")

# Set the preferred order for the genes in the heat map
goi_order <- c(
  "SELL",
  "CCR7",
  "TCF7",
  "CD69",
  "HLA-DPA1",
  "HLA-DPB1",
  "HLA-DRB1",
  "TIGIT",
  "TBX21",
  "STAT1",
  "MX1",
  "RORA",
  "GZMK",
  "GZMA",
  "CST7",
  "PRF1",
  "GZMM",
  "GNLY",
  "CD8A",
  "GZMB",
  "GZMH",
  "NCAM1"
)
goi <- goi %>% arrange(match(gene, goi_order))

# Average the expression of each gene across the cells in each cluster
avgexp <- AverageExpression(gex, assays = "RNA", return.seurat = T)

# Visualize the heat map
DoHeatmap(avgexp,
  features = goi$gene,
  raster = FALSE,
  draw.lines = F,
  group.bar.height = 0.06,
  label = T
) +
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) +
  theme(text = element_text(size = 10)) +
  guides(color = F)

# Manually annotate and merge the clusters according to their gene expression
new.cluster.ids <- c(
  "Naive Vd1", "Naive Vd1", "Naive Vd1", "GZMB Vd2", "CD69 Vd2",
  "CD69 Vd2", "Vd1 Teff", "Mixed 1", "Mixed 2",
  "Mixed 3", "STAT1-hi Vd2", "STAT1-hi Vd1"
)
names(new.cluster.ids) <- levels(gex)
gex <- RenameIdents(gex, new.cluster.ids)

# Confirm that the cluster IDs have been reassigned
DimPlot(gex,
  reduction = "umap",
  label = TRUE,
  label.color = "black",
  pt.size = 0.75,
  label.size = 4
)

# Set color palettes
my_cols9 <- viridis(9, begin = 0.15)
my_cols3 <- viridis(3, end = 0.9)

# Visualize the UMAPs
DimPlot(gex,
  reduction = "umap",
  label = TRUE,
  label.color = "black",
  pt.size = 0.75,
  label.size = 4,
  cols = my_cols9
)

new.cluster.ids <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
names(new.cluster.ids) <- levels(gex)
gex <- RenameIdents(gex, new.cluster.ids)

## FIGURE 1D
avgexp <- AverageExpression(gex, assays = "RNA", return.seurat = T)
DoHeatmap(avgexp,
  features = goi$gene,
  raster = FALSE,
  draw.lines = F,
  group.colors = my_cols9,
  group.bar.height = 0.06,
  label = F
) +
  scale_fill_gradientn(colors = c("#00A2DA", "#FEFFB2", "#E54500")) +
  theme(text = element_text(size = 10)) +
  guides(color = F)
ggsave("GEX/Out/Fig1D.pdf", width = 4, height = 4, units = "in")


## FIGURE 1A
DimPlot(gex,
  reduction = "umap",
  label = TRUE,
  label.color = "black",
  pt.size = 0.75,
  label.size = 4,
  cols = my_cols9
)
ggsave("GEX/Out/Fig1A.pdf", width = 5, height = 4, units = "in")

## FIGURE 1B
Idents(gex) <- "tcr_group"
gex %>%
  subset(idents = c("Vd1", "Vg9Vd2", "Vd3")) %>%
  DimPlot(reduction = "umap", cols = c("black", "black", "black"), split.by = "tcr_group")
ggsave("GEX/Out/Fig1B.pdf", width = 5.5, height = 2, units = "in")


###########################
## COMBINE WITH TCR DATA ##
###########################

tcr_grouped <- read.csv("TCR/Out/tcr_grouped.csv")
tcr <- read.csv("TCR/Out/tcr_combined.csv")

# Count how many times each TCR clonotype (as defined by paired-chain AA sequence) appears per ptid (sampe donor).
# Identify expanded TCRs
counts <- tcr %>%
  group_by() %>%
  pull(CTnt) %>%
  table() %>%
  as.data.frame()

counts <- tcr %>%
  group_by(ptid, CTaa) %>%
  summarize(Freq = n())

counts$expanded <- ifelse(counts$Freq > 1, "exp", "singlet")

# Quantify the percentage of expanded clonotypes per TCR group
tcr_grouped <- tcr_grouped %>% full_join(counts, by = c("CTaa" = "CTaa", "ptid" = "ptid"))
temp <- tcr_grouped %>%
  subset(select = c("CTaa", "tcr_group", "expanded")) %>%
  distinct()

group_counts <- temp %>%
  group_by(tcr_group, expanded) %>%
  summarize(Freq = n())


# Merge the counts with the original table of combined TCRs
tcr <- tcr %>% full_join(counts, by = c("CTaa" = "CTaa", "ptid" = "ptid"))

# Add to seurat object metadata
tcr_exp <- tcr %>% subset(select = c("barcode", "expanded"))
rownames(tcr_exp) <- tcr_exp$barcode
tcr_exp <- tcr_exp %>% subset(select = -c(barcode))

gex <- AddMetaData(gex, metadata = tcr_exp)


# Visualize the location of the expanded clonotypes
## FIGURE 1C
Idents(gex) <- "tcr_group"
gex_temp <- gex %>% subset(idents = c("Vd1", "Vd3", "Vg9Vd2"))

Idents(gex_temp) <- "expanded"
gex_temp %>%
  subset(idents = c("singlet", "exp")) %>%
  DimPlot(
    reduction = "umap", cols = c("lightgrey", "#440154FF"),
    order = "exp", split.by = "tcr_group"
  )
ggsave("GEX/Out/Fig1C.pdf", width = 7, height = 2.5, units = "in")
