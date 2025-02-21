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
library(tidyr)


###########################
## Set up Seurat Objects ##
###########################

#Read in the filtered barcode matrices and initialize the Seurat objects.
dg0d <- Read10X("Tube1_multi_vdjgenes_DG0D/per_sample_outs/results/count/sample_filtered_feature_bc_matrix")
dg0h <- Read10X("Tube1_multi_vdjgenes_DG0H/per_sample_outs/results/count/sample_filtered_feature_bc_matrix")
v015 <- Read10X("Tube1_multi_vdjgenes_V015/per_sample_outs/results/count/sample_filtered_feature_bc_matrix")
v113 <- Read10X("Tube1_multi_vdjgenes_V113/per_sample_outs/results/count/sample_filtered_feature_bc_matrix")

ptids <- c('dg0d', 'dg0h', 'v015', 'v113') 
initial_count = c(ncol(dg0d$`Gene Expression`), ncol(dg0h$`Gene Expression`), ncol(v015$`Gene Expression`), ncol(v113$`Gene Expression`))

qc_table <- data.frame(ptids, initial_count)

#Extract the ADT and HTO matrices
dg0d_adt <- !grepl("HTO", dg0d$'Antibody Capture'@Dimnames[[1]])
dg0d_adt_counts <- dg0d$'Antibody Capture'[dg0d_adt,]
dg0d_hto <- grepl("HTO", dg0d$'Antibody Capture'@Dimnames[[1]])
dg0d_hto_counts <- dg0d$'Antibody Capture'[dg0d_hto,]

dg0h_adt <- !grepl("HTO", dg0h$'Antibody Capture'@Dimnames[[1]])
dg0h_adt_counts <- dg0h$'Antibody Capture'[dg0h_adt,]
dg0h_hto <- grepl("HTO", dg0h$'Antibody Capture'@Dimnames[[1]])
dg0h_hto_counts <- dg0h$'Antibody Capture'[dg0h_hto,]

v015_adt <- !grepl("HTO", v015$'Antibody Capture'@Dimnames[[1]])
v015_adt_counts <- v015$'Antibody Capture'[v015_adt,]
v015_hto <- grepl("HTO", v015$'Antibody Capture'@Dimnames[[1]])
v015_hto_counts <- v015$'Antibody Capture'[v015_hto,]

v113_adt <- !grepl("HTO", v113$'Antibody Capture'@Dimnames[[1]])
v113_adt_counts <- v113$'Antibody Capture'[v113_adt,]
v113_hto <- grepl("HTO", v113$'Antibody Capture'@Dimnames[[1]])
v113_hto_counts <- v113$'Antibody Capture'[v113_hto,]


# Keep HTOs that appear in at least 50 cells (to prevent demux errors)
nonzero_cells <- rowSums(dg0d_hto_counts>=1)
keep_rows <- rownames(dg0d_hto_counts[(nonzero_cells) > 50, ])
dg0d_hto_counts <- dg0d_hto_counts[keep_rows,]
row.names(dg0d_hto_counts)

nonzero_cells <- rowSums(dg0h_hto_counts>=1)
keep_rows <- rownames(dg0h_hto_counts[(nonzero_cells) > 50, ])
dg0h_hto_counts <- dg0h_hto_counts[keep_rows,]
row.names(dg0h_hto_counts)

nonzero_cells <- rowSums(v113_hto_counts>=1)
keep_rows <- rownames(v113_hto_counts[(nonzero_cells) > 50, ])
v113_hto_counts <- v113_hto_counts[keep_rows,]
row.names(v113_hto_counts)

nonzero_cells <- rowSums(v015_hto_counts>=1)
keep_rows <- rownames(v015_hto_counts[(nonzero_cells) > 50, ])
v015_hto_counts <- v015_hto_counts[keep_rows,]
row.names(v015_hto_counts)


# Keep cells with >2 total HTO expression only (to prevent demux errors)
keep_cells <- colnames(dg0d_hto_counts[, colSums(dg0d_hto_counts) > 2])
dg0d_adt_counts <- dg0d_adt_counts[, keep_cells]
dg0d_hto_counts <- dg0d_hto_counts[, keep_cells]
dg0d_gex_counts <- dg0d$`Gene Expression`[,keep_cells]

keep_cells <- colnames(dg0h_hto_counts[, colSums(dg0h_hto_counts) > 2])
dg0h_adt_counts <- dg0h_adt_counts[, keep_cells]
dg0h_hto_counts <- dg0h_hto_counts[, keep_cells]
dg0h_gex_counts <- dg0h$`Gene Expression`[,keep_cells]

keep_cells <- colnames(v113_hto_counts[, colSums(v113_hto_counts) > 2])
v113_adt_counts <- v113_adt_counts[, keep_cells]
v113_hto_counts <- v113_hto_counts[, keep_cells]
v113_gex_counts <- v113$`Gene Expression`[,keep_cells]

keep_cells <- colnames(v015_hto_counts[, colSums(v015_hto_counts) > 2])
v015_adt_counts <- v015_adt_counts[, keep_cells]
v015_hto_counts <- v015_hto_counts[, keep_cells]
v015_gex_counts <- v015$`Gene Expression`[,keep_cells]

qc_table$hto_minimum <- c(ncol(dg0d_gex_counts), ncol(dg0h_gex_counts), ncol(v015_gex_counts), ncol(v113_gex_counts))



# Manually annotate the mitochondrial genes with an MT prefix
gtf_raw <- read.gtf("Macaca_mulatta.Mmul_10.107.filtered_with_mt.gtf", attr = "split") # This file came from running `cellranger mkgtf` on the Ensembl GTF. This was used to make our custom reference.

ensembl <- gtf_raw %>% # Load GTF, get gene names
  select(gene_id, gene_name, gene_biotype) %>%
  rename(id = gene_id,
         biotype = gene_biotype) %>%
  distinct() %>%
  # If gene_name is NA use id
  mutate(gene_name = coalesce(gene_name, id)) %>%
  # Add mitochondrial indicator
  mutate(name = case_when(
    biotype == "Mt_tRNA" | biotype == "Mt_rRNA" ~ paste0("MT-", gene_name),
    TRUE                                        ~ gene_name
  )) %>%
  melt(id.var = c('biotype', 'name'),
       variable.name = 'identifier') %>%
  select(value, biotype, name) %>%
  distinct()

rm(gtf_raw) # Remove unneeded object
gc()
ensembl %>% head(3)


# Confirm mitochondrial genes are "labelled"
ensembl %>% filter(grepl("^MT-", name)) %>% head(3)

# Update gene names in Seurat object
for (i in 1:nrow(dg0d_gex_counts)) {
  current <-  rownames(dg0d_gex_counts)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(dg0d_gex_counts)[i] <- updated
  } else {
    rownames(dg0d_gex_counts)[i] <- current
  }
}

for (i in 1:nrow(dg0h_gex_counts)) {
  current <-  rownames(dg0h_gex_counts)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(dg0h_gex_counts)[i] <- updated
  } else {
    rownames(dg0h_gex_counts)[i] <- current
  }
}

for (i in 1:nrow(v113_gex_counts)) {
  current <-  rownames(v113_gex_counts)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(v113_gex_counts)[i] <- updated
  } else {
    rownames(v113_gex_counts)[i] <- current
  }
}

for (i in 1:nrow(v015_gex_counts)) {
  current <-  rownames(v015_gex_counts)[i]
  updated <- ensembl[ensembl$value == current, ]$name
  if (length(updated) > 0) {
    rownames(v015_gex_counts)[i] <- updated
  } else {
    rownames(v015_gex_counts)[i] <- current
  }
}

rownames(v113_gex_counts)[grepl("^MT-", rownames(v113_gex_counts))] %>% head() #Confirm that the MT annotations were added


#Initialize the Seurat objects and attach the ADT and HTO counts
dg0d_s <- CreateSeuratObject(counts = dg0d_gex_counts, project = "iv_bcg_dg0d", min.cells = 3)
dg0d_s[["adt"]] <- CreateAssayObject(counts = dg0d_adt_counts)
dg0d_s[["hto"]] <- CreateAssayObject(counts = dg0d_hto_counts)

dg0h_s <- CreateSeuratObject(counts = dg0h_gex_counts, project = "iv_bcg_dg0h", min.cells = 3)
dg0h_s[["adt"]] <- CreateAssayObject(counts = dg0h_adt_counts)
dg0h_s[["hto"]] <- CreateAssayObject(counts = dg0h_hto_counts)

v113_s <- CreateSeuratObject(counts = v113_gex_counts, project = "iv_bcg_v113", min.cells = 3)
v113_s[["adt"]] <- CreateAssayObject(counts = v113_adt_counts)
v113_s[["hto"]] <- CreateAssayObject(counts = v113_hto_counts)

v015_s <- CreateSeuratObject(counts = v015_gex_counts, project = "iv_bcg_v015", min.cells = 3)
v015_s[["adt"]] <- CreateAssayObject(counts = v015_adt_counts)
v015_s[["hto"]] <- CreateAssayObject(counts = v015_hto_counts)


#Filter out the outlier cells
dg0d_s[["percent.mt"]] <- PercentageFeatureSet(dg0d_s, pattern= "^MT-") # Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
VlnPlot(dg0d_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # The Seurat Object also contains the automatically calculated metrics nCount_RNA and nFeature_RNA. These describe the total RNA count and the total number of unique features (genes) for each cell. Visualize all three of these QC metrics as a violin plot.
dg0d_s <- subset(dg0d_s, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10) # Use the violin plots to define the aberrent values. Exclude the aberrent values.
VlnPlot(dg0d_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Review the plots after filtering. 

dg0h_s[["percent.mt"]] <- PercentageFeatureSet(dg0h_s, pattern= "^MT-") # Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
VlnPlot(dg0h_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # The Seurat Object also contains the automatically calculated metrics nCount_RNA and nFeature_RNA. These describe the total RNA count and the total number of unique features (genes) for each cell. Visualize all three of these QC metrics as a violin plot.
dg0h_s <- subset(dg0h_s, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10) # Use the violin plots to define the aberrent values. Exclude the aberrent values.
VlnPlot(dg0h_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Review the plots after filtering. 

v113_s[["percent.mt"]] <- PercentageFeatureSet(v113_s, pattern= "^MT-") # Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
VlnPlot(v113_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # The Seurat Object also contains the automatically calculated metrics nCount_RNA and nFeature_RNA. These describe the total RNA count and the total number of unique features (genes) for each cell. Visualize all three of these QC metrics as a violin plot.
v113_s <- subset(v113_s, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10) # Use the violin plots to define the aberrent values. Exclude the aberrent values.
VlnPlot(v113_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Review the plots after filtering.

v015_s[["percent.mt"]] <- PercentageFeatureSet(v015_s, pattern= "^MT-") # Add a column to the Seurat object that describes the percent mitochondrial genes for each cell. All genes beginning with "MT-" are mitochondrial genes.
VlnPlot(v015_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # The Seurat Object also contains the automatically calculated metrics nCount_RNA and nFeature_RNA. These describe the total RNA count and the total number of unique features (genes) for each cell. Visualize all three of these QC metrics as a violin plot.
v015_s <- subset(v015_s, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 10) # Use the violin plots to define the aberrent values. Exclude the aberrent values.
VlnPlot(v015_s, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) # Review the plots after filtering.

dg0d_s@meta.data$nCount_RNA %>% length()
dg0h_s@meta.data$nCount_RNA %>% length()
v015_s@meta.data$nCount_RNA %>% length()
v113_s@meta.data$nCount_RNA %>% length()

qc_table$gex_filters <- c(ncol(dg0d_s), ncol(dg0h_s), ncol(v015_s), ncol(v113_s))

#Save the Seurat objects.
saveRDS(dg0d_s, file = "out/seurat_objects/dg0d_setup.rds")
saveRDS(dg0h_s, file = "out/seurat_objects/dg0h_setup.rds")
saveRDS(v113_s, file = "out/seurat_objects/v113_setup.rds")
saveRDS(v015_s, file = "out/seurat_objects/v015_setup.rds")



###################
## MULTIseqDemux ##
###################

#Read in the Seurat objects.
dg0d_s <- readRDS(file = "out/seurat_objects/dg0d_setup.rds")
dg0h_s <- readRDS(file = "out/seurat_objects/dg0h_setup.rds")
v113_s <- readRDS(file = "out/seurat_objects/v113_setup.rds")
v015_s <- readRDS(file = "out/seurat_objects/v015_setup.rds")

#Demultiplex the data using MULTIseqDemux
dg0d_s <- NormalizeData(dg0d_s, assay = "hto", normalization.method = "CLR")
RidgePlot(dg0d_s, assay = "hto", features = rownames(dg0d_s[["hto"]]), ncol = 3)
dg0d_dm <- MULTIseqDemux(dg0d_s, assay = "hto", autoThresh = TRUE)
RidgePlot(dg0d_dm, assay = "hto", features = rownames(dg0d_dm[["hto"]]), ncol = 3)
dg0d_dm_tbl <- as.data.frame(table(dg0d_dm$MULTI_ID))

dg0h_s <- NormalizeData(dg0h_s, assay = "hto", normalization.method = "CLR")
RidgePlot(dg0h_s, assay = "hto", features = rownames(dg0h_s[["hto"]]), ncol = 3)
dg0h_dm <- MULTIseqDemux(dg0h_s, assay = "hto", autoThresh = TRUE)
RidgePlot(dg0h_dm, assay = "hto", features = rownames(dg0h_dm[["hto"]]), ncol = 3)
dg0h_dm_tbl <- as.data.frame(table(dg0h_dm$MULTI_ID))

v113_s <- NormalizeData(v113_s, assay = "hto", normalization.method = "CLR")
RidgePlot(v113_s, assay = "hto", features = rownames(v113_s[["hto"]]), ncol = 3)
v113_dm <- MULTIseqDemux(v113_s, assay = "hto", autoThresh = TRUE)
RidgePlot(v113_dm, assay = "hto", features = rownames(v113_dm[["hto"]]), ncol = 3)
v113_dm_tbl <- as.data.frame(table(v113_dm$MULTI_ID))

v015_s <- NormalizeData(v015_s, assay = "hto", normalization.method = "CLR")
RidgePlot(v015_s, assay = "hto", features = rownames(v015_s[["hto"]]), ncol = 3)
v015_dm <- MULTIseqDemux(v015_s, assay = "hto", autoThresh = TRUE)
RidgePlot(v015_dm, assay = "hto", features = rownames(v015_dm[["hto"]]), ncol = 3)
v015_dm_tbl <- as.data.frame(table(v015_dm$MULTI_ID))


#Summarize the demultiplexing results
dg0d_dm_tbl <- dg0d_dm_tbl %>% pivot_wider(names_from = "Var1", values_from = "Freq")
dg0h_dm_tbl <- dg0h_dm_tbl %>% pivot_wider(names_from = "Var1", values_from = "Freq")
v113_dm_tbl <- v113_dm_tbl %>% pivot_wider(names_from = "Var1", values_from = "Freq")
v015_dm_tbl <- v015_dm_tbl %>% pivot_wider(names_from = "Var1", values_from = "Freq")

hto_tbl <- do.call(rbind, list(dg0d_dm_tbl, dg0h_dm_tbl, v113_dm_tbl, v015_dm_tbl))
hto_tbl$Sum <- rowSums(hto_tbl)
hto_tbl <- hto_tbl %>% mutate(Doublet_pct = Doublet/Sum,
                   HTO1_pct = HTO1/Sum,
                   HTO2_pct = HTO2/Sum,
                   HTO6_pct = HTO6/Sum,
                   HTO7_pct = HTO7/Sum,
                   HTO8_pct = HTO8/Sum,
                   Negative_pct = Negative/Sum)
rownames(hto_tbl) <- c("dg0d", "dg0h", "v113", "v015")
write.csv(hto_tbl, file = "out/hto_assignments_by_PTID_msdemux.csv")


#Visualize the RNA counts of singlets, doublets, and negatives for each PTID
Idents(dg0d_dm) <- "MULTI_ID"
VlnPlot(dg0d_dm, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

Idents(dg0h_dm) <- "MULTI_ID"
VlnPlot(dg0h_dm, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

Idents(v113_dm) <- "MULTI_ID"
VlnPlot(v113_dm, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

Idents(v015_dm) <- "MULTI_ID"
VlnPlot(v015_dm, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


#Cluster the data within each PTID based on HTO expression
#dg0d
Idents(dg0d_dm) <- "MULTI_ID"
DefaultAssay(dg0d_dm) <- "hto"
dg0d_dm <- ScaleData(dg0d_dm, features = rownames(dg0d_dm),
                     verbose = FALSE)
dg0d_dm <- RunPCA(dg0d_dm, features = rownames(dg0d_dm), approx = FALSE)
dg0d_dm <- RunTSNE(dg0d_dm, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(dg0d_dm)

#dg0h
Idents(dg0h_dm) <- "MULTI_ID"
DefaultAssay(dg0h_dm) <- "hto"
dg0h_dm <- ScaleData(dg0h_dm, features = rownames(dg0h_dm),
                     verbose = FALSE)
dg0h_dm <- RunPCA(dg0h_dm, features = rownames(dg0h_dm), approx = FALSE)
dg0h_dm <- RunTSNE(dg0h_dm, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(dg0h_dm)

#v113
Idents(v113_dm) <- "MULTI_ID"
DefaultAssay(v113_dm) <- "hto"
v113_dm <- ScaleData(v113_dm, features = rownames(v113_dm),
                     verbose = FALSE)
v113_dm <- RunPCA(v113_dm, features = rownames(v113_dm), approx = FALSE)
v113_dm <- RunTSNE(v113_dm, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(v113_dm)

#v015
Idents(v015_dm) <- "MULTI_ID"
DefaultAssay(v015_dm) <- "hto"
v015_dm <- ScaleData(v015_dm, features = rownames(v015_dm),
                     verbose = FALSE)
v015_dm <- RunPCA(v015_dm, features = rownames(v015_dm), approx = FALSE)
v015_dm <- RunTSNE(v015_dm, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(v015_dm)

#v015
Idents(v015_dm) <- "MULTI_ID"
DefaultAssay(v015_dm) <- "hto"
v015_dm <- ScaleData(v015_dm, features = rownames(v015_dm),
                     verbose = FALSE)
v015_dm <- RunPCA(v015_dm, features = rownames(v015_dm), approx = FALSE)
v015_dm <- RunTSNE(v015_dm, dims = 1:5, perplexity = 100, check_duplicates = FALSE)
DimPlot(v015_dm)



####################
## HTO Annotation ##
####################

#Annotate PTID, tissue, and timepoint for each HTO
#dg0d
i1 <- grepl("HTO1", dg0d_dm$MULTI_ID)
i2 <- grepl("HTO2", dg0d_dm$MULTI_ID)
i3 <- grepl("HTO6", dg0d_dm$MULTI_ID)
i4 <- grepl("HTO7", dg0d_dm$MULTI_ID) 
i5 <- grepl("HTO8", dg0d_dm$MULTI_ID) 
i6 <- grepl("Doublet", dg0d_dm$MULTI_ID) 
i7 <- grepl("Negative", dg0d_dm$MULTI_ID)

dg0d_dm@meta.data$PTID[i1] <- "DG0D" #PTID
dg0d_dm@meta.data$PTID[i2] <- "DG0D"
dg0d_dm@meta.data$PTID[i3] <- "DG0D"
dg0d_dm@meta.data$PTID[i4] <- "DG0D"
dg0d_dm@meta.data$PTID[i5] <- "DG0D"
dg0d_dm@meta.data$PTID[i6] <- "Doublet"
dg0d_dm@meta.data$PTID[i7] <- "Negative"

dg0d_dm@meta.data$Tissue[i1] <- "PBMC" #Tissue
dg0d_dm@meta.data$Tissue[i2] <- "PBMC"
dg0d_dm@meta.data$Tissue[i3] <- "PBMC"
dg0d_dm@meta.data$Tissue[i4] <- "BAL"
dg0d_dm@meta.data$Tissue[i5] <- "BAL"
dg0d_dm@meta.data$Tissue[i6] <- "Doublet"
dg0d_dm@meta.data$Tissue[i7] <- "Negative"

dg0d_dm@meta.data$Timepoint[i1] <- "Week_0" #Timepoint
dg0d_dm@meta.data$Timepoint[i2] <- "Week_4"
dg0d_dm@meta.data$Timepoint[i3] <- "Week_8"
dg0d_dm@meta.data$Timepoint[i4] <- "Week_4"
dg0d_dm@meta.data$Timepoint[i5] <- "Week_8"
dg0d_dm@meta.data$Timepoint[i6] <- "Doublet"
dg0d_dm@meta.data$Timepoint[i7] <- "Negative"


#dg0h
i1 <- grepl("HTO1", dg0h_dm$MULTI_ID)
i2 <- grepl("HTO2", dg0h_dm$MULTI_ID)
i3 <- grepl("HTO6", dg0h_dm$MULTI_ID)
i4 <- grepl("HTO7", dg0h_dm$MULTI_ID) 
i5 <- grepl("HTO8", dg0h_dm$MULTI_ID) 
i6 <- grepl("Doublet", dg0h_dm$MULTI_ID) 
i7 <- grepl("Negative", dg0h_dm$MULTI_ID)

dg0h_dm@meta.data$PTID[i1] <- "DG0H" #PTID
dg0h_dm@meta.data$PTID[i2] <- "DG0H"
dg0h_dm@meta.data$PTID[i3] <- "DG0H"
dg0h_dm@meta.data$PTID[i4] <- "DG0H"
dg0h_dm@meta.data$PTID[i5] <- "DG0H"
dg0h_dm@meta.data$PTID[i6] <- "Doublet"
dg0h_dm@meta.data$PTID[i7] <- "Negative"

dg0h_dm@meta.data$Tissue[i1] <- "PBMC" #Tissue
dg0h_dm@meta.data$Tissue[i2] <- "PBMC"
dg0h_dm@meta.data$Tissue[i3] <- "PBMC"
dg0h_dm@meta.data$Tissue[i4] <- "BAL"
dg0h_dm@meta.data$Tissue[i5] <- "BAL"
dg0h_dm@meta.data$Tissue[i6] <- "Doublet"
dg0h_dm@meta.data$Tissue[i7] <- "Negative"

dg0h_dm@meta.data$Timepoint[i1] <- "Week_0" #Timepoint
dg0h_dm@meta.data$Timepoint[i2] <- "Week_4"
dg0h_dm@meta.data$Timepoint[i3] <- "Week_8"
dg0h_dm@meta.data$Timepoint[i4] <- "Week_4"
dg0h_dm@meta.data$Timepoint[i5] <- "Week_8"
dg0h_dm@meta.data$Timepoint[i6] <- "Doublet"
dg0h_dm@meta.data$Timepoint[i7] <- "Negative"


#v113
i1 <- grepl("HTO1", v113_dm$MULTI_ID)
i2 <- grepl("HTO2", v113_dm$MULTI_ID)
i3 <- grepl("HTO6", v113_dm$MULTI_ID)
i4 <- grepl("HTO7", v113_dm$MULTI_ID) 
i5 <- grepl("HTO8", v113_dm$MULTI_ID) 
i6 <- grepl("Doublet", v113_dm$MULTI_ID) 
i7 <- grepl("Negative", v113_dm$MULTI_ID)

v113_dm@meta.data$PTID[i1] <- "V113" #PTID
v113_dm@meta.data$PTID[i2] <- "V113"
v113_dm@meta.data$PTID[i3] <- "V113"
v113_dm@meta.data$PTID[i4] <- "V113"
v113_dm@meta.data$PTID[i5] <- "V113"
v113_dm@meta.data$PTID[i6] <- "Doublet"
v113_dm@meta.data$PTID[i7] <- "Negative"

v113_dm@meta.data$Tissue[i1] <- "PBMC" #Tissue
v113_dm@meta.data$Tissue[i2] <- "PBMC"
v113_dm@meta.data$Tissue[i3] <- "PBMC"
v113_dm@meta.data$Tissue[i4] <- "BAL"
v113_dm@meta.data$Tissue[i5] <- "BAL"
v113_dm@meta.data$Tissue[i6] <- "Doublet"
v113_dm@meta.data$Tissue[i7] <- "Negative"

v113_dm@meta.data$Timepoint[i1] <- "Week_0" #Timepoint
v113_dm@meta.data$Timepoint[i2] <- "Week_4"
v113_dm@meta.data$Timepoint[i3] <- "Week_8"
v113_dm@meta.data$Timepoint[i4] <- "Week_4"
v113_dm@meta.data$Timepoint[i5] <- "Week_8"
v113_dm@meta.data$Timepoint[i6] <- "Doublet"
v113_dm@meta.data$Timepoint[i7] <- "Negative"


#v015
i1 <- grepl("HTO1", v015_dm$MULTI_ID)
i2 <- grepl("HTO2", v015_dm$MULTI_ID)
i3 <- grepl("HTO6", v015_dm$MULTI_ID)
i4 <- grepl("HTO7", v015_dm$MULTI_ID) 
i5 <- grepl("HTO8", v015_dm$MULTI_ID) 
i6 <- grepl("Doublet", v015_dm$MULTI_ID) 
i7 <- grepl("Negative", v015_dm$MULTI_ID)

v015_dm@meta.data$PTID[i1] <- "V015" #PTID
v015_dm@meta.data$PTID[i2] <- "V015"
v015_dm@meta.data$PTID[i3] <- "V015"
v015_dm@meta.data$PTID[i4] <- "V015"
v015_dm@meta.data$PTID[i5] <- "V015"
v015_dm@meta.data$PTID[i6] <- "Doublet"
v015_dm@meta.data$PTID[i7] <- "Negative"

v015_dm@meta.data$Tissue[i1] <- "PBMC" #Tissue
v015_dm@meta.data$Tissue[i2] <- "PBMC"
v015_dm@meta.data$Tissue[i3] <- "PBMC"
v015_dm@meta.data$Tissue[i4] <- "BAL"
v015_dm@meta.data$Tissue[i5] <- "BAL"
v015_dm@meta.data$Tissue[i6] <- "Doublet"
v015_dm@meta.data$Tissue[i7] <- "Negative"

v015_dm@meta.data$Timepoint[i1] <- "Week_0" #Timepoint
v015_dm@meta.data$Timepoint[i2] <- "Week_4"
v015_dm@meta.data$Timepoint[i3] <- "Week_8"
v015_dm@meta.data$Timepoint[i4] <- "Week_4"
v015_dm@meta.data$Timepoint[i5] <- "Week_8"
v015_dm@meta.data$Timepoint[i6] <- "Doublet"
v015_dm@meta.data$Timepoint[i7] <- "Negative"


#Exclude the doublets and negative cells
dg0d_dm_subset <- subset(dg0d_dm, idents = "Negative", invert = TRUE) %>% subset(idents = "Doublet", invert = TRUE)
dg0h_dm_subset <- subset(dg0h_dm, idents = "Negative", invert = TRUE) %>% subset(idents = "Doublet", invert = TRUE)
v015_dm_subset <- subset(v015_dm, idents = "Negative", invert = TRUE) %>% subset(idents = "Doublet", invert = TRUE)
v113_dm_subset <- subset(v113_dm, idents = "Negative", invert = TRUE) %>% subset(idents = "Doublet", invert = TRUE)


qc_table$hto_singlets <- c(ncol(dg0d_dm_subset), ncol(dg0h_dm_subset), ncol(v015_dm_subset), ncol(v113_dm_subset))

#Save the Seurat objects
saveRDS(dg0d_dm_subset, file = "out/seurat_objects/dg0d_demux.rds")
saveRDS(dg0h_dm_subset, file = "out/seurat_objects/dg0h_demux.rds")
saveRDS(v113_dm_subset, file = "out/seurat_objects/v113_demux.rds")
saveRDS(v015_dm_subset, file = "out/seurat_objects/v015_demux.rds")

write.csv(qc_table, 'qc_table.csv')
