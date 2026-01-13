library(ggplot2)
library(data.table)
library(tidyverse)
library(Seurat)

# Read in the necessary files
pbmc <- readRDS(file = "GEX/Out/pbmc_tcells_2024.rds")
fmux <- file("Infant_freemuxlet/unzipped/freemuxlet.clust1.samples", "rb")

# Convert the binary file from freemuxlet into a data frame
fmux_str <- readChar(fmux, 10^7, useBytes = FALSE)
fmux_df <- fread(fmux_str)

# Attach the donor assignments from freemuxlet to the seurat object
fmux_df <- fmux_df %>%
  filter(DROPLET.TYPE == "SNG") %>%
  subset(select = c(BARCODE, SNG.BEST.GUESS))
fmux_df <- fmux_df %>% rename("barcode" = "BARCODE", "ptid" = "SNG.BEST.GUESS")
fmux_df$ptid <- sub("^", "ptid", fmux_df$ptid)
rownames(fmux_df) <- fmux_df$barcode
pbmc <- AddMetaData(pbmc, metadata = fmux_df)

# Visualize the contribution of each donor to the seurat clusters
Idents(pbmc) <- "ptid"
DimPlot(pbmc, reduction = "umap")

ggsave("ptid_clust.pdf")

# Export the RDS file with PTID annotations included
saveRDS(pbmc, file = "GEX/Out/pbmc_tcells_2024.rds")


