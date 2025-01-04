library(ggplot2)
library(data.table)
library(tidyverse)
library(Seurat)

# Read in the necessary files
pbmc <- readRDS(file = 'GEX/out/pbmc_tcells_2024.rds')
fmux <- file("Infant_freemuxlet/unzipped/freemuxlet.clust1.samples", "rb")

# Convert the binary file into a data frame
fmux_str <- readChar(fmux, 10^7, useBytes = FALSE)
fmux_df <- fread(fmux_str)

# Attach the PTID assignments to the seurat object
fmux_df <- fmux_df %>% filter(DROPLET.TYPE == "SNG") %>% subset(select = c(BARCODE, SNG.BEST.GUESS))
fmux_df <- fmux_df %>% rename("barcode" = "BARCODE", "ptid" = "SNG.BEST.GUESS")
fmux_df$ptid <- sub("^", "ptid", fmux_df$ptid)
rownames(fmux_df) <- fmux_df$barcode
pbmc <- AddMetaData(pbmc, metadata = fmux_df)

#Visualize the contribution of each PTID to the seurat clusters
Idents(pbmc) <- "ptid"
DimPlot(pbmc, reduction = "umap")

ggsave("ptid_clust.pdf")

# Visualize the fraction of cells in C5 coming from each PTID
md <- pbmc@meta.data
md <- md %>% subset(select = c(seurat_clusters, ptid)) %>% table() %>% as.data.frame() %>% filter(seurat_clusters == "5")
md <- md %>% mutate(pct = Freq/sum(md$Freq)*100)

pie(md$pct, labels = c("PTID_0", "PTID_1", "PTID_2", "PTID_3", "PTID_4", "PTID_5"))

ggplot(md, aes(x = "", y = Freq, fill = ptid))+
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  theme_void()

# Export the RDS file with PTID annotations included
saveRDS(pbmc, file = "GEX/Out/pbmc_tcells_2024.rds")


