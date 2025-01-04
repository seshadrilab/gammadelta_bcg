#Load dplyr, Seurat, and patchwork
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(devtools)
library(scRepertoire)
library(tibble)

#Read in the filtered_contig_annotations files and the T cell seurat objects
dg0d_t <- read.csv(file = 'Tube1_vdj_v3.0.2_DG0D/filtered_contig_annotations.csv')
dg0h_t <- read.csv(file = 'Tube1_vdj_v3.0.2_DG0H/filtered_contig_annotations.csv')
v113_t <- read.csv(file = 'Tube1_vdj_v3.0.2_V113/filtered_contig_annotations.csv')
v015_t <- read.csv(file = 'Tube1_vdj_v3.0.2_V015/filtered_contig_annotations.csv')

dg0d_g <- readRDS(file = 'out/seurat_objects/dg0d_T.rds')
dg0h_g <- readRDS(file = 'out/seurat_objects/dg0h_T.rds')
v113_g <- readRDS(file = 'out/seurat_objects/v113_T.rds')
v015_g <- readRDS(file = 'out/seurat_objects/v015_T.rds')


#Filter out all TCR barcodes not present in the _analysis.rds files. This maintains only the barcodes that passed GEX QC.
filter_tcr <- function(gex_csv, tcr_csv) { 
  gex_bc <- colnames(gex_csv) %>% as.data.frame()
  colnames(gex_bc) <- "barcode"
  gex_bc$gex_pass_qc <- "true"
  tcr_qcd <- merge(x = gex_bc, y = tcr_csv, by = "barcode")
  return(tcr_qcd)
}

dg0d_qcd <- filter_tcr(dg0d_g, dg0d_t)
dg0h_qcd <- filter_tcr(dg0h_g, dg0h_t)
v113_qcd <- filter_tcr(v113_g, v113_t)
v015_qcd <- filter_tcr(v015_g, v015_t)


# Maintain TCRs satisfying all of these criteria:
# 1. A productive rearrangement
# 2. At least three UMIs
# 3. At least 20 reads


#Maintain productive TCRs
dg0d_qcd <- dg0d_qcd %>% filter(productive == "True")
dg0h_qcd <- dg0h_qcd %>% filter(productive == "True")
v113_qcd <- v113_qcd %>% filter(productive == "True")
v015_qcd <- v015_qcd %>% filter(productive == "True")

#Measure the proportion of TCRs with low UMI counts
low_umi <- function(df) {
  all <- df$barcode %>% length()
  one <- df %>% filter(umis == 1) %>% nrow()
  two <- df %>% filter(umis == 2) %>% nrow()
  three <- df %>% filter(umis == 3) %>% nrow()
  high <- df %>% filter(umis >3) %>% nrow()
  result <- c("UMI = 1" = one/all, "UMI = 2" = two/all, "UMI = 3" = three/all, "UMI >3" = high/all)
  return(result)
}

low_umi(dg0d_qcd)
low_umi(dg0h_qcd)
low_umi(v113_qcd)
low_umi(v015_qcd)


#Maintain TCRs with UMIs>2 and reads>20
dg0d_qcd <- dg0d_qcd %>% filter(umis >2, reads >20)
dg0h_qcd <- dg0h_qcd %>% filter(umis >2, reads >20)
v113_qcd <- v113_qcd %>% filter(umis >2, reads >20)
v015_qcd <- v015_qcd %>% filter(umis >2, reads >20)


#Use scRepertoire to define a single TRD and TRG per chain
contig_list <- list(dg0d_qcd, dg0h_qcd, v113_qcd, v015_qcd)
combined <- combineTCR(contig_list, samples = c("dg0d", "dg0h", "v113", "v015"), ID = NULL, cells = "T-GD", removeNA = TRUE, removeMulti = FALSE, filterMulti = TRUE) #Keeps only T cells with paired-chain data
combined_sch <- combineTCR(contig_list, samples = c("dg0d", "dg0h", "v113", "v015"), ID = NULL, cells = "T-GD", removeNA = FALSE, removeMulti = FALSE, filterMulti = TRUE) #Keeps T cells with at least one chain


#Export the combined TCR dfs
dg0d_pch <- combined$dg0d
dg0d_sch <- combined_sch$dg0d
write.csv(dg0d_pch, file= "out/tcr/dg0d_tcr_paired.csv")
write.csv(dg0d_sch, file = "out/tcr/dg0d_tcr_sch.csv")

dg0h_pch <- combined$dg0h
dg0h_sch <- combined_sch$dg0h
write.csv(dg0h_pch, file= "out/tcr/dg0h_tcr_paired.csv")
write.csv(dg0h_sch, file = "out/tcr/dg0h_tcr_sch.csv")

v113_pch <- combined$v113
v113_sch <- combined_sch$v113
write.csv(v113_pch, file= "out/tcr/v113_tcr_paired.csv")
write.csv(v113_sch, file = "out/tcr/v113_tcr_sch.csv")

v015_pch <- combined$v015
v015_sch <- combined_sch$v015
write.csv(v015_pch, file= "out/tcr/v015_tcr_paired.csv")
write.csv(v015_sch, file = "out/tcr/v015_tcr_sch.csv")


#List the cell barcodes with at least one chain passing QC and export.
barcodes <- c(dg0d_sch$barcode)
barcodes2 <- gsub("dg0d_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0d_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0d_sch <- subset(dg0d_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0d_sch <- dg0d_sch$barcode %>% as.data.frame()
names(dg0d_sch) <- c("barcode")
write.csv(dg0d_sch, "out/tcr/dg0d_gd_bc.csv")

barcodes <- c(dg0h_sch$barcode)
barcodes2 <- gsub("dg0h_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0h_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0h_sch <- subset(dg0h_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0h_sch <- dg0h_sch$barcode %>% as.data.frame()
names(dg0h_sch) <- c("barcode")
write.csv(dg0h_sch, "out/tcr/dg0h_gd_bc.csv")

barcodes <- c(v113_sch$barcode)
barcodes2 <- gsub("v113_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v113_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v113_sch <- subset(v113_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v113_sch <- v113_sch$barcode %>% as.data.frame()
names(v113_sch) <- c("barcode")
write.csv(v113_sch, "out/tcr/v113_gd_bc.csv")

barcodes <- c(v015_sch$barcode)
barcodes2 <- gsub("v015_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v015_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v015_sch <- subset(v015_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v015_sch <- v015_sch$barcode %>% as.data.frame()
names(v015_sch) <- c("barcode")
write.csv(v015_sch, "out/tcr/v015_gd_bc.csv")




###########
## NOTES ##
###########

## Generate new lists to allow comparisons between samples within a single PTID ##
#DG0D
dg0d_t <- combined$dg0d

barcodes <- c(dg0d_t$barcode) #Store the barcodes
barcodes2 <- gsub("dg0d_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0d_t$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0d_t <- subset(dg0d_t, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename

dg0d_g
dg0d_anno <- FetchData(object = dg0d_g, vars = c("Tissue", "Timepoint"))
dg0d_anno <- tibble::rownames_to_column(dg0d_anno, "barcode")

dg0d_t <- merge(dg0d_t, dg0d_anno, by = "barcode")
dg0d_t$sample <- paste(dg0d_t$Tissue, dg0d_t$Timepoint, sep = "_")
dg0d_split <- dg0d_t %>% split(dg0d_t, f = dg0d_t$sample)


#DG0H
dg0h_t <- combined$dg0h

barcodes <- c(dg0h_t$barcode) #Store the barcodes
barcodes2 <- gsub("dg0h_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0h_t$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0h_t <- subset(dg0h_t, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename

dg0h_anno <- FetchData(object = dg0h_g, vars = c("Tissue", "Timepoint"))
dg0h_anno <- tibble::rownames_to_column(dg0h_anno, "barcode")

dg0h_t <- merge(dg0h_t, dg0h_anno, by = "barcode")
dg0h_t$sample <- paste(dg0h_t$Tissue, dg0h_t$Timepoint, sep = "_")
dg0h_split <- dg0h_t %>% split(dg0h_t, f = dg0h_t$sample)


#V113
v113_t <- combined$v113

barcodes <- c(v113_t$barcode) #Store the barcodes
barcodes2 <- gsub("v113_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v113_t$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v113_t <- subset(v113_t, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename

v113_anno <- FetchData(object = v113_g, vars = c("Tissue", "Timepoint"))
v113_anno <- tibble::rownames_to_column(v113_anno, "barcode")

v113_t <- merge(v113_t, v113_anno, by = "barcode")
v113_t$sample <- paste(v113_t$Tissue, v113_t$Timepoint, sep = "_")
v113_split <- v113_t %>% split(v113_t, f = v113_t$sample)


#V015
v015_t <- combined$v015

barcodes <- c(v015_t$barcode) #Store the barcodes
barcodes2 <- gsub("v015_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v015_t$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v015_t <- subset(v015_t, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename

v015_anno <- FetchData(object = v015_g, vars = c("Tissue", "Timepoint"))
v015_anno <- tibble::rownames_to_column(v015_anno, "barcode")

v015_t <- merge(v015_t, v015_anno, by = "barcode")
v015_t$sample <- paste(v015_t$Tissue, v015_t$Timepoint, sep = "_")
v015_split <- v015_t %>% split(v015_t, f = v015_t$sample)


## Analyze the split TCR data using scRepertoire ##

#DG0D
quantContig(dg0d_split, cloneCall="nt", scale = T)
abundanceContig(dg0d_split, cloneCall = "nt", scale = T)
lengthContig(dg0d_split, cloneCall = "nt", chain = "both", scale = T)
compareClonotypes(dg0d_split, 
                  chain = "both", 
                  cloneCall = "aa", 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  graph = "alluvial")
clonalProportion(dg0d_split, cloneCall = "nt", 
                  split = c(1, 10, 1e+05))

#DG0H
quantContig(dg0h_split, cloneCall="nt", scale = T)
abundanceContig(dg0h_split, cloneCall = "nt", scale = T)
lengthContig(dg0h_split, cloneCall = "nt", chain = "both", scale = T)
compareClonotypes(dg0h_split,
                  numbers = 20,
                  chain = "both", 
                  cloneCall = "aa", 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  graph = "alluvial")
clonalProportion(dg0h_split, cloneCall = "nt", 
                 split = c(1, 10, 1e+05))

#V113
quantContig(v113_split, cloneCall="nt", scale = T)
abundanceContig(v113_split, cloneCall = "nt", scale = T)
lengthContig(v113_split, cloneCall = "nt", chain = "both", scale = T)
compareClonotypes(v113_split,
                  numbers = 20,
                  chain = "both", 
                  cloneCall = "aa", 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  graph = "alluvial")
clonalProportion(v113_split, cloneCall = "nt", 
                 split = c(1, 10, 1e+05))

#V015
quantContig(v015_split, cloneCall="nt", scale = T)
abundanceContig(v015_split, cloneCall = "nt", scale = T)
lengthContig(v015_split, cloneCall = "nt", chain = "both", scale = T)
compareClonotypes(v015_split,
                  numbers = 20,
                  chain = "both", 
                  cloneCall = "aa", 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  graph = "alluvial")
clonalProportion(v015_split, cloneCall = "nt", 
                 split = c(1, 10, 1e+05))


#Measure the proportion of barcodes with >2 contigs per chain
dg0d_rv <- dg0d_qcd %>% filter(!grepl("contig_1", contig_id)) %>% filter(!grepl("contig_2", contig_id))
dg0d_rv_bc <- dg0d_rv$barcode %>% as.data.frame()
colnames(dg0d_rv_bc) <- "barcode"
check1 <- merge(x = dg0d_qcd, y = dg0d_rv_bc, by = "barcode")

v113_rv <- v113_qcd %>% filter(!grepl("contig_1", contig_id)) %>% filter(!grepl("contig_2", contig_id))
v113_rv_bc <- v113_rv$barcode %>% as.data.frame()
colnames(v113_rv_bc) <- "barcode"
check2 <- merge(x = v113_qcd, y = v113_rv_bc, by = "barcode")
check2 <- check2 %>% distinct()
