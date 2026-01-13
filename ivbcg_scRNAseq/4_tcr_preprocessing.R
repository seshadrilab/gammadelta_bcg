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

# Read in the filtered_contig_annotations files and the T cell seurat objects
dg0d_t <- read.csv(file = 'Tube1_vdj_v3.0.2_DG0D/filtered_contig_annotations.csv')
dg0h_t <- read.csv(file = 'Tube1_vdj_v3.0.2_DG0H/filtered_contig_annotations.csv')
v113_t <- read.csv(file = 'Tube1_vdj_v3.0.2_V113/filtered_contig_annotations.csv')
v015_t <- read.csv(file = 'Tube1_vdj_v3.0.2_V015/filtered_contig_annotations.csv')

dg0d_g <- readRDS(file = 'out/seurat_objects/dg0d_T.rds')
dg0h_g <- readRDS(file = 'out/seurat_objects/dg0h_T.rds')
v113_g <- readRDS(file = 'out/seurat_objects/v113_T.rds')
v015_g <- readRDS(file = 'out/seurat_objects/v015_T.rds')


# Filter out all TCR barcodes not present in the RDS files. This maintains only the barcodes that passed GEX QC.
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


# Use scRepertoire to define a single TRD and TRG per chain
# Generate lists of (1) paired-chain TCRs and (2) barcodes from which only one productive chain was measured
contig_list <- list(dg0d_qcd, dg0h_qcd, v113_qcd, v015_qcd)
combined <- combineTCR(contig_list, 
                       samples = c("dg0d", "dg0h", "v113", "v015"), 
                       ID = NULL, 
                       cells = "T-GD", 
                       removeNA = TRUE, 
                       removeMulti = FALSE, 
                       filterMulti = TRUE) #Keeps only T cells with paired-chain data
combined_sch <- combineTCR(contig_list, 
                           samples = c("dg0d", "dg0h", "v113", "v015"), 
                           ID = NULL, 
                           cells = "T-GD", 
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = TRUE) #Keeps T cells with at least one chain


#Export the combined TCR dataframes
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
