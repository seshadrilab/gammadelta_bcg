library(dplyr)
library(scRepertoire)
library(Seurat)
library(circlize)
library(scales)
library(tidyverse)
library(ggpubr)

# Read in the TCR and RDS files. Files labeled _p only include cell barcodes with paired-chain data. Files labeled _sch include cell barcodes with at least one productive chain.
dg0d_p <- read.csv(file= "out/tcr/dg0d_tcr_paired.csv")
dg0d_sch <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_p <- read.csv(file= "out/tcr/dg0h_tcr_paired.csv")
dg0h_sch <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_p <- read.csv(file= "out/tcr/v113_tcr_paired.csv")
v113_sch <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_p <- read.csv(file= "out/tcr/v015_tcr_paired.csv")
v015_sch <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

dg0d_g <- readRDS(file = 'out/seurat_objects/dg0d_T.rds')
dg0h_g <- readRDS(file = 'out/seurat_objects/dg0h_T.rds')
v113_g <- readRDS(file = 'out/seurat_objects/v113_T.rds')
v015_g <- readRDS(file = 'out/seurat_objects/v015_T.rds')


# Filter the CSVs so they only include Vd1/3/4 data
dg0d_p$TRD = substr(dg0d_p$TCR2, 1, 5)
dg0d_p$TRG = substr(dg0d_p$TCR1, 1, 5)
dg0d_sch$TRD = substr(dg0d_sch$TCR2, 1, 5)
dg0d_sch$TRG = substr(dg0d_sch$TCR1, 1, 5)
dg0d_p <- dg0d_p %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")
dg0d_sch <- dg0d_sch %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")

dg0h_p$TRD = substr(dg0h_p$TCR2, 1, 5)
dg0h_p$TRG = substr(dg0h_p$TCR1, 1, 5)
dg0h_sch$TRD = substr(dg0h_sch$TCR2, 1, 5)
dg0h_sch$TRG = substr(dg0h_sch$TCR1, 1, 5)
dg0h_p <- dg0h_p %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")
dg0h_sch <- dg0h_sch %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")

v113_p$TRD = substr(v113_p$TCR2, 1, 5)
v113_p$TRG = substr(v113_p$TCR1, 1, 5)
v113_sch$TRD = substr(v113_sch$TCR2, 1, 5)
v113_sch$TRG = substr(v113_sch$TCR1, 1, 5)
v113_p <- v113_p %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")
v113_sch <- v113_sch %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")

v015_p$TRD = substr(v015_p$TCR2, 1, 5)
v015_p$TRG = substr(v015_p$TCR1, 1, 5)
v015_sch$TRD = substr(v015_sch$TCR2, 1, 5)
v015_sch$TRG = substr(v015_sch$TCR1, 1, 5)
v015_p <- v015_p %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")
v015_sch <- v015_sch %>% filter(TRD == "TRDV1" | TRD == "TRDV3" | TRD == "TRDV4")

## Create the combined lists for analysis in scRepertoire ##
#DG0D
barcodes <- c(dg0d_sch$barcode) #Store the barcodes
barcodes2 <- gsub("dg0d_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0d_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0d_sch <- subset(dg0d_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0d_anno <- FetchData(object = dg0d_g, vars = c("Tissue", "Timepoint"))
dg0d_anno <- tibble::rownames_to_column(dg0d_anno, "barcode")
dg0d_sch <- merge(dg0d_sch, dg0d_anno, by = "barcode")
dg0d_sch$sample <- paste(dg0d_sch$Tissue, dg0d_sch$Timepoint, sep = "_")
dg0d_split_sch <- dg0d_sch %>% split(dg0d_sch, f = dg0d_sch$sample)

barcodes <- c(dg0d_p$barcode) #Store the barcodes
barcodes2 <- gsub("dg0d_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0d_p$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0d_p <- subset(dg0d_p, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0d_anno <- FetchData(object = dg0d_g, vars = c("Tissue", "Timepoint"))
dg0d_anno <- tibble::rownames_to_column(dg0d_anno, "barcode")
dg0d_p <- merge(dg0d_p, dg0d_anno, by = "barcode")
dg0d_p$sample <- paste(dg0d_p$Tissue, dg0d_p$Timepoint, sep = "_")
dg0d_split_p <- dg0d_p %>% split(dg0d_p, f = dg0d_p$sample)


#DG0H
barcodes <- c(dg0h_sch$barcode) #Store the barcodes
barcodes2 <- gsub("dg0h_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0h_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0h_sch <- subset(dg0h_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0h_anno <- FetchData(object = dg0h_g, vars = c("Tissue", "Timepoint"))
dg0h_anno <- tibble::rownames_to_column(dg0h_anno, "barcode")
dg0h_sch <- merge(dg0h_sch, dg0h_anno, by = "barcode")
dg0h_sch$sample <- paste(dg0h_sch$Tissue, dg0h_sch$Timepoint, sep = "_")
dg0h_split_sch <- dg0h_sch %>% split(dg0h_sch, f = dg0h_sch$sample)

barcodes <- c(dg0h_p$barcode) #Store the barcodes
barcodes2 <- gsub("dg0h_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
dg0h_p$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
dg0h_p <- subset(dg0h_p, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
dg0h_anno <- FetchData(object = dg0h_g, vars = c("Tissue", "Timepoint"))
dg0h_anno <- tibble::rownames_to_column(dg0h_anno, "barcode")
dg0h_p <- merge(dg0h_p, dg0h_anno, by = "barcode")
dg0h_p$sample <- paste(dg0h_p$Tissue, dg0h_p$Timepoint, sep = "_")
dg0h_split_p <- dg0h_p %>% split(dg0h_p, f = dg0h_p$sample)


#V113
barcodes <- c(v113_sch$barcode) #Store the barcodes
barcodes2 <- gsub("v113_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v113_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v113_sch <- subset(v113_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v113_anno <- FetchData(object = v113_g, vars = c("Tissue", "Timepoint"))
v113_anno <- tibble::rownames_to_column(v113_anno, "barcode")
v113_sch <- merge(v113_sch, v113_anno, by = "barcode")
v113_sch$sample <- paste(v113_sch$Tissue, v113_sch$Timepoint, sep = "_")
v113_split_sch <- v113_sch %>% split(v113_sch, f = v113_sch$sample)

barcodes <- c(v113_p$barcode) #Store the barcodes
barcodes2 <- gsub("v113_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v113_p$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v113_p <- subset(v113_p, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v113_anno <- FetchData(object = v113_g, vars = c("Tissue", "Timepoint"))
v113_anno <- tibble::rownames_to_column(v113_anno, "barcode")
v113_p <- merge(v113_p, v113_anno, by = "barcode")
v113_p$sample <- paste(v113_p$Tissue, v113_p$Timepoint, sep = "_")
v113_split_p <- v113_p %>% split(v113_p, f = v113_p$sample)


#V015
barcodes <- c(v015_sch$barcode) #Store the barcodes
barcodes2 <- gsub("v015_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v015_sch$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v015_sch <- subset(v015_sch, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v015_anno <- FetchData(object = v015_g, vars = c("Tissue", "Timepoint"))
v015_anno <- tibble::rownames_to_column(v015_anno, "barcode")
v015_sch <- merge(v015_sch, v015_anno, by = "barcode")
v015_sch$sample <- paste(v015_sch$Tissue, v015_sch$Timepoint, sep = "_")
v015_split_sch <- v015_sch %>% split(v015_sch, f = v015_sch$sample)

barcodes <- c(v015_p$barcode) #Store the barcodes
barcodes2 <- gsub("v015_", "", barcodes, fixed = TRUE) #Remove the prefix from the stored barcodes
v015_p$barcode2 <- c(barcodes2) #Replace the barcodes column with the cleaned-up barcodes
v015_p <- subset(v015_p, select=-c(barcode)) %>% rename(barcode = barcode2) #Drop the replicate column and rename
v015_anno <- FetchData(object = v015_g, vars = c("Tissue", "Timepoint"))
v015_anno <- tibble::rownames_to_column(v015_anno, "barcode")
v015_p <- merge(v015_p, v015_anno, by = "barcode")
v015_p$sample <- paste(v015_p$Tissue, v015_p$Timepoint, sep = "_")
v015_split_p <- v015_p %>% split(v015_p, f = v015_p$sample)



######################################
## Analyze the data in scRepertoire ##
######################################

## DG0D ##
dg0d_split_sch_pbmc <- dg0d_split_sch[c(3, 4, 5)]
dg0d_split_sch <- dg0d_split_sch[c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8", "BAL_Week_4", "BAL_Week_8")]

quantContig(dg0d_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRD", scale = T)
quantContig(dg0d_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRG", scale = T)
quantContig(dg0d_split_sch_pbmc, cloneCall = "gene+nt", chain = "both", scale = T)

abundanceContig(dg0d_split_sch_pbmc, cloneCall = "nt", scale = T)

lengthContig(dg0d_split_sch_pbmc, cloneCall = "nt", chain = "TRD", scale = T)
lengthContig(dg0d_split_sch_pbmc, cloneCall = "nt", chain = "TRG", scale = T)

compareClonotypes(dg0d_split_sch, 
                  numbers = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial")

compareClonotypes(dg0d_split_sch, 
                  numbers = 10, 
                  samples = c("BAL_Week_4", "BAL_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial")

 compareClonotypes(dg0d_split_sch, 
                  numbers = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "BAL_Week_4"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial") +
    scale_x_discrete(limits = c("PBMC_Week_0", "PBMC_Week_4", "BAL_Week_4"))

compareClonotypes(dg0d_split_sch, 
                  numbers = 10, 
                  samples = c("BAL_Week_4", "PBMC_Week_4", "BAL_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial") +
  scale_x_discrete(limits = c("BAL_Week_4", "PBMC_Week_4", "BAL_Week_8"))

clonalProportion(dg0d_split_sch_pbmc, cloneCall = "nt", 
                 split = c(10, 1e+05))

clonalDiversity(dg0d_split_sch_pbmc, cloneCall = "nt",
                group.by = "sample")



## DG0H ##
dg0h_split_sch_pbmc <- dg0h_split_sch[c(3, 4, 5)]

quantContig(dg0h_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRD", scale = T)
quantContig(dg0h_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRG", scale = T)
quantContig(dg0h_split_sch_pbmc, cloneCall = "gene+nt", chain = "both", scale = T)

abundanceContig(dg0h_split_sch_pbmc, cloneCall = "nt", scale = T)

lengthContig(dg0h_split_sch_pbmc, cloneCall = "nt", chain = "TRD", scale = T)
lengthContig(dg0h_split_sch_pbmc, cloneCall = "nt", chain = "TRG", scale = T)

compareClonotypes(dg0h_split_sch, 
                  numbers = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial")

compareClonotypes(dg0h_split_sch, 
                  numbers = 10, 
                  samples = c("BAL_Week_4", "BAL_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial")

compareClonotypes(dg0h_split_sch, 
                  numbers = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "BAL_Week_4"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial") +
  scale_x_discrete(limits = c("PBMC_Week_0", "PBMC_Week_4", "BAL_Week_4"))

compareClonotypes(dg0h_split_sch, 
                  numbers = 10, 
                  samples = c("BAL_Week_4", "PBMC_Week_4", "BAL_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  graph = "alluvial") +
  scale_x_discrete(limits = c("BAL_Week_4", "PBMC_Week_4", "BAL_Week_8"))

clonalProportion(dg0h_split_sch_pbmc, cloneCall = "nt", 
                 split = c(10, 1e+05))

clonalDiversity(dg0h_split_sch_pbmc, cloneCall = "nt",
                group.by = "sample")

## V113 ##
v113_split_sch_pbmc <- v113_split_sch[c(3, 4, 5)]

quantContig(v113_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRD", scale = T)
quantContig(v113_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRG", scale = T)
quantContig(v113_split_sch_pbmc, cloneCall = "gene+nt", chain = "both", scale = T)

abundanceContig(v113_split_sch_pbmc, cloneCall = "nt", scale = T)

lengthContig(v113_split_sch_pbmc, cloneCall = "nt", chain = "TRD", scale = T)
lengthContig(v113_split_sch_pbmc, cloneCall = "nt", chain = "TRG", scale = T)


v113_p <- clonalCompare(v113_split_sch, 
                  top.clones = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  palette = 'viridis')



clonalProportion(v113_split_sch_pbmc, cloneCall = "nt", 
                 split = c(10, 1e+05))






## V015 ##
v015_split_sch_pbmc <- v015_split_sch[c(3, 4, 5)]

quantContig(v015_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRD", scale = T)
quantContig(v015_split_sch_pbmc, cloneCall = "gene+nt", chain = "TRG", scale = T)
quantContig(v015_split_sch_pbmc, cloneCall = "gene+nt", chain = "both", scale = T)

abundanceContig(v015_split_sch_pbmc, cloneCall = "nt", scale = T)

lengthContig(v015_split_sch_pbmc, cloneCall = "nt", chain = "TRD", scale = T)
lengthContig(v015_split_sch_pbmc, cloneCall = "nt", chain = "TRG", scale = T)


v015_p <- clonalCompare(v015_split_sch, 
                  top.clones = 10, 
                  samples = c("PBMC_Week_0", "PBMC_Week_4", "PBMC_Week_8"), 
                  cloneCall="aa",
                  chain = "TRD",
                  palette = 'viridis')


clonalProportion(v015_split_sch_pbmc, cloneCall = "nt", 
                 split = c(10, 1e+05))



## FIGURE 6B
ggarrange(v113_p, v015_p, ncol = 2, nrow = 1)
ggsave('out/images/tcr_alluvial.pdf', width = 9, height = 3, units = 'in')




##############################################
## Analyze the data outside of scRepertoire ##
##############################################

#Read in the combined RDS file and tcr files
ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
ivbcg_all <- readRDS(file = "out/seurat_objects/ivbcg_integrated.rds")
dg0d_tcr <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_tcr <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_tcr <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_tcr <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

## Generate an annotated table of all GD TCRs
ivbcg_tcr_all <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr) 
barcodes <- FetchData(object = ivbcg_all, vars = c("seurat_clusters", "Timepoint", "Tissue"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr_all <- merge(x = barcodes, y = ivbcg_tcr_all, by = "barcode")

ivbcg_tcr_all$TRD = substr(ivbcg_tcr_all$TCR1, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
ivbcg_tcr_all$TRG = substr(ivbcg_tcr_all$TCR2, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)
ivbcg_tcr_all$paired_tcr <- paste(ivbcg_tcr_all$TRD, ivbcg_tcr_all$TRG, sep = "_")

write.csv(ivbcg_tcr_all, 'out/tcr/ivbcg_tcr_Supp2.csv')


#Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr)

#Add the metadata to each barcode
barcodes <- FetchData(object = ivbcg_gd, vars = c("seurat_clusters", "Timepoint", "Tissue"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr <- merge(x = barcodes, y = ivbcg_tcr, by = "barcode")






#Unmerge the TCR CSVs
dg0d_t <- ivbcg_tcr %>% filter(sample == "dg0d")
dg0h_t <- ivbcg_tcr %>% filter(sample == "dg0h")
v113_t <- ivbcg_tcr %>% filter(sample == "v113")
v015_t <- ivbcg_tcr %>% filter(sample == "v015")


## For each PTID, calculate the frequency of each TCR chain pairing in each group ##
#DG0D
dg0d_t$TRD = substr(dg0d_t$TCR2, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
dg0d_t$TRG = substr(dg0d_t$TCR1, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)

dg0d_t <- dg0d_t %>% #Within each cluster, calculate the % and absolute number of cells expressing each TCR gene pairing.
  group_by(Timepoint) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  group_by(Timepoint, TRD, TRG) %>%
  mutate(
    TCR_size = n(),
    TCR_pct = TCR_size/cluster_size
  ) %>%
  ungroup () %>%
  arrange(Timepoint)

dg0d_freq <- dg0d_t %>% subset(select = c(Timepoint, TRD, TRG, TCR_pct)) #Generate a new df containing only the TCR gene assignments, counts, and frequencies for each cluster.
dg0d_freq$paired_tcr <- paste(dg0d_freq$TRD, dg0d_freq$TRG, sep = "_") #Generate a new column listing the paired TCR genes
dg0d_freq <- distinct(dg0d_freq) #Remove the redundant rows
dg0d_freq$TCR_pct <- dg0d_freq$TCR_pct*100 # Multiply by 100 to get percentages


#DG0H
dg0h_t$TRD = substr(dg0h_t$TCR2, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
dg0h_t$TRG = substr(dg0h_t$TCR1, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)

dg0h_t <- dg0h_t %>% #Within each cluster, calculate the % and absolute number of cells expressing each TCR gene pairing.
  group_by(Timepoint) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  group_by(Timepoint, TRD, TRG) %>%
  mutate(
    TCR_size = n(),
    TCR_pct = TCR_size/cluster_size
  ) %>%
  ungroup () %>%
  arrange(Timepoint)

dg0h_freq <- dg0h_t %>% subset(select = c(Timepoint, TRD, TRG, TCR_pct)) #Generate a new df containing only the TCR gene assignments, counts, and frequencies for each cluster.
dg0h_freq$paired_tcr <- paste(dg0h_freq$TRD, dg0h_freq$TRG, sep = "_") #Generate a new column listing the paired TCR genes
dg0h_freq <- distinct(dg0h_freq) #Remove the redundant rows
dg0h_freq$TCR_pct <- dg0h_freq$TCR_pct*100 # Multiply by 100 to get percentages


#v113
v113_t$TRD = substr(v113_t$TCR2, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
v113_t$TRG = substr(v113_t$TCR1, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)

v113_t <- v113_t %>% #Within each cluster, calculate the % and absolute number of cells expressing each TCR gene pairing.
  group_by(Timepoint) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  group_by(Timepoint, TRD, TRG) %>%
  mutate(
    TCR_size = n(),
    TCR_pct = TCR_size/cluster_size
  ) %>%
  ungroup () %>%
  arrange(Timepoint)

v113_freq <- v113_t %>% subset(select = c(Timepoint, TRD, TRG, TCR_pct)) #Generate a new df containing only the TCR gene assignments, counts, and frequencies for each cluster.
v113_freq$paired_tcr <- paste(v113_freq$TRD, v113_freq$TRG, sep = "_") #Generate a new column listing the paired TCR genes
v113_freq <- distinct(v113_freq) #Remove the redundant rows
v113_freq$TCR_pct <- v113_freq$TCR_pct*100 # Multiply by 100 to get percentages


#V015
v015_t$TRD = substr(v015_t$TCR2, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
v015_t$TRG = substr(v015_t$TCR1, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)

v015_t <- v015_t %>% #Within each cluster, calculate the % and absolute number of cells expressing each TCR gene pairing.
  group_by(Timepoint) %>%
  mutate(cluster_size = n()) %>%
  ungroup() %>%
  group_by(Timepoint, TRD, TRG) %>%
  mutate(
    TCR_size = n(),
    TCR_pct = TCR_size/cluster_size
  ) %>%
  ungroup () %>%
  arrange(Timepoint)

v015_freq <- v015_t %>% subset(select = c(Timepoint, TRD, TRG, TCR_pct)) #Generate a new df containing only the TCR gene assignments, counts, and frequencies for each cluster.
v015_freq$paired_tcr <- paste(v015_freq$TRD, v015_freq$TRG, sep = "_") #Generate a new column listing the paired TCR genes
v015_freq <- distinct(v015_freq) #Remove the redundant rows
v015_freq$TCR_pct <- v015_freq$TCR_pct*100 # Multiply by 100 to get percentages


## For each PTID, generate pie charts describing the TRD and TRG gene usage at each timepoint ##
#DG0D gamma
dg0d_fg <- dg0d_freq
dg0d_fg <- dg0d_fg %>% subset(select = c(Timepoint, TRG, TCR_pct)) # Drop the unnecessary columns

pct_vg1 <- dg0d_fg %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV1") %>%
  summarize(Vg1 = sum(TCR_pct))

pct_vg2 <- dg0d_fg %>% # Calculate the percentage of Vg2 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV2") %>%
  summarize(Vg2 = sum(TCR_pct))

pct_vg3 <- dg0d_fg %>% # Calculate the percentage of Vg3 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV3") %>%
  summarize(Vg3 = sum(TCR_pct))

pct_vg8 <- dg0d_fg %>% # Calculate the percentage of Vg8 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV8") %>%
  summarize(Vg8 = sum(TCR_pct))

pct_vg9 <- dg0d_fg %>% # Calculate the percentage of Vg9 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV9") %>%
  summarize(Vg9 = sum(TCR_pct))

tcr_g_abbrev <- pct_vg1 %>% full_join(pct_vg2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vg3, by = "Timepoint") %>%
  full_join(pct_vg8, by = "Timepoint") %>%
  full_join(pct_vg9, by = "Timepoint")

tcr_g_abbrev[is.na(tcr_g_abbrev)] <- 0 # Convert NA to 0

tcr_g_abbrev <- tcr_g_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_g_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_g_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_g_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)


#DG0D delta
dg0d_fd <- dg0d_freq
dg0d_fd <- dg0d_fd %>% subset(select = c(Timepoint, TRD, TCR_pct)) # Drop the unnecessary columns

pct_vd1 <- dg0d_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV1") %>%
  summarize(Vd1 = sum(TCR_pct))

pct_vd2 <- dg0d_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV2") %>%
  summarize(Vd2 = sum(TCR_pct))

pct_vd3 <- dg0d_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV3") %>%
  summarize(Vd3 = sum(TCR_pct))

tcr_d_abbrev <- pct_vd1 %>% full_join(pct_vd2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vd3, by = "Timepoint")

tcr_d_abbrev[is.na(tcr_d_abbrev)] <- 0 # Convert NA to 0

tcr_d_abbrev <- tcr_d_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_d_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_d_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_d_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)



#DG0H gamma
dg0h_fg <- dg0h_freq
dg0h_fg <- dg0h_fg %>% subset(select = c(Timepoint, TRG, TCR_pct)) # Drop the unnecessary columns

pct_vg1 <- dg0h_fg %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV1") %>%
  summarize(Vg1 = sum(TCR_pct))

pct_vg2 <- dg0h_fg %>% # Calculate the percentage of Vg2 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV2") %>%
  summarize(Vg2 = sum(TCR_pct))

pct_vg3 <- dg0h_fg %>% # Calculate the percentage of Vg3 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV3") %>%
  summarize(Vg3 = sum(TCR_pct))

pct_vg8 <- dg0h_fg %>% # Calculate the percentage of Vg8 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV8") %>%
  summarize(Vg8 = sum(TCR_pct))

pct_vg9 <- dg0h_fg %>% # Calculate the percentage of Vg9 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV9") %>%
  summarize(Vg9 = sum(TCR_pct))

tcr_g_abbrev <- pct_vg1 %>% full_join(pct_vg2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vg3, by = "Timepoint") %>%
  full_join(pct_vg8, by = "Timepoint") %>%
  full_join(pct_vg9, by = "Timepoint")

tcr_g_abbrev[is.na(tcr_g_abbrev)] <- 0 # Convert NA to 0

tcr_g_abbrev <- tcr_g_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_g_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_g_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_g_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)


#DG0H delta
dg0h_fd <- dg0h_freq
dg0h_fd <- dg0h_fd %>% subset(select = c(Timepoint, TRD, TCR_pct)) # Drop the unnecessary columns

pct_vd1 <- dg0h_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV1") %>%
  summarize(Vd1 = sum(TCR_pct))

pct_vd2 <- dg0h_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV2") %>%
  summarize(Vd2 = sum(TCR_pct))

pct_vd3 <- dg0h_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV3") %>%
  summarize(Vd3 = sum(TCR_pct))

pct_vd4 <- dg0h_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV4") %>%
  summarize(Vd4 = sum(TCR_pct))

tcr_d_abbrev <- pct_vd1 %>% full_join(pct_vd2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vd3, by = "Timepoint") %>%
  full_join(pct_vd4, by = "Timepoint")

tcr_d_abbrev[is.na(tcr_d_abbrev)] <- 0 # Convert NA to 0

tcr_d_abbrev <- tcr_d_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_d_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_d_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_d_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)




#V113 gamma
v113_fg <- v113_freq
v113_fg <- v113_fg %>% subset(select = c(Timepoint, TRG, TCR_pct)) # Drop the unnecessary columns

pct_vg1 <- v113_fg %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV1") %>%
  summarize(Vg1 = sum(TCR_pct))

pct_vg2 <- v113_fg %>% # Calculate the percentage of Vg2 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV2") %>%
  summarize(Vg2 = sum(TCR_pct))

pct_vg3 <- v113_fg %>% # Calculate the percentage of Vg3 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV3") %>%
  summarize(Vg3 = sum(TCR_pct))

pct_vg8 <- v113_fg %>% # Calculate the percentage of Vg8 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV8") %>%
  summarize(Vg8 = sum(TCR_pct))

pct_vg9 <- v113_fg %>% # Calculate the percentage of Vg9 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV9") %>%
  summarize(Vg9 = sum(TCR_pct))

tcr_g_abbrev <- pct_vg1 %>% full_join(pct_vg2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vg3, by = "Timepoint") %>%
  full_join(pct_vg8, by = "Timepoint") %>%
  full_join(pct_vg9, by = "Timepoint")

tcr_g_abbrev[is.na(tcr_g_abbrev)] <- 0 # Convert NA to 0

tcr_g_abbrev <- tcr_g_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_g_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_g_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_g_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)


#V113 delta
v113_fd <- v113_freq
v113_fd <- v113_fd %>% subset(select = c(Timepoint, TRD, TCR_pct)) # Drop the unnecessary columns

pct_vd1 <- v113_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV1") %>%
  summarize(Vd1 = sum(TCR_pct))

pct_vd2 <- v113_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV2") %>%
  summarize(Vd2 = sum(TCR_pct))

pct_vd3 <- v113_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV3") %>%
  summarize(Vd3 = sum(TCR_pct))

pct_vd4 <- v113_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV4") %>%
  summarize(Vd4 = sum(TCR_pct))

tcr_d_abbrev <- pct_vd1 %>% full_join(pct_vd2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vd3, by = "Timepoint") %>%
  full_join(pct_vd4, by = "Timepoint")

tcr_d_abbrev[is.na(tcr_d_abbrev)] <- 0 # Convert NA to 0

tcr_d_abbrev <- tcr_d_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_d_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_d_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_d_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)




#V015 gamma
v015_fg <- v015_freq
v015_fg <- v015_fg %>% subset(select = c(Timepoint, TRG, TCR_pct)) # Drop the unnecessary columns

pct_vg1 <- v015_fg %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV1") %>%
  summarize(Vg1 = sum(TCR_pct))

pct_vg2 <- v015_fg %>% # Calculate the percentage of Vg2 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV2") %>%
  summarize(Vg2 = sum(TCR_pct))

pct_vg3 <- v015_fg %>% # Calculate the percentage of Vg3 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV3") %>%
  summarize(Vg3 = sum(TCR_pct))

pct_vg8 <- v015_fg %>% # Calculate the percentage of Vg8 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV8") %>%
  summarize(Vg8 = sum(TCR_pct))

pct_vg9 <- v015_fg %>% # Calculate the percentage of Vg9 in each cluster
  group_by(Timepoint) %>%
  filter(TRG == "TRGV9") %>%
  summarize(Vg9 = sum(TCR_pct))

tcr_g_abbrev <- pct_vg1 %>% full_join(pct_vg2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vg3, by = "Timepoint") %>%
  full_join(pct_vg8, by = "Timepoint") %>%
  full_join(pct_vg9, by = "Timepoint")

tcr_g_abbrev[is.na(tcr_g_abbrev)] <- 0 # Convert NA to 0

tcr_g_abbrev <- tcr_g_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_g_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_g_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_g_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)


#V015 delta
v015_fd <- v015_freq
v015_fd <- v015_fd %>% subset(select = c(Timepoint, TRD, TCR_pct)) # Drop the unnecessary columns

pct_vd1 <- v015_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV1") %>%
  summarize(Vd1 = sum(TCR_pct))

pct_vd2 <- v015_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV2") %>%
  summarize(Vd2 = sum(TCR_pct))

pct_vd3 <- v015_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV3") %>%
  summarize(Vd3 = sum(TCR_pct))

pct_vd4 <- v015_fd %>% # Calculate the percentage of Vg1 in each cluster
  group_by(Timepoint) %>%
  filter(TRD == "TRDV4") %>%
  summarize(Vd4 = sum(TCR_pct))

tcr_d_abbrev <- pct_vd1 %>% full_join(pct_vd2, by = "Timepoint") %>% # Join the tables
  full_join(pct_vd3, by = "Timepoint") %>%
  full_join(pct_vd4, by = "Timepoint")

tcr_d_abbrev[is.na(tcr_d_abbrev)] <- 0 # Convert NA to 0

tcr_d_abbrev <- tcr_d_abbrev %>% pivot_longer(!Timepoint, names_to = "tcr", values_to = "freq")

wk0 <- tcr_d_abbrev %>% filter(Timepoint == "Week_0") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk4 <- tcr_d_abbrev %>% filter(Timepoint == "Week_4") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

wk8 <- tcr_d_abbrev %>% filter(Timepoint == "Week_8") %>%
  ggplot(aes(x = "", y = freq, fill = tcr)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar ("y", start = 0) +
  theme_void()

ggarrange(wk0, wk4, wk8, ncol = 3, nrow = 1, common.legend = TRUE)



####################################################################
## Define the top expanded clones in each PTID at each time point ##
####################################################################

ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
dg0d_p <- read.csv(file= "out/tcr/dg0d_tcr_paired.csv")
dg0d_sch <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_p <- read.csv(file= "out/tcr/dg0h_tcr_paired.csv")
dg0h_sch <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_p <- read.csv(file= "out/tcr/v113_tcr_paired.csv")
v113_sch <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_p <- read.csv(file= "out/tcr/v015_tcr_paired.csv")
v015_sch <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

#Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_sch, dg0h_sch, v113_sch, v015_sch)

#Add the metadata to each barcode
barcodes <- FetchData(object = ivbcg_gd, vars = c("Timepoint", "Tissue", "seurat_clusters"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr <- merge(x = barcodes, y = ivbcg_tcr, by = "barcode")

#Add a new column with the paired-chain gene usage
ivbcg_tcr$TRD = substr(ivbcg_tcr$TCR1, 1, 5) #List the first five characters from the TCR1 column (TRD gene assignment)
ivbcg_tcr$TRG = substr(ivbcg_tcr$TCR2, 1, 5) #List the first five characters from the TCR2 column (TRG gene assignment)
ivbcg_tcr$paired_tcr <- paste(ivbcg_tcr$TRD, ivbcg_tcr$TRG, sep = "_")

#Enumerate the TCRs analyzed in PBMC and BAL
ivbcg_tcr %>% filter(Tissue == 'PBMC') %>% pull(barcode) %>% length()
ivbcg_tcr %>% filter(Tissue == 'BAL') %>% pull(barcode) %>% length()

## Measure the counts of each unique CDR3 in PBMC at each timepoint ##

# DG0D
counts <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_0", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0d_0 <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_0", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0d_0 <- dg0d_0 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct() 

counts <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_4", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0d_4 <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_4", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0d_4 <- dg0d_4 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_8", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0d_8 <- ivbcg_tcr %>% filter(sample == "dg0d", Timepoint == "Week_8", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0d_8 <- dg0d_8 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

dg0d_counts <- bind_rows(dg0d_0, dg0d_4, dg0d_8)


# DG0H
counts <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_0", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0h_0 <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_0", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0h_0 <- dg0h_0 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct() 

counts <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_4", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0h_4 <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_4", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0h_4 <- dg0h_4 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_8", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
dg0h_8 <- ivbcg_tcr %>% filter(sample == "dg0h", Timepoint == "Week_8", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
dg0h_8 <- dg0h_8 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

dg0h_counts <- bind_rows(dg0h_0, dg0h_4, dg0h_8)


# V113
counts <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_0", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v113_0 <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_0", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v113_0 <- v113_0 %>% subset(select = c(Timepoint, Tissue, CTaa, cdr3_aa2, Freq, paired_tcr)) %>% distinct() 

counts <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_4", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v113_4 <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_4", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v113_4 <- v113_4 %>% subset(select = c(Timepoint, Tissue, CTaa, cdr3_aa2, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_8", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v113_8 <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_8", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v113_8 <- v113_8 %>% subset(select = c(Timepoint, Tissue, CTaa, cdr3_aa2, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_4", Tissue == "BAL") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v113_4b <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_4", Tissue == "BAL") %>% full_join(counts, by = "CTaa")
v113_4b <- v113_4 %>% subset(select = c(Timepoint, Tissue, CTaa, cdr3_aa2, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_8", Tissue == "BAL") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v113_8b <- ivbcg_tcr %>% filter(sample == "v113", Timepoint == "Week_8", Tissue == "BAL") %>% full_join(counts, by = "CTaa")
v113_8b <- v113_8 %>% subset(select = c(Timepoint, Tissue, CTaa, cdr3_aa2, Freq, paired_tcr)) %>% distinct()

v113_counts <- bind_rows(v113_0, v113_4, v113_8, v113_4b, v113_8b)


# V015
counts <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_0", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v015_0 <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_0", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v015_0 <- v015_0 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct() 

counts <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_4", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v015_4 <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_4", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v015_4 <- v015_4 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

counts <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_8", Tissue == "PBMC") %>% pull(CTaa) %>% table() %>% as.data.frame()
counts <- rename(counts, CTaa = .)
v015_8 <- ivbcg_tcr %>% filter(sample == "v015", Timepoint == "Week_8", Tissue == "PBMC") %>% full_join(counts, by = "CTaa")
v015_8 <- v015_8 %>% subset(select = c(Timepoint, CTaa, Freq, paired_tcr)) %>% distinct()

v015_counts <- bind_rows(v015_0, v015_4, v015_8)


# Export the data

write.csv(dg0d_counts, file = "out/tcr/dg0d_counts.csv")
write.csv(dg0h_counts, file = "out/tcr/dg0h_counts.csv")
write.csv(v113_counts, file = "out/tcr/v113_counts.csv")
write.csv(v015_counts, file = "out/tcr/v015_counts.csv")

dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")

dg0d_ct <- dg0d_ct %>% subset(select = -c(X))
dg0h_ct <- dg0h_ct %>% subset(select = -c(X))
v113_ct <- v113_ct %>% subset(select = -c(X))
v015_ct <- v015_ct %>% subset(select = -c(X))



#####################################################

library(scales)

## Define the TCRs that are expanded after IV-BCG ##
dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")

#DG0H
dg0h_ct <- dg0h_ct %>% subset(select = -c(X))
dg0h_ct <- dg0h_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
dg0h_ct[is.na(dg0h_ct)] <- 0
dg0h_ct <- dg0h_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct = Week_0_count/sum(dg0h_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(dg0h_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(dg0h_ct$Week_8_count)*100)
dg0h_ct <- dg0h_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)

dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1,
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)

dg0h_exp <- dg0h_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf") %>% 
  filter(Week_4_count>2 | Week_8_count>2)

dg0h_pre <- dg0h_ct %>% filter(Week_0_count>1) %>%
  anti_join(dg0h_exp, by = 'CTaa')

dg0h_ct_new <- dg0h_ct %>%
  mutate(expanded = case_when(
    CTaa %in% dg0h_exp$CTaa ~ "Resp",
    CTaa %in% dg0h_pre$CTaa ~ "NR_exp",
    .default="NR_nonexp"
  ))

dg0h_ct_new %>% ggplot(aes(x = Week_0_pct_plus1, y = Week_4_pct_plus1, color = expanded)) +
  geom_point(size = 2) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)



#V113
v113_ct <- v113_ct %>% subset(select = -c(X))
v113_ct <- v113_ct %>% distinct()

v113_ct <- v113_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v113_ct$Week_0[is.na(v113_ct$Week_0)] <- 0
v113_ct$Week_4[is.na(v113_ct$Week_4)] <- 0
v113_ct$Week_8[is.na(v113_ct$Week_8)] <- 0

v113_ct <- v113_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
v113_ct <- v113_ct %>% mutate(Week_0_pct = Week_0_count/sum(v113_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v113_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v113_ct$Week_8_count)*100)
v113_ct <- v113_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)


v113_ct <- v113_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1,
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)

v113_ct %>% ggplot(aes(x = Week_0_pct_plus1, y = Week_4_pct_plus1)) +
  geom_point(size = 2) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans())



v113_exp <- v113_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf")
v113_exp <- v113_exp %>% filter(Week_4_count>2 | Week_8_count>2)

#V015
v015_ct <- v015_ct %>% subset(select = -c(X))
v015_ct <- v015_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v015_ct[is.na(v015_ct)] <- 0
v015_ct <- v015_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
v015_ct <- v015_ct %>% mutate(Week_0_pct = Week_0_count/sum(v015_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v015_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v015_ct$Week_8_count)*100)
v015_ct <- v015_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)

v015_ct %>% ggplot(aes(x = Week_0_pct, y = Week_4_pct)) +
  geom_point(size = 2) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans())

v015_exp <- v015_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf")
v015_exp <- v015_exp %>% filter(Week_4_count>2 | Week_8_count>2)

# Save the list of TCRs
write.csv(dg0h_exp, "out/tcr/dg0h_expanded.csv")
write.csv(v015_exp, "out/tcr/v015_expanded.csv")
write.csv(v113_exp, "out/tcr/v113_expanded.csv")



## Define expanded TCRs that did not increase after IV-BCG ##
# Filter to clones that were counted more than once before IV-BCG
dg0h_pre <- dg0h_ct %>% filter(Week_0_count>1)
v015_pre <- v015_ct %>% filter(Week_0_count>1)
v113_pre <- v113_ct %>% filter(Week_0_count>1)

# Exclude TCRs that expand after IV-BCG
dg0h_pre <- dg0h_pre %>% anti_join(dg0h_exp, by = 'CTaa')
v015_pre <- v015_pre %>% anti_join(v015_exp, by = 'CTaa')
v113_pre <- v113_pre %>% anti_join(v113_exp, by = 'CTaa')

# Save the list of TCRs
write.csv(dg0h_pre, "out/tcr/dg0h_pre_exp.csv")
write.csv(v015_pre, "out/tcr/v015_pre_exp.csv")
write.csv(v113_pre, "out/tcr/v113_pre_exp.csv")


## See whether expanded TCRs are present in BAL ##
#Read in the combined RDS file and tcr files
ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
dg0d_tcr <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_tcr <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_tcr <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_tcr <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

#Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr)

#Add the metadata to each barcode
barcodes <- FetchData(object = ivbcg_gd, vars = c("seurat_clusters", "Timepoint", "Tissue"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr <- merge(x = barcodes, y = ivbcg_tcr, by = "barcode")

#Unmerge the TCR CSVs
dg0d_t <- ivbcg_tcr %>% filter(sample == "dg0d")
dg0h_t <- ivbcg_tcr %>% filter(sample == "dg0h")
v113_t <- ivbcg_tcr %>% filter(sample == "v113")
v015_t <- ivbcg_tcr %>% filter(sample == "v015")


dg0h_exp <- dg0h_t %>% filter(CTaa == "CAGWDSTGWIKIF_CALWSLWGGYFTAQLFF" | CTaa == "CAGWDVRYKKLF_CALWTYPVGPIKLIF")
v113_exp <- v113_t %>% filter(CTaa == "CANWGNYYKKLF_CATFNGRDYWWNTDKLIF" | CTaa == "CAGWGSSTWWIKKF_CALRVRNSYWWDPDKLIF")
v015_exp <- v015_t %>% filter(CTaa == "CANWDRRVNYYKKLF_CALIGLHQVLGVYPTAQLFF" | CTaa == "CANWDRNYYKKLF_CALWELEEYWGYWEFWYTDKLIF")



###################################
## Combine TCR and Seurat Object ##
###################################

# Read in the combined RDS file and tcr files
ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")
dg0d_tcr <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_tcr <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_tcr <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_tcr <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

# Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr)

# Generate a df containing barcodes and CDR3s
cdr3 <- ivbcg_tcr %>% subset(select = c("barcode", "CTaa", "sample"))

# Annotate the expanded and non-expanded TCRs using the CDR3s
exp1 <- read.csv("out/tcr/v113_expanded.csv")
exp1$sample <- "v113"
exp2 <- read.csv("out/tcr/v015_expanded.csv")
exp2$sample <- "v015"
exp3 <- read.csv("out/tcr/dg0h_expanded.csv")
exp3$sample <- "dg0h"
exp_all <- bind_rows(exp1, exp2, exp3) %>%
  subset(select = c("CTaa", "sample"))
exp_all$expanded <- "postvax"

pre1 <- read.csv("out/tcr/dg0h_pre_exp.csv")
pre1$sample <- "dg0h"
pre2 <- read.csv("out/tcr/v015_pre_exp.csv")
pre2$sample <- "v015"
pre3 <- read.csv("out/tcr/v113_pre_exp.csv")
pre3$sample <- "v113"
pre_all <- bind_rows(pre1, pre2, pre3) %>%
  subset(select = c("CTaa", "sample"))
pre_all$expanded <- "prevax"

# Check that the pre- and post- expanded TCRs are non-overlapping
test <- anti_join(exp_all, pre_all, by = "CTaa")
test2 <- anti_join(pre_all, exp_all, by = "CTaa")

# Annotate the list of all TCRs as either expanded pre-vax, post-vax, or neither
cdr3_anno <- merged_df <- cdr3 %>%
  left_join(exp_all, by = c("CTaa", "sample")) %>%
  left_join(pre_all, by = c("CTaa", "sample"), suffix = c("_post", "_pre")) %>%
  mutate(expanded = coalesce(expanded_post, expanded_pre),
         expanded = ifelse(is.na(expanded), "no", expanded)) %>%
  select(barcode, CTaa, sample, expanded)
cdr3_anno <- cdr3_anno %>% remove_rownames %>% column_to_rownames(var="barcode")

# Combine with the RDS file and save
ivbcg_gd <- AddMetaData(ivbcg_gd, cdr3_anno)
saveRDS(ivbcg_gd, file = "out/seurat_objects/ivbcg_vd1_integrated_2024.rds")

# Visualize the expanded clonotypes
Idents(ivbcg_gd) <- "Tissue"
ivbcg_gd_p <- ivbcg_gd %>% subset(idents = c('PBMC'))

Idents(ivbcg_gd_p) <- "expanded"
DimPlot(ivbcg_gd_p, reduction = "umap", label = F, 
        cols = c('no' = "azure2", 'prevax' = "black", 'postvax' = "cornflowerblue"), 
        order = c("postvax", "prevax", "no")
        )
ggsave('out/images/Fig5D.pdf', width = 4.5, height = 3, units = 'in')

summary <- ivbcg_gd_p@meta.data %>% 
  group_by(expanded, seurat_clusters) %>%
  summarize(n= n()) %>%
  mutate(prop = n/sum(n)*100) %>%
  ungroup()

summary2 <- ivbcg_gd_p@meta.data %>%
  group_by(expanded) %>%
  summarize(n = n())

ggplot(ivbcg_gd_p@meta.data, aes(x=expanded, fill=seurat_clusters)) + 
  geom_bar(position = 'fill') +
  scale_fill_viridis_d(option = "turbo")
ggsave('out/images/Fig5  E.pdf', width = 4, height = 4, units = 'in')


## Generate a volcano plot displaying the DE genes between pre-expanded and BCG-expanded clonotypes
# Find DE genes between Pre and Post
Idents(ivbcg_gd) <- "expanded"
ivbcg_gd <- ivbcg_gd %>% subset(idents = c('prevax', 'postvax'))
markers <- FindAllMarkers(ivbcg_gd)
markers <- markers %>% filter(!grepl(c('ENSMMUG|TRA|TRB'), gene)) %>% filter(cluster == "postvax")

# Make the volcano plot
markers$diffexp <- "NO"
markers$diffexp[markers$avg_log2FC> 1 & markers$p_val_adj < 0.05] <- "UP"
markers$diffexp[markers$avg_log2FC< -1 & markers$p_val_adj < 0.05] <- "DOWN"

markers$delabel <- NA
markers$delabel[markers$diffexp != "NO"] <- markers$gene[markers$diffexp != "NO"]

markers %>% ggplot(aes(x = avg_log2FC, y = -log10(p_val_adj), col = diffexp, label = delabel)) + 
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), col = "azure4") +
  geom_hline(yintercept = -log10(0.05), col = "black") +
  scale_color_manual(values = c("blue", "gray", "red")) +
  geom_text_repel(label = markers$delabel,
                  point.padding = 0.25)
ggsave('out/images/Fig6F.pdf', width = 4.5, height = 3, units = 'in')



######################################
## Look for matches with known TCRs ##
######################################

# Read in the merged TCR df
#Read in the combined RDS file and tcr files
ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
dg0d_tcr <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_tcr <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_tcr <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_tcr <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

#Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr)

#Add the metadata to each barcode
barcodes <- FetchData(object = ivbcg_gd, vars = c("seurat_clusters", "Timepoint", "Tissue"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr <- merge(x = barcodes, y = ivbcg_tcr, by = "barcode")

#Search the data set for known TCRs
known_tcr_ref <- read.csv(file = "known_tcrs_for_conga_alleles.csv")
known_tcr_ref <- known_tcr_ref %>% rename (cdr3_aa1 = cdr3a, cdr3_aa2 = cdr3b, CTaa = CDR3_paired)
known_tcr_ref$CTaa <- paste(known_tcr_ref$cdr3_aa1, known_tcr_ref$cdr3_aa2, sep = "_")

tcr_known_delta <- inner_join(ivbcg_tcr, known_tcr_ref, by = "cdr3_aa2")
tcr_known_gamma <- inner_join(ivbcg_tcr, known_tcr_ref, by = "cdr3_aa1")
tcr_known_paired <- inner_join(ivbcg_tcr, known_tcr_ref, by = "CTaa") #Note that only gamma matches are found.


tcr_known_gamma <- tcr_known_gamma %>% subset(select = c(barcode, seurat_clusters, CTaa.x, CTaa.y, cdr3_aa1, dataset, epitope))
tcr_known_gamma <- tcr_known_gamma %>% rename(cdr3_infant = CTaa.x, cdr3_ref = CTaa.y) #Trim down the df and rename the columns for clarity

known_gamma_summary <- tcr_known_gamma %>% group_by(cdr3_aa1, epitope, seurat_clusters) %>% summarize(count = n())
known_gamma_summary <- known_gamma_summary %>% group_by(cdr3_aa1) %>% mutate(sum = sum(count)) #df with the number of matches per cluster

known_gamma_summary2 <- known_gamma_summary %>% subset(select = c(cdr3_aa1, epitope, sum)) %>% distinct() #df with the total number of matches

write.csv(known_gamma_summary, file = "R_files/TCR/Out/gamma_matches.csv")
write.csv(known_gamma_summary2, file = "R_files/TCR/Out/gamma_matches_sums.csv")






##############################################
## Compare the Gini Index pre- and post-Vax ##
##############################################

library(DescTools)
library(ggpubr)
library(vegan)
library(dplyr)
library(viridis)


#Import the TCR counts, annotate by PTID, and merge
dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")

dg0d_ct <- dg0d_ct %>% subset(select = -c(X))
dg0h_ct <- dg0h_ct %>% subset(select = -c(X))
v113_ct <- v113_ct %>% subset(select = -c(X))
v015_ct <- v015_ct %>% subset(select = -c(X))

dg0d_ct$PTID <- "DG0D"
dg0h_ct$PTID <- "DG0H"
v113_ct$PTID <- "V113"
v015_ct$PTID <- "V015"

counts <- bind_rows(dg0d_ct, dg0h_ct, v113_ct, v015_ct)


#Calculate the Gini Coefficient at each timepoint for each PTID
gini <- counts %>% 
  group_by(PTID, Timepoint) %>%
  summarize(gini = Gini(Freq, unbiased = FALSE))

# Plot the results
ggpaired(gini, x = "Timepoint", y = "gini", id = "PTID", 
         line.color = "black", 
         point.size = 2,
         line.size = 0.4, 
         fill = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
    ggtitle("Gini Coefficient over Time")
  

#Calculate and plot the Shannon Diversity Index
shannon <- counts %>% group_by(PTID, Timepoint) %>% summarize(shannon = diversity(Freq))

ggpaired(shannon, x = "Timepoint", y = "shannon", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Week_0", "Week_4"), c("Week_0", "Week_8"), c("Week_4", "Week_8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("Shannon Index over Time")


#Calculate and plot the Inverse Simpson Index
invsimp <- counts %>% group_by(PTID, Timepoint) %>% summarize(invsimp = diversity(Freq, "invsimpson"))

ggpaired(invsimp, x = "Timepoint", y = "invsimp", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Week_0", "Week_4"), c("Week_0", "Week_8"), c("Week_4", "Week_8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("Inverse Simpson Index over Time")


# Calculate the percent expanded clonotypes at each timepoint
expanded <- counts %>% group_by(PTID, Timepoint) %>%
  mutate(total_clonos = n()) %>%
  dplyr::filter(Freq >1) %>%
  mutate(total_exp = n(),
            pct_exp = total_exp/total_clonos*100)

expanded <- expanded %>% subset(select = c("Timepoint", "PTID", "total_clonos", "total_exp", "pct_exp")) %>% distinct()

expanded[nrow(expanded) + 1,] <- list("Week_8", "DG0D", 0,0,0)
expanded[nrow(expanded) + 1,] <- list("Week_8", "DG0H", 0,0,0)

pct_exp <- ggpaired(expanded, x = "Timepoint", y = "pct_exp", id = "PTID", 
         line.color = "black", 
         line.size = 1, 
         point.size = 1,
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("Percent Expanded Clonotypes over Time") 

pct_exp$layers <- pct_exp$layers[-1]
pct_exp
ggsave("out/images/fig6c.pdf", width = 4, height = 4, units = 'in')
