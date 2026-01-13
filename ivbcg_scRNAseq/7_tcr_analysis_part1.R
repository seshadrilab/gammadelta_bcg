library(dplyr)
library(scRepertoire)
library(Seurat)
library(circlize)
library(scales)
library(ggpubr)
library(DescTools)
library(vegan)
library(viridis)
library(ggrepel)
library(fgsea)
library(presto)
library(msigdbr)
library(tidyverse)


####################################################
## Define the TCRs that are expanded after IV-BCG ##
####################################################

# Read in the TCR counts for each donor
dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")


## For each donor, define the TCRs that expand in response to IV-BCG (expanded responding)
## The clonotype is expanded in response to vaccination if any of the following criteria are satisfied:
## 1. The TCR is expanded more than 3-fold from pre-vaccination to week 4
## 2. The TCR is expanded more than 3-fold from pre-vaccination to week 8
## 3. The TCR is absent pre-vaccination and detected at least 3 times at week 4
## 4. The TCR is absent pre-vaccination and detected at least 3 times at week 8

# Note that donor DG0D has insufficient cell counts to define these TCR groups

## DG0H
dg0h_ct <- dg0h_ct %>% subset(select = -c(X))
dg0h_ct <- dg0h_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
dg0h_ct[is.na(dg0h_ct)] <- 0
dg0h_ct <- dg0h_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
# Calculate the frequency (%) of each clonotype at each time point
dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct = Week_0_count/sum(dg0h_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(dg0h_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(dg0h_ct$Week_8_count)*100)
# Calculate the fold-change (ratio) of each clonotype at week 4 or 8 compared to pre-vaccination
dg0h_ct <- dg0h_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)
dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1, # Create a +1 column to allow log transformation later on
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)
# Define the expanded TCRs according to the criteria above
dg0h_exp <- dg0h_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf") 
dg0h_exp <- dg0h_exp %>% filter(Week_4_count>2 | Week_8_count>2)


## V113
v113_ct <- v113_ct %>% subset(select = -c(X))
v113_ct <- v113_ct %>% distinct()
v113_ct <- v113_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v113_ct[is.na(v113_ct)] <- 0
v113_ct <- v113_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
# Calculate the frequency (%) of each clonotype at each time point
v113_ct <- v113_ct %>% mutate(Week_0_pct = Week_0_count/sum(v113_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v113_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v113_ct$Week_8_count)*100)
# Calculate the fold-change (ratio) of each clonotype at week 4 or 8 compared to pre-vaccination
v113_ct <- v113_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)
v113_ct <- v113_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1, # Create a +1 column to allow log transformation later on
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)
# Define the expanded TCRs according to the criteria above
v113_exp <- v113_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf")
v113_exp <- v113_exp %>% filter(Week_4_count>2 | Week_8_count>2)


## V015
v015_ct <- v015_ct %>% subset(select = -c(X))
v015_ct <- v015_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v015_ct[is.na(v015_ct)] <- 0
v015_ct <- v015_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
# Calculate the frequency (%) of each clonotype at each time point
v015_ct <- v015_ct %>% mutate(Week_0_pct = Week_0_count/sum(v015_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v015_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v015_ct$Week_8_count)*100)
# Calculate the fold-change (ratio) of each clonotype at week 4 or 8 compared to pre-vaccination
v015_ct <- v015_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)
v015_ct <- v015_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1, # Create a +1 column to allow log transformation later on
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)
# Define the expanded TCRs according to the criteria above
v015_exp <- v015_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf")
v015_exp <- v015_exp %>% filter(Week_4_count>2 | Week_8_count>2)


# Save the lists of TCRs that expand in response to vaccination
write.csv(dg0h_exp, "out/tcr/dg0h_expanded_3x.csv")
write.csv(v015_exp, "out/tcr/v015_expanded_3x.csv")
write.csv(v113_exp, "out/tcr/v113_expanded_3x.csv")


## For each donor, define TCRs that were counted more than once before vaccination and did not respond to IV-BCG (expanded non-responding)
# Filter to clones that were counted more than once before IV-BCG
dg0h_pre <- dg0h_ct %>% filter(Week_0_count>1)
v015_pre <- v015_ct %>% filter(Week_0_count>1)
v113_pre <- v113_ct %>% filter(Week_0_count>1)

# Exclude TCRs that expand after IV-BCG
dg0h_pre <- dg0h_pre %>% anti_join(dg0h_exp, by = 'CTaa')
v015_pre <- v015_pre %>% anti_join(v015_exp, by = 'CTaa')
v113_pre <- v113_pre %>% anti_join(v113_exp, by = 'CTaa')

# Save the list of TCRs
write.csv(dg0h_pre, "out/tcr/dg0h_pre_exp_3x.csv")
write.csv(v015_pre, "out/tcr/v015_pre_exp_3x.csv")
write.csv(v113_pre, "out/tcr/v113_pre_exp_3x.csv")



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

# Generate a df containing barcodes and CDR3 sequences
cdr3 <- ivbcg_tcr %>% subset(select = c("barcode", "CTaa", "sample"))

# Read in and concatenate the CSVs containing the expanded_responding (expR) TCRs
exp1 <- read.csv("out/tcr/v113_expanded_3x.csv")
exp1$sample <- "v113"
exp2 <- read.csv("out/tcr/v015_expanded_3x.csv")
exp2$sample <- "v015"
exp3 <- read.csv("out/tcr/dg0h_expanded_3x.csv")
exp3$sample <- "dg0h"
exp_all <- bind_rows(exp1, exp2, exp3) %>%
  subset(select = c("CTaa", "sample"))
exp_all$expanded <- "postvax"

# Read in and concatenate the CSVs containing the expanded_NONresponding (expNR) TCRs
pre1 <- read.csv("out/tcr/dg0h_pre_exp_3x.csv")
pre1$sample <- "dg0h"
pre2 <- read.csv("out/tcr/v015_pre_exp_3x.csv")
pre2$sample <- "v015"
pre3 <- read.csv("out/tcr/v113_pre_exp_3x.csv")
pre3$sample <- "v113"
pre_all <- bind_rows(pre1, pre2, pre3) %>%
  subset(select = c("CTaa", "sample"))
pre_all$expanded <- "prevax"

# Check that the expR and expNR TCRs are non-overlapping
test <- anti_join(exp_all, pre_all, by = "CTaa")
test2 <- anti_join(pre_all, exp_all, by = "CTaa")

# Annotate the list of all CDR3s as either expanded pre-vax, expanded post-vax, or not expanded
# Note than in the manuscript, 'expanded_post'=expR, 'expanded_pre'=expNR, and 'no'= not expanded
cdr3_anno <- cdr3 %>%
  left_join(exp_all, by = c("CTaa", "sample")) %>%
  left_join(pre_all, by = c("CTaa", "sample"), suffix = c("_post", "_pre")) %>%
  mutate(expanded = coalesce(expanded_post, expanded_pre),
         expanded = ifelse(is.na(expanded), "no", expanded)) %>%
  select(barcode, CTaa, sample, expanded)
cdr3_anno <- cdr3_anno %>% remove_rownames %>% column_to_rownames(var="barcode")

# Add these metadata to the Seurat object
ivbcg_gd <- AddMetaData(ivbcg_gd, cdr3_anno)
saveRDS(ivbcg_gd, file = "out/seurat_objects/ivbcg_vd1_integrated_exp_3x.rds")


# Visualize the expanded/unexpanded clonotypes on the PBMC UMAP
Idents(ivbcg_gd) <- "Tissue"
ivbcg_gd_p <- ivbcg_gd %>% subset(idents = c('PBMC'))

Idents(ivbcg_gd_p) <- "expanded"
DimPlot(ivbcg_gd_p, reduction = "umap", label = F, 
        cols = c('no' = "azure2", 'prevax' = "black", 'postvax' = "cornflowerblue"), 
        order = c("postvax", "prevax", "no")
        )
ggsave('out/images/Fig5D_3x_reorder.pdf', width = 4.5, height = 3, units = 'in')


# Visualize a stacked bar plot showing the cluster frequencies for expR, expNR, and unexpanded clonotypea
ggplot(ivbcg_gd_p@meta.data, aes(x=expanded, fill=seurat_clusters)) + 
  geom_bar(position = 'fill') +
  scale_fill_viridis_d(option = "turbo")
ggsave('out/images/Fig5E_3x.pdf', width = 4, height = 4, units = 'in')


# Identify the differentially-expressed genes between expR and expNR clonotypes in PBMC
# Exclude genes with no gene symbol and alpha-beta TCR genes
Idents(ivbcg_gd_p) <- "expanded"
ivbcg_gd_p <- ivbcg_gd_p %>% subset(idents = c('prevax', 'postvax'))
markers <- FindAllMarkers(ivbcg_gd_p)
markers <- markers %>% filter(!grepl(c('ENSMMUG|TRA|TRB'), gene)) %>% filter(cluster == "postvax")

# Visualize the genes on a volcano plot
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
ggsave('out/images/Fig6F_3x_p.pdf', width = 4.5, height = 3, units = 'in')



################################################
## OVERALL % OF EXPANDED CLONOTYPES OVER TIME ##
################################################

#Import the TCR counts, annotate by sample donor, and merge
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
