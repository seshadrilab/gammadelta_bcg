#Load dplyr, Seurat, and patchwork
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(scRepertoire)


# Read in the files
tcr_annotated <- read.csv(file = 'TCR/In/all_contig_annotations.csv')
prod <- read.csv(file = 'TCR/In/productive_contigs.csv')
gex <- readRDS(file = 'GEX/Out/pbmc_tcells_2024.rds')


###################################
## APPLY QC FILTERS TO TCR READS ##
###################################

# Starting counts
tcr_annotated$barcode %>% unique() %>% length ()
# 26052 cell barcodes and 49512 TCR contigs

# Create a new TCR dataframe containing TCR data from only the cells that passed GEX QC.
# Add the PTID metadata to the TCR df.
gex_md <- gex@meta.data %>%
  subset(select = c('ptid')) %>%
  rownames_to_column('barcode')
gex_md$gex_pass_qc <- 'true'
tcr_final <- merge(x=gex_md, y=tcr_annotated, by="barcode") #Merge the GEX QC'd barcodes with the TCR df using inner join.
tcr_final$barcode %>% unique() %>% length()
# 19,002 cell barcodes and 41,190 TCR contigs

# Keep only full length TCRs with a CDR3
tcr_final <- tcr_final %>% subset(!cdr3 == "None") #Remove entires with CDR3 = "None"
tcr_final <- tcr_final %>% subset(full_length == "True") #Remove entries that are not full length TCRs.
tcr_final$barcode %>% unique() %>% length()
# 18,968 cell barcodes and 40,011 TCR contigs

# Exclude TCRs with <3 UMIs or <20 reads
tcr_final <- tcr_final %>% filter(reads >19 & umis >2)
tcr_final$barcode %>% unique() %>% length()
# 17,321 cell barcodes and 28,522 TCR contigs

# Keep only productive TCRs
prod <- prod %>% subset(select = c('contig_id', 'productive'))
tcr_final <- tcr_final %>% subset(select = -c(productive))
tcr_final <- merge(x = tcr_final, y = prod, by = 'contig_id')
tcr_final$barcode %>% unique() %>% length()
# 14,302 cell barcodes and 21,591 TCR contigs

#Save the resulting dataframe
write.csv(tcr_final, "TCR/out/tcr_QCd_productive.csv")



###################################################
## DEFINE PAIRED-CHAIN TCR FOR EACH CELL BARCODE ##
###################################################

# Read in the TCR df
# Drop any entries with an unknown sample donor (ptid)
tcr <- read.csv(file = 'TCR/out/tcr_QCd_productive.csv')
tcr <- tcr %>% drop_na(ptid)
tcr$barcode %>% unique() %>% length()
# 12,197 cell barcodes and 17,588 TCR contigs

# Store the ptid corresponding to each barcode
ptids <- tcr %>% subset(select = c('barcode', 'ptid')) %>% distinct()

# Pull the column names from the built-in example data
data("contig_list")
ex <- (contig_list[[1]])
colorder <- names(ex)

# Reorder the columns of the tcr DF
# Note that this deletes the PTID column which breaks combineTCR
tcr <- tcr %>% select(colorder)

# Generate the contig list
contig_list_mdm <- list(infant = tcr)

## Using scRepertoire, reconfigure the TCR df such that each barcode appears in only one row and is assigned only one contig per chain.
# Use removeNA=TRUE to exclude barcodes without paired-chain data.
# Use removeMulti=FALSE to keep barcodes with multiple contigs per chain 
# Use filterMulti=TRUE to select the TCRg and TCRd with the highest expression for each barcode.
# combineTCR requires the TCRs in list format with each sample constituting an item in the list. Here, we only have one sample, so we assign a sham sample and ID to make our data compatible with combineTCR.
tcr_combined <- combineTCR(
  contig_list_mdm,
  removeNA = TRUE,
  removeMulti = FALSE,
  filterMulti = TRUE)

tcr_combined <- as.data.frame(tcr_combined)

#Clean-up the tcr_combined dataframe and save the file
tcr_combined <- tcr_combined %>% rename_all(~stringr::str_replace(.,"S1.","")) #Remove the prefix from the columns

#Re-annotate the barcodes in the tcr_combined df with the ptid numbers
tcr_combined <- left_join(tcr_combined, ptids, by = 'barcode')

## Save the final file.
write.csv(tcr_combined, "TCR/Out/tcr_combined.csv")


## For each cell in the "tcr" dataframe, determine whether it is Vd1, Vd3, Vg9Vd2, or Other_Vd2 (Vd2 not paired with Vg9).
tcr <- read.csv('TCR/Out/tcr_combined.csv')

# Create a new column with the paired TCR genes
tcr$TRG = substr(tcr$TCR1, 1, 5)
tcr$TRD = substr(tcr$TCR2, 1, 5)
tcr$paired_tcr <- paste(tcr$TRD, tcr$TRG, sep = "_")

# Annotate the TCR groups
tcr_vd1 <- tcr %>% filter(TRD == "TRDV1") %>% mutate(tcr_group = "Vd1")
tcr_vd3 <- tcr %>% filter(TRD == "TRDV3") %>% mutate(tcr_group = "Vd3") 
tcr_vg9vd2 <- tcr %>% filter(paired_tcr == "TRDV2_TRGV9") %>% mutate(tcr_group = "Vg9Vd2")
tcr_other_vd2 <- tcr %>% filter(TRD == "TRDV2" & !(paired_tcr == "TRDV2_TRGV9")) %>% mutate(tcr_group = "Other_Vd2")
tcr_grouped <- rbind(tcr_vd1, tcr_vg9vd2, tcr_other_vd2, tcr_vd3)


# Enumerate the percentage of cell barcodes in each tcr group
group_nos <- table(tcr_grouped$tcr_group) %>% as.data.frame()
group_nos <- group_nos %>% mutate(
  group_pct = Freq/sum(group_nos$Freq)
)

write.csv(tcr_grouped, file = "TCR/Out/tcr_grouped.csv") #Export

# Enumerate the percentage of Vg9-expressing T cells expressing Vd1, Vd2, and Vd3 chains
tcr_grouped <- read.csv('TCR/Out/tcr_grouped.csv')
tcr_vg9 <- tcr_grouped %>% filter(TRG == 'TRGV9')
tcr_not_vg9 <- tcr_grouped %>% filter(TRG != 'TRGV9')

table <- table(tcr_not_vg9$TRD) %>% as.data.frame()
sum = sum(table$Freq)
table <- table %>% mutate(
  pct = Freq/sum*100
)

