# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)

####################
## SET UP TABLES ###
####################

## NOTE: The study animal IDs have been altered in the public SRA files.
## Please substitute the anonymized IDs (left) in the SRA downloads for the following animal IDs (right):
## SEDXX2: V113
## SEDXX5: V015
## SEDXX6: DG0D
## SEDXX7: DG0H
## 20221: 36327
## 4219: DF4P
## 5075: D12L
## 19821: MC30
## 4119: DF2C
## 22119: A14V139
## 5076: DIC4
## 5077: 17C231
## 8321: MF46
## 21519: HRP
## 21819: DGKM
## 4319: O8M
## 5078: P599
## 5079: 18C062
## 20021:36852
## 4519: DF1R
## 5078:P599
## 5079: 18C062
## 20021: 36852
## 4519: DF1R
## 21219: DGPR
## 5219: O4A
## 5080: DHZI
## 5081: 16C192
## 22619: DGHL
## 21919: DGFM
## 4819: OCC
## 5082: DIAV
## 5083: DI4R
## 19421: 36818
## 19221: MI12
## 19521: MB92
## 2618: 2618
## 2718: 2718
## 2918: 2918

# Read in main CSV
df <- read.csv("vd1_cytof.csv")

# Note, PTIDs DIC7, 17C333, and 18C136 include too few Vg9neg cells and must be excluded
df <- df %>% filter(!((PTID == "DIC7") | (PTID == "17C333") | (PTID == "18C136")))

# Filter the dataframe by unstimulated cells
us <- df %>% filter(Stim == "NS")

# Calculate the background-subtracted responses to stimulation
df_NS <- df %>% filter(Stim == "NS")
df_Mtb <- df %>% filter(Stim == "MtbL")
metadata <- df_Mtb[, 1:6]
df_bgsub <- df_Mtb[, 7:37] - df_NS[, 7:37]
df_bgsub <- cbind(metadata, df_bgsub)
df_bgsub <- df_bgsub %>% arrange(-desc(Tissue))
df_bgsub[df_bgsub < 0] <- 0

# Arrange the timepoints in the correct order
us$Timepoint <- as.factor(us$Timepoint)
levels(us$Timepoint) <- c("pre", "4", "8")

df_bgsub$Timepoint <- as.factor(df_bgsub$Timepoint)
levels(df_bgsub$Timepoint) <- c("pre", "4", "8")


######################
## SET UP FUNCTION ###
######################

# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the trajectory of all PTIDs over time as well as the median trajectory.
# p-values represent comparisons between wk0/wk4, and wk0/wk8

make_time_plt <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(Cell_type == subset_use)
  out <-
    ggpaired(
      new_df,
      x = "Timepoint",
      y = stat_use,
      id = "PTID",
      facet.by = "Tissue",
      point.size = 1.5,
      line.size = 0.5,
      short.panel.labs = FALSE
    ) +
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cells") +
    xlab("Timepoint") +
    rremove("legend")

  return(out)
}


########################
## unstimulated cells ##
########################
## Display the % expression of each marker over time
## Supplementary Figure 6A-F
cd69 <- make_time_plt(us, "Vd1", "CD69_Freq", "% CD69+")
ifng <- make_time_plt(us, "Vd1", "IFNg_Freq", "% IFN-gamma+")
tnf <- make_time_plt(us, "Vd1", "TNF_Freq", "% TNF+")
grb <- make_time_plt(us, "Vd1", "GranB_Freq", "% Granzyme B+")
grk <- make_time_plt(us, "Vd1", "GranK_Freq", "% Granzyme K+")
prf <- make_time_plt(us, "Vd1", "Perforin_Freq", "% Perforin+")

ggarrange(cd69, ifng, tnf, grb, grk, prf, nrow = 2, ncol = 3)

ggsave("cytof_unstim.pdf")


#############################
## RESPONSE TO STIMULATION ##
#############################
## Display the % expression of each marker over time
## Supplementary Figure 6G-L
cd69 <- make_time_plt(df_bgsub, "Vd1", "CD69_Freq", "% CD69+")
ifng <- make_time_plt(df_bgsub, "Vd1", "IFNg_Freq", "% IFN-gamma+")
tnf <- make_time_plt(df_bgsub, "Vd1", "TNF_Freq", "% TNF+")
grb <- make_time_plt(df_bgsub, "Vd1", "GranB_Freq", "% Granzyme B+")
grk <- make_time_plt(df_bgsub, "Vd1", "GranK_Freq", "% Granzyme K+")
prf <- make_time_plt(df_bgsub, "Vd1", "Perforin_Freq", "% Perforin+")

ggarrange(cd69, ifng, tnf, grb, grk, prf, nrow = 2, ncol = 3)

ggsave("cytof_mtbl.pdf")
