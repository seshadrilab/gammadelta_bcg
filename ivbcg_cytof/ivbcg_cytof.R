
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

# Read in main CSV
df <- read.csv("table_231207_vd1.csv")

# Note, PTIDs DIC7 and 17C333 include too few Vg9neg cells and must be excluded
df <- df %>% filter(!((PTID == "DIC7")|(PTID == "17C333")|(PTID == "18C136")))

# Filter the dataframe by unstimulated cells
us <- df %>% filter(Stim == 'NS')

# Calculate the background-subtracted responses to stimulation
df_NS <- df %>% filter(Stim == "NS")
df_Mtb <- df %>% filter(Stim == "MtbL")
metadata <- df_Mtb[,1:6]
df_bgsub <- df_Mtb[,7:37]-df_NS[,7:37]
df_bgsub <- cbind(metadata, df_bgsub)
df_bgsub <- df_bgsub %>% arrange(-desc(Tissue))
df_bgsub[df_bgsub < 0] <- 0

# Arrange the timepoints in the correct order
us$Timepoint <- as.factor(us$Timepoint)
levels(us$Timepoint) <- c('pre', '4', '8')

df_bgsub$Timepoint <- as.factor(df_bgsub$Timepoint)
levels(df_bgsub$Timepoint) <- c('pre', '4', '8')

#######################
## SET UP FUNCTIONS ###
#######################

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

cd69 <- make_time_plt(us, "Vd1", "CD69_Freq", "% CD69+")
ifng <- make_time_plt(us, "Vd1", "IFNg_Freq", "% IFN-gamma+")
tnf <- make_time_plt(us, "Vd1", "TNF_Freq", "% TNF+")
grb <- make_time_plt(us, "Vd1", "GranB_Freq", "% Granzyme B+")
grk <- make_time_plt(us, "Vd1", "GranK_Freq", "% Granzyme K+")
prf <- make_time_plt(us, "Vd1", "Perforin_Freq", "% Perforin+")

ggarrange(cd69, ifng, tnf, grb, grk, prf, nrow = 2, ncol = 3)

ggsave('cytof_unstim.pdf')

#############################
## RESPONSE TO STIMULATION ##
#############################

cd69 <- make_time_plt(df_bgsub, "Vd1", "CD69_Freq", "% CD69+")
ifng <- make_time_plt(df_bgsub, "Vd1", "IFNg_Freq", "% IFN-gamma+")
tnf <- make_time_plt(df_bgsub, "Vd1", "TNF_Freq", "% TNF+")
grb <- make_time_plt(df_bgsub, "Vd1", "GranB_Freq", "% Granzyme B+")
grk <- make_time_plt(df_bgsub, "Vd1", "GranK_Freq", "% Granzyme K+")
prf <- make_time_plt(df_bgsub, "Vd1", "Perforin_Freq", "% Perforin+")

ggarrange(cd69, ifng, tnf, grb, grk, prf, nrow = 2, ncol = 3)

ggsave('cytof_mtbl.pdf')


####################################################################
## HOW ENRICHED ARE MTB-RESPONSIVE CELLS IN BAL COMPARED TO PBMC? ##
####################################################################

cd69 <- pbmc_v_bal(df_bgsub, "Vd1", "CD69_Freq", "% CD69+")
hladr <- pbmc_v_bal(df_bgsub, "Vd1", "HLADR_Freq", "% HLA-DR+")
cd154 <- pbmc_v_bal(df_bgsub, "Vd1", "CD154_Freq", "% CD154+")
cd107a <- pbmc_v_bal(df_bgsub, "Vd1", "CD107a_Freq", "% CD107a+")
ccl4 <- pbmc_v_bal(df_bgsub, "Vd1", "CCL4_Freq", "% CCL4+")
cxcr3 <- pbmc_v_bal(df_bgsub, "Vd1", "CXCR3_Freq", "% CXCR3+")
ifng <- pbmc_v_bal(df_bgsub, "Vd1", "IFNg_Freq", "% IFN-gamma+")
tnf <- pbmc_v_bal(df_bgsub, "Vd1", "TNF_Freq", "% TNF+")
il2 <- pbmc_v_bal(df_bgsub, "Vd1", "IL2_Freq", "% IL-2+")
il10 <- pbmc_v_bal(df_bgsub, "Vd1", "IL10_Freq", "% IL-10+")
il17 <- pbmc_v_bal(df_bgsub, "Vd1", "IL17A_Freq", "% IL-17+")
grb <- pbmc_v_bal(df_bgsub, "Vd1", "GranB_Freq", "% Granzyme B+")
grk <- pbmc_v_bal(df_bgsub, "Vd1", "GranK_Freq", "% Granzyme K+")
prf <- pbmc_v_bal(df_bgsub, "Vd1", "Perforin_Freq", "% Perforin+")
nkg2a <- pbmc_v_bal(df_bgsub, "Vd1", "NKG2A_Freq", "% NKG2A+")
mr1 <- pbmc_v_bal(df_bgsub, "Vd1", "MR1_Freq", "% MR1-5-OPRU+")

ggarrange(cd69, hladr, cd154, cd107a, ccl4, cxcr3, ifng, tnf, il2, il10, il17, grb, grk, prf, nkg2a, mr1, nrow = 4, ncol = 4)
ggsave("out/vd1_pbmc_v_bal.pdf", units = c("in"), dpi = 300, width = 16, height = 12)

####################################################
## Does IV-BC affect the Frequency of GD T cells? ##
####################################################

V9neg_ns <- V9neg %>% filter(Stim == "NS")
V9neg_ns <- V9neg_ns %>% arrange(Timepoint)

ggpaired(V9neg_ns, x = "Timepoint", y = "Freq", id = "PTID", 
         line.color = "black",
         fill = "Timepoint",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Tissue",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  xlab("Timepoint") + ylab("Percent of Parent") +
  
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% TCRgd+/Vg9-neg of Total T cell")



#############################################################################
## How does IV-BCG alter the ex-vivo phenotype of Vg9-negative GD T cells? ##
#############################################################################

# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the trajectory of all PTIDs over time as well as the median trajectory.
# p-values represent paired comparisons between wk0/wk2, wk0/wk4, and wk0/wk8

make_pbmc_plt <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use)
  new_df2 <- new_df %>% filter(median != "yes")
  out <-
    ggpaired(
      new_df,
      x = "Timepoint",
      y = stat_use,
      id = "PTID",
      color = "median",
      line.color = "median",
      palette = c("grey", "black"),
      point.size = 1.5,
      line.size = 1.5,
      short.panel.labs = FALSE
    ) +
    ggtitle(title_use) +
    ylab("Freq of Parent") +
    xlab("Timepoint") +
    rremove("legend") +
    stat_compare_means(
      data = new_df2,
      comparisons = list(c("wk0", "wk2"), c("wk0", "wk4"), c("wk0", "wk8")),
      method = "wilcox.test",
      paired = T,
      show.legend = T
    ) +
    stat_summary(
      aes(group = timepoint),
      fun = median,
      geom = 'line',
      color = 'black'
    )
  rremove("legend")
  
  out$layers <- out$layers[-1]
  return(out)
}


V9neg_ns %>% ggpaired(x = "Timepoint", y = "CD154_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD154+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)





V9neg %>% filter(Tissue == "PBMCs") %>% ggpaired(x = "Timepoint", y = "Vd1_CD69_Freq", id = "PTID", 
                                                 line.color = "black",
                                                 fill = "Timepoint",
                                                 line.size = 0.4, 
                                                 point.size = 2,
                                                 facet.by = "Tissue",
                                                 short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD69+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_Naive_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Naive of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_CD28_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD28+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_CM_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CM of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_EM_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% EM of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_TEMRA_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% TEMRA of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "CD4_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% CD4 of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "CD8_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% CD8 of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "CD4CD8dn_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% CD4-/CD8- of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_GranB_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Granzyme B+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_GranK_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Granzyme K+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_Perforin_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Perforin+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_CD107a_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD107a+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_MR1_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% MR1+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_CXCR3_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CXCR3+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg %>% ggpaired(x = "Timepoint", y = "Vd1_NKG2A_Freq", id = "PTID", 
                   line.color = "black",
                   fill = "Timepoint",
                   line.size = 0.4, 
                   point.size = 2,
                   facet.by = "Tissue",
                   short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% NKG2A+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


#############################################
## Background Subtracted Responses to Stim ##
#############################################

# Generate the background-subtracted df
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% subset(select = -c(Sample.))

V9neg_NS <- V9neg %>% filter(Stim == "NS")
V9neg_Mtb <- V9neg %>% filter(Stim == "MtbL")

metadata <- V9neg_Mtb[,1:4]

V9neg_bgsub <- V9neg_Mtb[,5:39]-V9neg_NS[,5:39]
V9neg_bgsub <- cbind(metadata, V9neg_bgsub)

V9neg_bgsub <- V9neg_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint)
V9neg_bgsub[V9neg_bgsub < 0] <- 0


# Generate the plots

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_CD69_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD69+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_CD154_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD154+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_IFNg_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% IFNg+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_TNF_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% TNF+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_CXCR3_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CXCR3+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_CCL4_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CCL4+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_CD107a_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD107a+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_GranB_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% GranB+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_GranK_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% GranK+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_Perforin_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Perforin+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_IL17A_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% IL17A+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "Vd1_NKG2A_Freq", id = "PTID", 
                         line.color = "black",
                         fill = "Timepoint",
                         line.size = 0.4, 
                         point.size = 2,
                         facet.by = "Tissue",
                         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% NKG2A+ of Vg9-")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

###################################
## Responsiveness in BAL vs PBMC ##
###################################

# Generate the background-subtracted df
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% subset(select = -c(Sample.))

V9neg_NS <- V9neg %>% filter(Stim == "NS")
V9neg_Mtb <- V9neg %>% filter(Stim == "MtbL")

metadata <- V9neg_Mtb[,1:4]

V9neg_bgsub <- V9neg_Mtb[,5:39]-V9neg_NS[,5:39]
V9neg_bgsub <- cbind(metadata, V9neg_bgsub)

V9neg_bgsub <- V9neg_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint)
V9neg_bgsub[V9neg_bgsub < 0] <- 0


# Comparisons
V9neg_bgsub <- V9neg_bgsub %>% filter(!(Timepoint == "wk0")) %>% arrange(desc(Tissue))

ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_IFNg_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) +
  ggtitle("IFN-g")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_GranB_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) + 
  ggtitle("Granzyme B")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_TNF_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) + 
  ggtitle("TNF-a")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_CD154_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) + 
  ggtitle("CD154")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_IL2_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) +  
  ggtitle("IL-2")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(V9neg_bgsub, x = "Tissue", y = "Vd1_IL17A_Freq", id = "PTID", 
         line.color = "black",
         fill = "Tissue",
         line.size = 0.4, 
         point.size = 2,
         facet.by = "Timepoint",
         short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  stat_compare_means(comparisons = list(c("BAL", "PBMCs")), method = "wilcox.test", paired = TRUE) +
  ggtitle("IL-17a")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)



####################################################################
## HOW ENRICHED ARE MTB-RESPONSIVE CELLS IN BAL COMPARED TO PBMC? ##
####################################################################


# Generate the background-subtracted df
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% subset(select = -c(Sample.))

V9neg_NS <- V9neg %>% filter(Stim == "NS")
V9neg_Mtb <- V9neg %>% filter(Stim == "MtbL")

metadata <- V9neg_Mtb[,1:4]

V9neg_bgsub <- V9neg_Mtb[,5:40]-V9neg_NS[,5:40]
V9neg_bgsub <- cbind(metadata, V9neg_bgsub)

V9neg_bgsub <- V9neg_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint)
V9neg_bgsub[V9neg_bgsub < 0] <- 0


# Comparisons
V9neg_bgsub <- V9neg_bgsub %>% filter(!(Timepoint == "wk0")) %>% arrange(desc(Tissue))

pbmc <- V9neg_bgsub %>% filter(Tissue == "PBMCs")
bal <- V9neg_bgsub %>% filter(Tissue == "BAL")

ratio <- bal[,5:40]/pbmc[,5:40]
ratio[sapply(ratio, is.infinite)] <- NA
ratio[sapply(ratio, is.nan)] <- NA
ratio <- ratio %>% subset(select = c(Vd1_IFNg_Freq, Vd1_TNF_Freq, Vd1_CD154_Freq, Vd1_IL2_Freq))

avg <- colMeans(ratio, na.rm = TRUE)
mean(avg)




##################################
## MR1-POSITIVE CELL PHENOTYPES ##
##################################

## Unstim PBMC
MR1 <- read.csv("MR1vsVg9.csv")
MR1 <- MR1 %>% filter(Stim == "NS", Tissue == "PBMCs") %>% arrange(desc(Subset)) %>% filter(!(PTID == "A14V044" | PTID == "A14V075"))

ggpaired(MR1, x = "Subset", y = "CD154_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD154+, PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "GranB_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% GranB+, PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "GranK_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% GranK+, PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "Perforin_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Perforin+, PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD28_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD28 PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD8a_CD8aa_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD8AA PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD4_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD4 PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(MR1, x = "Subset", y = "Naive_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Naive, PBMC")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


## Unstim BAL

MR1 <- read.csv("MR1vsVg9.csv")
MR1 <- MR1 %>% filter(Stim == "NS", Tissue == "BAL") %>% arrange(desc(Subset))

ggpaired(MR1, x = "Subset", y = "CD154_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD154 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "GranB_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("GranB BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "GranK_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("GranK BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "Perforin_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("Perforin BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD28_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD28 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CCL4_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CCL4 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD8a_CD8aa_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD8AA BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CD4_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CD4 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CXCL10_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CXCL10 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1, x = "Subset", y = "CXCR3_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Timepoint", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("CXCR3 BAL")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


#####################################################################
## Phenotypes in MR1 vs Total Vg9-    ALL TISSUES AND TIMEPOINTS   ##
#####################################################################

MR1 <- read.csv("MR1vsVg9.csv")
MR1 <- MR1 %>% filter(Stim == "NS") %>% arrange(desc(Subset)) %>% filter(!(PTID == "A14V044" | PTID == "A14V075"))


ggpaired(MR1, x = "Subset", y = "CD154_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CD154+")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(MR1, x = "Subset", y = "Naive_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("Vg9-", "MR1")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Naive")+
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)





#########################################
## Responsiveness in MR1 vs Total Vg9- ##
#########################################


# Generate the background-subtracted df
MR1 <- read.csv("MR1vsVg9.csv")
MR1 <- MR1 %>% subset(select = -c(Sample., Count, Freq, X))

MR1_NS <- MR1 %>% filter(Stim == "NS")
MR1_Mtb <- MR1 %>% filter(Stim == "MtbL")

metadata <- MR1_Mtb[,1:5]
MR1_bgsub <- MR1_Mtb[,6:36]-MR1_NS[,6:36]
MR1_bgsub <- cbind(metadata, MR1_bgsub)
MR1_bgsub <- MR1_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint) %>% filter(Subset == "MR1")

MR1_bgsub[MR1_bgsub < 0] <- 0

ggpaired(MR1_bgsub, x = "Timepoint", y = "IFNg_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% IFNg+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(MR1_bgsub, x = "Timepoint", y = "TNFa_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% TNFa+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(MR1_bgsub, x = "Timepoint", y = "IL6_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% IL-6+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggpaired(MR1_bgsub, x = "Timepoint", y = "IL17A_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% IL-17A+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1_bgsub, x = "Timepoint", y = "CCL4_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CCL4+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1_bgsub, x = "Timepoint", y = "CXCL10_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% CXCL10+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggpaired(MR1_bgsub, x = "Timepoint", y = "Perforin_Freq", id = "PTID", 
         line.color = "gray", 
         line.size = 0.4, 
         facet.by = "Tissue", 
         short.panel.labs = FALSE) +
  stat_compare_means(comparisons = list(c("wk0", "wk4"), c("wk0", "wk8"), c("wk4", "wk8")), method = "wilcox.test", paired = TRUE) +
  ggtitle("% Perforin+") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)




###########################################
## Figures for TB Symposium Presentation ##
###########################################

# Slide 7
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% filter(Stim == "NS", Tissue == "PBMCs") %>% arrange(Timepoint)

V9neg <- V9neg %>% filter(Stim == "NS", Tissue == "BAL") %>% arrange(Timepoint)

a <- V9neg %>% ggpaired(x = "Timepoint", y = "GranB_Freq", id = "PTID", 
                        line.color = "black",
                        fill = "Timepoint",
                        line.size = 0.4, 
                        point.size = 2,
                        short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("Granzyme B")+
  ylab("Frequency of Total Vg9- T cells") +
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

b <- V9neg %>% ggpaired(x = "Timepoint", y = "GranK_Freq", id = "PTID", 
                        line.color = "black",
                        fill = "Timepoint",
                        line.size = 0.4, 
                        point.size = 2,
                        short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("Granzyme K")+
  ylab("Frequency of Total Vg9- T cells") +
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

c <- V9neg %>% ggpaired(x = "Timepoint", y = "Perforin_Freq", id = "PTID", 
                        line.color = "black",
                        fill = "Timepoint",
                        line.size = 0.4, 
                        point.size = 2,
                        short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("Perforin")+
  ylab("Frequency of Total Vg9- T cells") +
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggarrange(a,
          b +  theme(axis.title.y = element_blank()),
          c + theme(axis.title.y = element_blank()),
          ncol = 3, nrow = 1, common.legend = TRUE, align = "v")

ggsave("out/images/cytotoxic_markers_pbmc.png", bg = "white")



# Slide 8
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% subset(select = -c(Sample.)) %>% filter(Tissue == "PBMCs")

V9neg <- V9neg %>% subset(select = -c(Sample.)) %>% filter(Tissue == "BAL")

V9neg_NS <- V9neg %>% filter(Stim == "NS")
V9neg_Mtb <- V9neg %>% filter(Stim == "MtbL")

metadata <- V9neg_Mtb[,1:5]

V9neg_bgsub <- V9neg_Mtb[,6:36]-V9neg_NS[,6:36]
V9neg_bgsub <- cbind(metadata, V9neg_bgsub)

V9neg_bgsub <- V9neg_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint)
V9neg_bgsub[V9neg_bgsub < 0] <- 0


a <- V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "CD69_Freq", id = "PTID", 
                              line.color = "black",
                              fill = "Timepoint",
                              line.size = 0.4, 
                              point.size = 2,
                              short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("CD69")+
  ylab("Frequency of Total Vg9- T cells") + 
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

b <- V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "IFNg_Freq", id = "PTID", 
                              line.color = "black",
                              fill = "Timepoint",
                              line.size = 0.4, 
                              point.size = 2,
                              short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("IFN-gamma")+
  ylab("Frequency of Total Vg9- T cells") +
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

c <- V9neg_bgsub %>% ggpaired(x = "Timepoint", y = "TNF_Freq", id = "PTID", 
                              line.color = "black",
                              fill = "Timepoint",
                              line.size = 0.4, 
                              point.size = 2,
                              short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("TNF")+
  ylab("Frequency of Total Vg9- T cells") +
  xlab("Timepoint") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)


ggarrange(a + theme(plot.margin = unit(c(1,1,1,1), "mm")),
          b + theme(axis.title.y = element_blank(),
                    plot.margin = unit(c(1,1,1,1), "mm")),
          c + theme(axis.title.y = element_blank(),
                    plot.margin = unit(c(1,1,1,1), "mm")), 
          ncol = 3, 
          nrow = 1, 
          common.legend = TRUE,
          align = "v")

ggsave("out/images/responses.png", bg = "white")


# Slide 8 contd
V9neg <- read.csv("regate_Vd1.csv")
V9neg <- V9neg %>% subset(select = -c(Sample.))

V9neg_NS <- V9neg %>% filter(Stim == "NS")
V9neg_Mtb <- V9neg %>% filter(Stim == "MtbL")

metadata <- V9neg_Mtb[,1:5]

V9neg_bgsub <- V9neg_Mtb[,6:39]-V9neg_NS[,6:39]
V9neg_bgsub <- cbind(metadata, V9neg_bgsub)

V9neg_bgsub <- V9neg_bgsub %>% arrange(-desc(Tissue)) %>% arrange(Timepoint)
V9neg_bgsub[V9neg_bgsub < 0] <- 0

V9neg_bgsub <- V9neg_bgsub %>% arrange(desc(Tissue))


a <- ggpaired(V9neg_bgsub, x = "Tissue", y = "IFNg_Freq", id = "PTID", 
              line.color = "black",
              fill = "Tissue",
              line.size = 0.4, 
              point.size = 2,
              facet.by = "Timepoint",
              short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("IFN-gamma")+
  ylab("Frequency of Total Vd1 T cells")
geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

b <- ggpaired(V9neg_bgsub, x = "Tissue", y = "TNF_Freq", id = "PTID", 
              line.color = "black",
              fill = "Tissue",
              line.size = 0.4, 
              point.size = 2,
              facet.by = "Timepoint",
              short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("TNF")+
  ylab("Frequency of Total Vd1 T cells")
geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 1.2)

ggarrange(a,
          b +  theme(axis.title.y = element_blank()),
          ncol = 2, nrow = 1, common.legend = TRUE, align = "v")


