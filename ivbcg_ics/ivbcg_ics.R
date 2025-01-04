
# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(readxl)


####################
## SET UP TABLES ###
####################

# Read in the main table
raw <- read.csv("concat_regate.csv")
raw <- raw %>% filter(ptid != '36852') # Exclude ptid #36852 which is missing wk2 data
raw$median <- "no"
raw <- raw %>% select('median', everything()) %>%
  arrange(timepoint)

# Get names of unvax BAL donors
temp <- raw %>% filter(tissue == 'BAL' & timepoint == 'pre')
temp$ptid %>% unique()

# Calculate median values of readouts from each condition group
medians <- raw %>% 
  arrange(timepoint) %>%
  group_by(tissue, timepoint, stim, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:cytokine), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians and SDs with df
df <- dplyr::bind_rows(raw, medians)

# Filter the dataframe by unstimulated pbmc
us_p <- df %>% filter(stim == "unstim" & tissue == "PBMC") %>%
  arrange(timepoint)

# Filter the dataframe by unstimulated bal
us_b <- df %>% filter(stim == "unstim" & tissue == "BAL") %>%
  arrange(timepoint)

# Calculate the background-subtracted responses to stimulation
df_NS <- raw %>% filter(stim == "unstim")
df_Mtb <- raw %>% filter(stim == "mtbl")
metadata <- df_Mtb[,1:14]
df_bgsub <- df_Mtb[,15:29]-df_NS[,15:29]
df_bgsub <- cbind(metadata, df_bgsub)
df_bgsub <- df_bgsub %>% arrange(-desc(tissue)) %>% arrange(timepoint)
df_bgsub[df_bgsub < 0] <- 0

# Calculate median of readouts from each condition group
medians <- df_bgsub %>% 
  group_by(tissue, timepoint, stim, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:cytokine), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians with df
df_bgsub <- dplyr::bind_rows(df_bgsub, medians) #combine the tables

# Export to use in prism
write.csv(df_bgsub, "bgsub.csv")

# Save new tables containing only PBMC or BAL
df_bgs_p <- df_bgsub %>% filter(tissue == 'PBMC')
df_bgs_b <- df_bgsub %>% filter(tissue == "BAL")

# Generate a table containing only Week 0 and Week 8 for PBMC vs. BAL comparisons
df_pvb <- df %>% filter(stim == 'unstim' & (timepoint == 'pre' | timepoint == 'unvax' | timepoint == 'wk.8'))

# Generate a table containing background-subtracted responses at wk0 and wk8
df_bgsub_pvb <- df_bgsub %>% filter(timepoint == 'pre' | timepoint == 'unvax' | timepoint == 'wk.8') %>% 
  arrange(desc(tissue))

# Export tables
write.csv(us_p, 'out/us_p.csv')
write.csv(us_b, 'out/us_b.csv')
write.csv(df_bgs_p, 'out/df_bgs_p.csv')
write.csv(df_bgs_b, 'out/df_bgs_b.csv')
write.csv(df_pvb, 'out/df_pvb.csv')
write.csv(df_bgsub_pvb, 'out/df_bgsub_pvb.csv')

#######################
## SET UP FUNCTIONS ###
#######################

# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the trajectory of all PTIDs over time as well as the median trajectory.
# p-values represent paired comparisons between wk0/wk2, wk0/wk4, and wk0/wk8
make_pbmc_plt <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use)
  new_df2 <- new_df %>% filter(median != "yes")
  out <-
    ggpaired(
      new_df,
      x = "timepoint",
      y = stat_use,
      id = "ptid",
      color = "median",
      line.color = "median",
      palette = c("grey", "black"),
      point.size = 'median',
      line.size = "median",
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Timepoint") +
    rremove("legend") +
    stat_compare_means(
      data = new_df2,
      comparisons = list(c("pre", "wk.2"), c("pre", "wk.4"), c("pre", "wk.8")),
      method = "wilcox.test",
      paired = T,
      show.legend = T,
      size = 5,
      label = 'p.format'
    )
  rremove("legend")
  
  out$layers <- out$layers[-1]
  return(out)
}


# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the trajectory of all PTIDs over time as well as the median trajectory.
# p-values represent unpaired comparisons between wk0/wk8
make_bal_plt <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use & median != 'yes')
  out <-
    ggpaired(
      new_df,
      x = "timepoint",
      y = stat_use,
      id = "ptid",
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Timepoint") +
    rremove("legend") +
    stat_compare_means(
      comparisons = list(c("unvax", "wk.8")),
      method = "wilcox.test",
      paired = F,
      size = 5,
      show.legend = T,
      label = 'p.format'
    ) 
  rremove("legend")
  
  return(out)
}

# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the difference between PBMC and BAL at Week 0.
# p-values represent unpaired comparisons between PBMC and BAL
pbmc_v_bal_0 <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use & median != 'yes')
  out <-
    ggpaired(
      new_df,
      x = "tissue",
      y = stat_use,
      id = "ptid",
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Tissue") +
    rremove("legend") +
    stat_compare_means(
      comparisons = list(c("PBMC", "BAL")),
      method = "wilcox.test",
      paired = F,
      size = 5,
      show.legend = T,
      label = 'p.format'
    ) 
  rremove("legend")
  
  return(out)
}


# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the difference between PBMC and BAL at Week 8.
# p-values represent unpaired comparisons between PBMC and BAL
pbmc_v_bal_8 <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use)
  new_df2 <- new_df %>% filter(median != "yes")
  out <-
    ggpaired(
      new_df,
      x = "tissue",
      y = stat_use,
      id = "ptid",
      color = "median",
      line.color = "median",
      palette = c("grey", "black"),
      point.size = 'median',
      line.size = 'median',
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Tissue") +
    rremove("legend") +
    stat_compare_means(
      data = new_df2,
      comparisons = list(c("PBMC", "BAL")),
      method = "wilcox.test",
      paired = T,
      size = 5,
      show.legend = T,
      label = 'p.format'
    ) 
  rremove("legend")
  
  out$layers <- out$layers[-1]
  return(out)
}


# Create a function generating a ggplot which shows a given statistic across all T cell subsets.
# The plot displays the kinetics over time for a given statistic in PBMC.
# The plots are faceted to keep the y-axis constant across subsets.
# p-values represent paired comparisons between timepoints.
comp <- function(in_df, stat_use, title_use) {
  in_df <- in_df %>% filter(subset != "tcrab")
  new_df <- in_df %>% filter(median != "yes")
  out <-
    ggpaired(
      in_df,
      x = "timepoint",
      y = stat_use,
      id = "ptid",
      facet.by = "subset",
      color = "median",
      line.color = "median",
      palette = c("grey", "black"),
      point.size = 'median',
      line.size = 'median',
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Timepoint") +
    rremove("legend") +
    stat_compare_means(
      data = new_df,
      comparisons = list(c("pre", "wk.2"), c("pre", "wk.4", c("pre", "wk.8"))),
      method = "wilcox.test",
      paired = T,
      size = 5,
      show.legend = T,
      label = 'p.format'
    )
  
  out$layers <- out$layers[-1]
  return(out)
}


# Create a function generating a ggplot which shows a given statistic across all T cell subsets.
# The plot displays the kinetics over time for a given statistic in BAL.
# The plots are faceted to keep the y-axis constant across subsets.
# p-values represent paired comparisons between week 0 and week 8.
comp_b <- function(in_df, stat_use, title_use) {
  in_df <- in_df %>% filter(subset != "tcrab")
  new_df <- in_df %>% filter(median != "yes")
  out <-
    ggpaired(
      in_df,
      x = "timepoint",
      y = stat_use,
      id = "ptid",
      facet.by = "subset",
      color = "median",
      line.color = "median",
      palette = c("grey", "black"),
      point.size = 'median',
      line.size = 'median',
      short.panel.labs = FALSE
    ) +
    scale_size_manual(values = c(0.5, 1.5)) +
    
    ggtitle(title_use) +
    ylab("Freq of Vd1 T cell") +
    xlab("Timepoint") +
    rremove("legend") +
    stat_compare_means(
      data = new_df,
      comparisons = list(c("pre", "wk.8")),
      method = "wilcox.test",
      paired = F,
      size = 5,
      show.legend = T,
      label = 'p.format'
    )
  
  
  out$layers <- out$layers[-1]
  return(out)
}


########################
## EX VIVO PHENOTYPES ##
########################

### PBMC ###
# Visualize the frequency of all T cell subsets over time
vg9 <- make_pbmc_plt(us_p, "vg9pos", "subset_pct_of_tcell", "% Vg9+ of T cell")
vd1 <- make_pbmc_plt(us_p, "vg9neg", "subset_pct_of_tcell", "% Vg9neg of T cell")
cd4 <- make_pbmc_plt(us_p, "cd4", "subset_pct_of_tcell", "% CD4 of T cell")
cd8 <- make_pbmc_plt(us_p, "cd8", "subset_pct_of_tcell", "% CD8 of T cell")
ggarrange(vg9, vd1, cd4, cd8)
ggsave("out/subset_frq.pdf", units = c("in"), dpi = 300, width = 5, height = 7.5)


# Visualize the analytes in unstimulated Vg9neg over time
cd4 <- make_pbmc_plt(us_p, "vg9neg", "cd4", "% CD4+")
cd8 <- make_pbmc_plt(us_p, "vg9neg", "cd8", "% CD8+")
dn <- make_pbmc_plt(us_p, "vg9neg", "dn", "% DN+")
activ <- make_pbmc_plt(us_p, "vg9neg", "activated", "% CD69+ or HLA-DR+")
cd137 <- make_pbmc_plt(us_p, "vg9neg", "cd137", "% CD137+")
grz <- make_pbmc_plt(us_p, "vg9neg", "granzyme", "% Granzyme+")
cd107a <- make_pbmc_plt(us_p, "vg9neg", "cd107a", "% CD107a+")
cytok <- make_pbmc_plt(us_p, "vg9neg", "cytokine", "% Cytokine+")
ggarrange(cd4, cd8, dn, activ, cd137, grz, cd107a, cytok)
ggsave("out/vd1_pheno.pdf", units = c("in"), dpi = 300, width = 12, height = 12)


ggarrange(cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok + rremove ('xlab'), 
          ncol = 1, nrow = 4)
ggsave("out/vd1_unstim_pbmc.pdf", units = c('in'), dpi = 300, width = 4, height = 19)



### BAL ###
# Visualize the frequency of all T cell subsets over time
vg9 <- make_bal_plt(us_b, "vg9pos", "subset_pct_of_tcell", "% Vg9+ of T cell")
vd1 <- make_bal_plt(us_b, "vg9neg", "subset_pct_of_tcell", "% Vg9neg of T cell")
cd4 <- make_bal_plt(us_b, "cd4", "subset_pct_of_tcell", "% CD4 of T cell")
cd8 <- make_bal_plt(us_b, "cd8", "subset_pct_of_tcell", "% CD8 of T cell")
ggarrange(vg9, vd1, cd4, cd8)
ggsave("out/subset_frq_bal.pdf", units = c("in"), dpi = 300, width = 5, height = 7.5)

# Visualize the analytes in unstim Vg9neg over time
us_b$timepoint <- as.character(us_b$timepoint)
us_b <- us_b %>% arrange(timepoint)

cd4 <- make_bal_plt(us_b, "vg9neg", "cd4", "% CD4+")
cd8 <- make_bal_plt(us_b, "vg9neg", "cd8", "% CD8+")
dn <- make_bal_plt(us_b, "vg9neg", "dn", "% DN+")
hladr <- make_bal_plt(us_b, "vg9neg", "hladr", "% HLA-DR+")
cd137 <- make_bal_plt(us_b, "vg9neg", "cd137", "% CD137+")
grz <- make_bal_plt(us_b, "vg9neg", "granzyme", "% Granzyme+")
cd107a <- make_bal_plt(us_b, "vg9neg", "cd107a", "% CD107a+")
cytok <- make_bal_plt(us_b, "vg9neg", "cytokine", "% Cytokine+")
ggarrange(cd4, cd8, dn, activ, cd137, grz, cd107a, cytok)
ggsave("out/vd1_pheno_bal.pdf", units = c("in"), dpi = 300, width = 12, height = 12)

ggarrange(cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok, 
          ncol = 1, nrow = 4)
ggsave("out/vd1_unstim_bal.pdf", units = c('in'), dpi = 300, width = 4, height = 19)



##############################
## RESPONSES TO STIMULATION ##
##############################

## PBMC ###
#Visualize responses to stimulation in Vg9neg
activ <- make_pbmc_plt(df_bgs_p, "vg9neg", "activated", "% CD69/HLA-DR+") + theme(text=element_text(size=16))
cd137 <- make_pbmc_plt(df_bgs_p, "vg9neg", "cd137", "% CD137+") + theme(text=element_text(size=16))
grz <- make_pbmc_plt(df_bgs_p, "vg9neg", "granzyme", "% Granzyme+") + theme(text=element_text(size=16))
cd107a <- make_pbmc_plt(df_bgs_p, "vg9neg", "cd107a", "% CD107a+") + theme(text=element_text(size=16))
cytok <- make_pbmc_plt(df_bgs_p, "vg9neg", "cytokine", "% Cytokine+") + theme(text=element_text(size=16))

ggarrange(activ + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok, 
          ncol = 1, nrow = 5)
ggsave("out/vd1_stim_pformat.pdf", units = c('in'), dpi = 300, width = 4, height = 19)


## BAL ##
#Visualize responses to stimulation in Vg9neg
hladr <- make_bal_plt(df_bgs_b, "vg9neg", "hladr", "% HLADR+") + theme(text=element_text(size=16))
cd137 <- make_bal_plt(df_bgs_b, "vg9neg", "cd137", "% CD137+") + theme(text=element_text(size=16))
grz <- make_bal_plt(df_bgs_b, "vg9neg", "granzyme", "% Granzyme+") + theme(text=element_text(size=16))
cd107a <- make_bal_plt(df_bgs_b, "vg9neg", "cd107a", "% CD107a+") + theme(text=element_text(size=16))
cytok <- make_bal_plt(df_bgs_b, "vg9neg", "cytokine", "% Cytokine+") + theme(text=element_text(size=16))

ggarrange(hladr + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok, 
          ncol = 1, nrow = 5)
ggsave("out/vd1_stim_bal_pformat.pdf", units = c('in'), dpi = 300, width = 4, height = 19)



###################################
## PBMC vs BAL RESPONSES TO MTBL ##
###################################

df_bgsub_pvb_0 <- df_bgsub_pvb %>% filter(timepoint == 'pre' | timepoint == 'unvax')
df_bgsub_pvb_8 <- df_bgsub_pvb %>% filter(timepoint == 'wk.8') %>%
  filter(is.na(ptid) | ptid != 'D12L') %>%
  filter(is.na(ptid) | ptid != 'DGHL')


# Week 0
hladr <- pbmc_v_bal_0(df_bgsub_pvb_0, "vg9neg", "hladr", "% HLADR+")
cd137 <- pbmc_v_bal_0(df_bgsub_pvb_0, "vg9neg", "cd137", "% CD137+")
grz <- pbmc_v_bal_0(df_bgsub_pvb_0, "vg9neg", "granzyme", "% Granzyme+")
cd107a <- pbmc_v_bal_0(df_bgsub_pvb_0, "vg9neg", "cd107a", "% CD107a+")
cytok <- pbmc_v_bal_0(df_bgsub_pvb_0, "vg9neg", "cytokine", "% Cytokine+")

ggarrange(hladr + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok + rremove ('xlab'), 
          ncol = 5, nrow = 1)
ggsave("out/vd1_pvb_bgsub_wk0_pformat.pdf", units = c('in'), dpi = 300, width = 18, height = 4)


# Week 8
hladr <- pbmc_v_bal_8(df_bgsub_pvb_8, "vg9neg", "hladr", "% HLADR+") + theme(text=element_text(size=16))
cd137 <- pbmc_v_bal_8(df_bgsub_pvb_8, "vg9neg", "cd137", "% CD137+") + theme(text=element_text(size=16))
grz <- pbmc_v_bal_8(df_bgsub_pvb_8, "vg9neg", "granzyme", "% Granzyme+") + theme(text=element_text(size=16))
cd107a <- pbmc_v_bal_8(df_bgsub_pvb_8, "vg9neg", "cd107a", "% CD107a+") + theme(text=element_text(size=16))
cytok <- pbmc_v_bal_8(df_bgsub_pvb_8, "vg9neg", "cytokine", "% Cytokine+") + theme(text=element_text(size=16))

ggarrange(hladr + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok + rremove('xlab'), 
          ncol = 5, nrow = 1)
ggsave("out/vd1_pvb_bgsub_wk8_pformat.pdf", units = c('in'), dpi = 300, width = 18, height = 4)



##############################
## PBMC vs BAL UNSTIMULATED ##
##############################

df_pvb_0 <- df_pvb %>% filter(timepoint == 'pre' | timepoint == 'unvax')
df_pvb_8 <- df_pvb %>% filter(timepoint == 'wk.8') %>%
  filter(is.na(ptid) | ptid != 'D12L') %>%
  filter(is.na(ptid) | ptid != 'DGHL') %>%
  arrange(desc(tissue))
# These two ptids are removed because they each only have 
# PBMC or BAL samples, so PBMC vs. BAL is not possible


# Week 0
hladr <- pbmc_v_bal_0(df_pvb_0, "vg9neg", "hladr", "% HLADR+")
cd137 <- pbmc_v_bal_0(df_pvb_0, "vg9neg", "cd137", "% CD137+")
grz <- pbmc_v_bal_0(df_pvb_0, "vg9neg", "granzyme", "% Granzyme+")
cd107a <- pbmc_v_bal_0(df_pvb_0, "vg9neg", "cd107a", "% CD107a+")
cytok <- pbmc_v_bal_0(df_pvb_0, "vg9neg", "cytokine", "% Cytokine+")

ggarrange(hladr + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok + rremove ('xlab'), 
          ncol = 5, nrow = 1)
ggsave("out/vd1_pvb_wk0_pformat.pdf", units = c('in'), dpi = 300, width = 18, height = 4)


# Week 8
hladr <- pbmc_v_bal_8(df_pvb_8, "vg9neg", "hladr", "% HLADR+") + theme(text=element_text(size=16))
cd137 <- pbmc_v_bal_8(df_pvb_8, "vg9neg", "cd137", "% CD137+") + theme(text=element_text(size=16))
grz <- pbmc_v_bal_8(df_pvb_8, "vg9neg", "granzyme", "% Granzyme+") + theme(text=element_text(size=16))
cd107a <- pbmc_v_bal_8(df_pvb_8, "vg9neg", "cd107a", "% CD107a+") + theme(text=element_text(size=16))
cytok <- pbmc_v_bal_8(df_pvb_8, "vg9neg", "cytokine", "% Cytokine+") + theme(text=element_text(size=16))

ggarrange(hladr + rremove ('xlab'), 
          cd137 + rremove ('xlab'), 
          grz + rremove ('xlab'), 
          cd107a + rremove ('xlab'), 
          cytok + rremove('xlab'), 
          ncol = 5, nrow = 1)
ggsave("out/vd1_pvb_wk8_pformat.pdf", units = c('in'), dpi = 300, width = 18, height = 4)




#######################################
## COMPARE SUBSET PHENOTYPES IN PBMC ##
#######################################

# Visualize responses to stimulation by subset, this time all in one figure
hladr <- comp(us_p, "hladr", "% HLA-DR+")
cd69 <- comp(us_p, "cd69", "% CD69+")
cd137 <- comp(us_p, "cd137", "% CD137+")
grb <- comp(us_p, "granzymeb", "% Granzyme B+")
grk <- comp(us_p, "granzymek", "% Granzyme K+")
cd107a <- comp(us_p, "cd107a", "% CD107a+")
ifng <- comp(us_p, "ifng", "% IFN-g+")
tnf <- comp(us_p, "tnf", "% TNF+")
ggarrange(hladr, cd69, cd137, grb, grk, cd107a, ifng, tnf, ncol = 4, nrow = 2)
ggsave("out/compare_all_pheno.pdf", units = c("in"), dpi = 300, width = 16, height = 17)



#######################################
## COMPARE SUBSET PHENOTYPES IN BAL ##
#######################################

# Visualize responses to stimulation by subset, this time all in one figure
hladr <- comp_b(us_b, "hladr", "% HLA-DR+")
cd69 <- comp_b(us_b, "cd69", "% CD69+")
cd137 <- comp_b(us_b, "cd137", "% CD137+")
grb <- comp_b(us_b, "granzymeb", "% Granzyme B+")
grk <- comp_b(us_b, "granzymek", "% Granzyme K+")
cd107a <- comp_b(us_b, "cd107a", "% CD107a+")
ifng <- comp_b(us_b, "ifng", "% IFN-g+")
tnf <- comp_b(us_b, "tnf", "% TNF+")
ggarrange(hladr, cd69, cd137, grb, grk, cd107a, ifng, tnf, ncol = 4, nrow = 2)
ggsave("out/compare_all_pheno_bal.pdf", units = c("in"), dpi = 300, width = 16, height = 17)



###################################
## COMPARE SUBSETS -- LINE PLOTS ##
###################################

my_cols3 <- viridis(3, end = 0.75)
df_bgs_p$subset <- factor(df_bgs_p$subset, levels = c("vg9neg", "vg9pos", "cd8"))
df_bgs_b$subset <- factor(df_bgs_b$subset, levels = c("vg9neg", "vg9pos", "cd8"))

## PBMC over time
cd137 <- df_bgs_p %>% 
  filter(median == 'yes' & subset != 'tcrab' & subset != 'cd4') %>% 
  ggplot(aes(x = timepoint, y = cd137, group = subset, color = subset)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = my_cols3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0,0.3)))

grz <- df_bgs_p %>% 
  filter(median == 'yes' & subset != 'tcrab' & subset != 'cd4') %>% 
  ggplot(aes(x = timepoint, y = granzyme, group = subset, color = subset)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = my_cols3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0,0.3)))

cd107a <- df_bgs_p %>% 
  filter(median == 'yes' & subset != 'tcrab' & subset != 'cd4') %>% 
  ggplot(aes(x = timepoint, y = cd107a, group = subset, color = subset)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = my_cols3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0,0.3)))

cytok <- df_bgs_p %>% 
  filter(median == 'yes' & subset != 'tcrab' & subset != 'cd4') %>% 
  ggplot(aes(x = timepoint, y = cytokine, group = subset, color = subset)) +
  geom_line(linewidth = 1.4) +
  scale_color_manual(values = my_cols3) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0,0.3)))

ggarrange(cd137, grz, cd107a, cytok, nrow = 4, ncol = 1, common.legend = T)

ggsave("out/comp_subset_p_line_pformat.pdf", units = c("in"), dpi = 300, width = 3.5, height = 12)

 
