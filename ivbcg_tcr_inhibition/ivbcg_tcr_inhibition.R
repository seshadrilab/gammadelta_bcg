
# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggedit)

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

# Read in the CSV
df <- read.csv("tcr_inhibition.csv")

## NOTE: In the dataframe, US = unstimulated,Cytok = cytokine-stimulated, CsA = Cyclosporin A, NI = no inhibitor
## The experiments also included a MAPK inhibitor, but we did not include these data due to broad effects of MAPK inhibition

####################################
## Background subtraction -- PBMC ##
####################################

# Filter out the AB T cells and only include PBMC
df <- df %>% filter(subset != 'tcrab' & tissue == 'pbmc')

# Set the levels for stimulation and cell subset (US = Unstimulated)
df$stim <- factor(df$stim, levels = c('us', 'pha', 'cytok', 'mtbl'))
df$subset <- factor(df$subset, levels = c('vg9neg', 'vg9pos', 'cd8', 'cd4'))

# Subset out the stimulation conditions
df <- df %>% arrange(tissue, timepoint, stim, inhib)
df_us <- df %>% filter(stim == 'us', tissue == 'pbmc')
df_pha <- df %>% filter(stim == 'pha', tissue == 'pbmc')
df_cyt <- df %>% filter(stim == 'cytok', tissue == 'pbmc')
df_mtb <- df %>% filter(stim == 'mtbl', tissue == 'pbmc')

# For each stimulation condition, extract the metadata and merge with 
# background-subtracted measurements
metadata <- df_pha[,1:15]
pha_bgsub <- df_pha[,16:27]-df_us[,16:27]
pha_bgsub <- cbind(metadata, pha_bgsub)
pha_bgsub[pha_bgsub < 0] <- 0

metadata <- df_cyt[,1:15]
cyt_bgsub <- df_cyt[,16:27]-df_us[,16:27]
cyt_bgsub <- cbind(metadata, cyt_bgsub)
cyt_bgsub[cyt_bgsub < 0] <- 0

metadata <- df_mtb[,1:15]
mtb_bgsub <- df_mtb[,16:27]-df_us[,16:27]
mtb_bgsub <- cbind(metadata, mtb_bgsub)
mtb_bgsub[mtb_bgsub < 0] <- 0


#######################
## Calculate medians ##
#######################

## MtbL-stimulated
# Calculate median values of readouts from each condition group
medians <- mtb_bgsub %>% 
  group_by(tissue, timepoint, stim, inhib, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:tnf), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians with df
mtb_bgsub <- dplyr::bind_rows(mtb_bgsub, medians) #combine the tables


## PHA-stimulated
# Calculate median values of readouts from each condition group
medians <- pha_bgsub %>% 
  group_by(tissue, timepoint, stim, inhib, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:tnf), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians with df
pha_bgsub <- dplyr::bind_rows(pha_bgsub, medians) #combine the tables



## Cytokine-stimulated
# Calculate median values of readouts from each condition group
medians <- cyt_bgsub %>% 
  group_by(tissue, timepoint, stim, inhib, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:tnf), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians with df
cyt_bgsub <- dplyr::bind_rows(cyt_bgsub, medians) #combine the tables



#################################
## Are the inhibitors working? ##
#################################

## Repeated code blocks visualize CD137 and IFN-g

pha_bgsub$inhib <- factor(pha_bgsub$inhib, levels = c('ni', 'csa', 'mapk'))
cyt_bgsub$inhib <- factor(cyt_bgsub$inhib, levels = c('ni', 'csa', 'mapk'))
mtb_bgsub$inhib <- factor(mtb_bgsub$inhib, levels = c('ni', 'csa', 'mapk'))


## PHA plus inhibitors
pha_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 8 pbmc + PHA")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none')

pha_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 8 pbmc + PHA")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none')


## Cytokine plus inhibitors
cyt_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 8 pbmc + cytok")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none')

cyt_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 8 pbmc + cytok")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none')



#########################################################
## How do the inhibitors affect MtbL-stimulated cells? ##
#########################################################

## Set up the function
make_plot <- function(in_df, subset_use, stat_use, title_use) {
  new_df <- in_df %>% filter(subset == subset_use)
  new_df2 <- new_df %>% filter(median != "yes")
  out <-
    ggpaired(
      new_df,
      x = "inhib",
      y = stat_use,
      id = "ptid",
      color = "median",
      line.color = "median",
      facet.by = 'timepoint',
      nrow = 1,
      ncol = 4,
      palette = c("black", "grey"),
      point.size = 1.5,
      line.size = 1.5,
      short.panel.labs = FALSE
    ) +
    ggtitle(title_use) +
    ylab("Freq of Parent") +
    xlab("Inhibitor") +
    rremove("legend") +
    stat_compare_means(
      data = new_df2,
      comparisons = list(c("ni", "csa")),
      method = "wilcox.test",
      paired = T,
      show.legend = T
    ) +
    stat_summary(
      aes(group = inhib),
      fun = median,
      geom = 'line',
      color = 'black'
    )
  rremove("legend")
  
  out$layers <- out$layers[-1]
  return(out)
}


## VG9-NEG ##
a <-  mtb_bgsub %>% filter(inhib != 'mapk') %>% make_plot("vg9neg", "cd137", "% CD137+")
b <- mtb_bgsub %>% filter(inhib != 'mapk') %>% make_plot("vg9neg", "ifng", "% IFN-γ+")
ggarrange(a,b, nrow = 1, ncol = 2)
ggsave('Fig5G.pdf', width = 10, height = 3, units = 'in')



#################################################################
## How does CsA affect responses in cytokine-stimulated cells? ##
#################################################################

#Vg9neg
a <-  cyt_bgsub %>% filter(inhib != 'mapk') %>% make_plot("vg9neg", "cd137", "% CD137+")
a <- a + expand_limits(y = 0)
b <- cyt_bgsub %>% filter(inhib != 'mapk') %>% make_plot("vg9neg", "ifng", "% IFN-γ+")
b <- b + expand_limits(y = 0)
ggarrange(a,b, nrow = 1, ncol = 2)
ggsave(filename = 'SuppFig14.pdf', width = 10, height = 3, units = 'in')



####################################################################
## Are the cells more responsive to PHA after IV-BCG? ##
####################################################################

## PHA
pha_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc', subset == 'vg9neg') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, pbmc + pha")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none')
ggsave(filename= 'SuppFig13A.pdf')

pha_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc', subset == 'vg9neg') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, pbmc + pha")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') 
ggsave(filename= 'SuppFig13B.pdf')
