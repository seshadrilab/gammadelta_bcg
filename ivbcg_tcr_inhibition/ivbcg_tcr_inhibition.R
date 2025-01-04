
# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggedit)

# Read in the CSVs
df <- read.csv("table_concat2.csv")


####################################
## Background subtraction -- PBMC ##
####################################

df <- df %>% filter(subset != 'tcrab' & tissue == 'pbmc')

df$stim <- factor(df$stim, levels = c('us', 'pha', 'cytok', 'mtbl'))
df$subset <- factor(df$subset, levels = c('vg9neg', 'vg9pos', 'cd8', 'cd4'))

df <- df %>% arrange(tissue, timepoint, stim, inhib)

df_us <- df %>% filter(stim == 'us', tissue == 'pbmc')
df_pha <- df %>% filter(stim == 'pha', tissue == 'pbmc')
df_cyt <- df %>% filter(stim == 'cytok', tissue == 'pbmc')
df_mtb <- df %>% filter(stim == 'mtbl', tissue == 'pbmc')

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

# Calculate median values of readouts from each condition group
medians <- mtb_bgsub %>% 
  group_by(tissue, timepoint, stim, inhib, subset) %>% 
  summarize_at(vars(subset_pct_of_tcell:tnf), median, na.rm=TRUE)
medians$median <- "yes" #mark each value as the median

# Merge medians with df
mtb_bgsub <- dplyr::bind_rows(mtb_bgsub, medians) #combine the tables



#######################
## SET UP FUNCTIONS ###
#######################

# Create a function generating a ggplot which shows a given statistic within a T cell subset.
# The plot displays the trajectory of all PTIDs over time as well as the median trajectory.
# p-values represent paired comparisons between wk0/wk2, wk0/wk4, and wk0/wk8
-
##################################################
## How are readouts affected by the inhibitors? ##
##################################################

a <- df_us %>% filter(timepoint == "wk0" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk0")
b <- df_us %>% filter(timepoint == "wk2" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk2")
c <- df_us %>% filter(timepoint == "wk4" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk4")
d <- df_us %>% filter(timepoint == "wk8" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk8")
ggarrange(a,b,c,d, ncol = 2, nrow = 2)

a <- pha_bgsub %>% filter(timepoint == "wk0" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk0")
b <- pha_bgsub %>% filter(timepoint == "wk2" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk2")
c <- pha_bgsub %>% filter(timepoint == "wk4" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk4")
d <- pha_bgsub %>% filter(timepoint == "wk8" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk8")
ggarrange(a,b,c,d, ncol = 2, nrow = 2)

a <- cyt_bgsub %>% filter(timepoint == "wk0" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk0")
b <- cyt_bgsub %>% filter(timepoint == "wk2" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk2")
c <- cyt_bgsub %>% filter(timepoint == "wk4" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk4")
d <- cyt_bgsub %>% filter(timepoint == "wk8" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk8")
ggarrange(a,b,c,d, ncol = 2, nrow = 2)

a <- mtb_bgsub %>% filter(timepoint == "wk0" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk0")
b <- mtb_bgsub %>% filter(timepoint == "wk2" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk2")
c <- mtb_bgsub %>% filter(timepoint == "wk4" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9neg, wk4")
d <- mtb_bgsub %>% filter(timepoint == "wk8" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "cd137", "% cd137 of vg9negw, wk8")
ggarrange(a,b,c,d, ncol = 2, nrow = 2)

a <- mtb_bgsub %>% filter(timepoint == "wk0" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "ifng", "% ifng of vg9neg, wk0")
b <- mtb_bgsub %>% filter(timepoint == "wk2" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "ifng", "% ifng of vg9neg, wk2")
c <- mtb_bgsub %>% filter(timepoint == "wk4" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "ifng", "% ifng of vg9neg, wk4")
d <- mtb_bgsub %>% filter(timepoint == "wk8" & inhib != 'mapk') %>% 
  make_pbmc_plt("vg9neg", "ifng", "% ifng of vg9neg, wk8")
ggarrange(a,b,c,d, ncol = 2, nrow = 2)



###############################################################
## How are readouts induced by cytokine and PHA stimulation? ##
###############################################################

df_temp <- df %>% filter(inhib == 'ni' & tissue == 'pbmc')

# Week 0
a <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = cd69)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

b <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

c <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = cd107a)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

d <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = granzymeb)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

e <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = granzymek)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

f <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

g <- df_temp %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = stim, y = tnf)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

ggarrange(a, b, c, d, e, f, g, ncol = 2, nrow = 4)


# Week 8
a <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = cd69)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  facet_wrap(~subset, ncol = 5)

b <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = cd137)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

c <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = cd107a)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

d <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = granzymeb)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

e <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = granzymek)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

f <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = ifng)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

g <- df_temp %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = stim, y = tnf)) +
  geom_dotplot(position = 'identity') +
  facet_wrap(~subset, ncol = 5)

ggarrange(a, b, c, d, e, f, g, ncol = 2, nrow = 4)




a <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "cd69", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% cd69, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim")

b <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% cd137, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

c <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "cd107a", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% cd107a, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

d <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "granzymeb", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% granzyme b, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

e <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "granzymek", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% granzyme k, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

f <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% ifn-g, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

g <- df_temp %>% filter(timepoint == 'wk8') %>% arrange(stim) %>%
  ggpaired(x = "stim", y = "tnf", id = "ptid", 
           line.color = "black",
           fill = "stim",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d()+
  ggtitle("% tnf, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("stim") 

ggarrange(a, b, c, d, e, f, g, ncol = 4, nrow = 2, common.legend = T)

#################################
## Are the inhibitors working? ##
#################################

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

ggsave('csa_mtb.pdf', width = 10, height = 3, units = 'in')

mtb_bgsub %>% filter(subset == 'vg9neg' & inhib != 'mapk') %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "timepoint",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d(begin = 0.5, end = 0.8) +
  ggtitle("% CD137+")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )



## COMPARE SUBSETS ##
## CD137
a <- mtb_bgsub %>% filter(timepoint == 'wk0', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 0 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

b <- mtb_bgsub %>% filter(timepoint == 'wk2', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 2 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

c <- mtb_bgsub %>% filter(timepoint == 'wk4', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 4 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

d <- mtb_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

ggarrange(a, b, c, d, ncol = 4, nrow = 1)



## IFNg
a <- mtb_bgsub %>% filter(timepoint == 'wk0', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 0 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

b <- mtb_bgsub %>% filter(timepoint == 'wk2', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 2 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

c <- mtb_bgsub %>% filter(timepoint == 'wk4', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 4 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

d <- mtb_bgsub %>% filter(timepoint == 'wk8', inhib != 'mapk') %>% arrange(stim) %>%
  ggpaired(x = "inhib", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "inhib",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, week 8 pbmc")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c("ni", "csa")),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

ggarrange(a, b, c, d, ncol = 4, nrow = 1)


####################################################################
## Are the cells more responsive to PHA or cytokine after IV-BCG? ##
####################################################################

pha_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "cd137", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd137, pbmc + pha")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c('wk0', 'wk8')),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

ggsave(filename= 'CD137_pha.pdf')

pha_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, pbmc + pha")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c('wk0', 'wk8')),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

ggsave(filename= 'ifng_pha.pdf')

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
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c('wk0', 'wk8')),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

ggsave(filename= 'ifng_pha.pdf')

cyt_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "ifng", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% ifng, pbmc + cytok")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c('wk0', 'wk8')),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )

cyt_bgsub %>% filter(inhib == 'ni', tissue == 'pbmc') %>% arrange(stim) %>%
  ggpaired(x = "timepoint", y = "cd69", id = "ptid", 
           line.color = "black",
           fill = "timepoint",
           facet.by = "subset",
           line.size = 0.4, 
           point.size = 2,
           short.panel.labs = FALSE) +
  scale_fill_viridis_d() +
  ggtitle("% cd69, pbmc + cytok")+
  ylab("Frequency within subset") +
  xlab("inhib") +
  theme(legend.position = 'none') +
  stat_compare_means(
    comparisons = list(c('wk0', 'wk8')),
    method = "wilcox.test",
    paired = T,
    show.legend = T
  )


####################################
## How does CsA affect BAL cells? ##
####################################

bal <- df %>% filter(tissue == 'bal' & subset != 'tcrab' & ptid == '36818')
bal$inhib <- factor(bal$inhib, levels = c('ni', 'csa'))

bal_us <- bal %>% filter(stim == 'us')
bal_pha <- bal %>% filter(stim == 'pha')
bal_cyt <- bal %>% filter(stim == 'cytok')
bal_mtb <- bal %>% filter(stim == 'mtbl')

metadata <- bal_pha[,1:15]
bal_pha_bgs <- bal_pha[,16:27]-bal_us[,16:27]
bal_pha_bgs <- cbind(metadata, bal_pha_bgs)
bal_pha_bgs[bal_pha_bgs < 0] <- 0

metadata <- bal_cyt[,1:15]
bal_cyt_bgs <- bal_cyt[,16:27]-bal_us[,16:27]
bal_cyt_bgs <- cbind(metadata, bal_cyt_bgs)
bal_cyt_bgs[bal_cyt_bgs < 0] <- 0

metadata <- bal_mtb[,1:15]
bal_mtb_bgs <- bal_mtb[,16:27]-bal_us[,16:27]
bal_mtb_bgs <- cbind(metadata, bal_mtb_bgs)
bal_mtb_bgs[bal_mtb_bgs < 0] <- 0


## Is CsA working in BAL?
bal_pha_bgs %>% arrange(desc(inhib)) %>% ggplot(aes(x = inhib, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

bal_cyt_bgs %>% arrange(desc(inhib)) %>% ggplot(aes(x = inhib, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

bal_mtb_bgs %>% arrange(desc(inhib)) %>% ggplot(aes(x = inhib, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)

bal_mtb_bgs %>% arrange(desc(inhib)) %>% ggplot(aes(x = inhib, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset, ncol = 5)




#################################################################
## What happens to unstimulated cells treated with inhibitors? ##
#################################################################

df_us$inhib <- factor(df_us$inhib, levels = c('ni', 'csa', 'mapk'))

df_us %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = inhib, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

df_us %>% filter(timepoint == 'wk0') %>% ggplot(aes(x = inhib, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

df_us %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = inhib, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

df_us %>% filter(timepoint == 'wk8') %>% ggplot(aes(x = inhib, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)



################################
## Percent reduction with CsA ##
################################

csa <- mtb_bgsub %>% filter(inhib == 'csa')
md <- csa[,1:15]
csa <- csa[,16:27]
ni <- mtb_bgsub %>% filter(inhib == 'ni')
ni <- ni[,16:27]
freq_down <- ni-csa
pct_down <- (freq_down/ni)*100
pct_down <- cbind(md, pct_down)


pct_down %>% ggplot(aes(x = timepoint, y = hladr)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

pct_down  %>% ggplot(aes(x = timepoint, y = cd69)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

pct_down  %>% ggplot(aes(x = timepoint, y = cd137)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

pct_down  %>% ggplot(aes(x = timepoint, y = cd107a)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

pct_down  %>% ggplot(aes(x = timepoint, y = ifng)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)

pct_down  %>% ggplot(aes(x = timepoint, y = tnf)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~subset)
