
# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggedit)

# Read in the CSVs
df <- read.csv("table_concat.csv")


############################
## Background subtraction ##
############################

df <- df %>% filter(subset != 'tcrab')

df$stim <- factor(df$stim, levels = c('unstim', 'mtbl', 'pha'))
df$subset <- factor(df$subset, levels = c('vg9neg', 'vg9', 'cd8', 'cd4'))
df$batch <- as.character(df$batch)

df_us <- df %>% filter(stim == 'unstim')
df_mtb <- df %>% filter(stim == 'mtbl')
df_us23 <- df_us %>% filter(batch != '1')
df_pha <- df %>% filter(stim == 'pha')

metadata <- df_mtb[,1:15]
mtb_bgsub <- df_mtb[,16:27]-df_us[,16:27]
mtb_bgsub <- cbind(metadata, mtb_bgsub)
mtb_bgsub[mtb_bgsub < 0] <- 0

metadata <- df_pha[,1:15]
pha_bgsub <- df_pha[,16:27]-df_us23[,16:27]
pha_bgsub <- cbind(metadata, pha_bgsub)
pha_bgsub[pha_bgsub < 0] <- 0


## Calculate the median frequencies in infants and cord blood ##

medians <- mtb_bgsub %>%
  filter(subset == 'vg9neg') %>%
  group_by(tissue) %>%
  summarize(across(19:26, median, na.rm = TRUE))



######################
## REPONSES TO MTBL ##
######################

# Generate plots comparing the upregulation of activation markers, cytotoxic proteins, and cytokines
# between infant pbmc and cord blood

a <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% HLA-DR+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.signif') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

b <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD69+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.signif') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

c <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD137+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.signif', label.y = 1.25) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,1.5)

d <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD107a+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.signif', label.y = 0.8) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,1)

e <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme B+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,8)

f <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme K+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.signif', label.y = 4) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,4.5)

g <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% IFN-g+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 0.7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,0.8)

h <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% TNF+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 0.7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  ylim(0,0.8)

ggarrange(c,d,e,f,g,h, nrow = 2, ncol = 3, common.legend = T)
ggsave("Fig1E_pformat.pdf", width = 6, height = 4, units = 'in')




########################
## UNSTIMULATED CELLS ##
########################

# Generate plots comparing the baseline expression of activation markers, cytotoxic proteins, and cytokines
# between infant pbmc and cord blood

a <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% HLA-DR+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

b <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD69+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

c <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD137+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 1.1) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

d <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD107a+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 1.2) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

e <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme B+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 70) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

f <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme K+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 14) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% IFN-g+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

h <- df_us %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% TNF+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format', label.y = 0.7) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggarrange(c,d,e,f,g,h, nrow = 2, ncol = 3, common.legend = T)
ggsave("Fig1F_pformat.pdf", width = 6, height = 4, units = 'in')







##########################
## NOT USED FOR FIGURES ##
##########################


# Generate plots comparing T cell populations
a <- mtb_bgsub %>% ggplot(aes(x = tissue, y = hladr, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% HLA-DR up') +
  stat_compare_means(method = 't.test', label = 'p.format')

b <- mtb_bgsub %>% ggplot(aes(x = tissue, y = cd69, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD69 up') +
  stat_compare_means(method = 't.test', label = 'p.format')

c <- mtb_bgsub %>% ggplot(aes(x = tissue, y = cd137, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD137 up') +
  stat_compare_means(method = 't.test', label = 'p.format')

d <- mtb_bgsub %>% ggplot(aes(x = tissue, y = cd107a, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD107a up') +
  stat_compare_means(method = 't.test', label = 'p.format')

e <- mtb_bgsub %>% ggplot(aes(x = tissue, y = granzymeb, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% Granzyme B up') +
  stat_compare_means(method = 't.test', label = 'p.format')

f <- mtb_bgsub %>% ggplot(aes(x = tissue, y = granzymek, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% GranzymeK up') +
  stat_compare_means(method = 't.test', label = 'p.format')

g <- mtb_bgsub %>% ggplot(aes(x = tissue, y = ifng, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% IFN-g up') +
  stat_compare_means(method = 't.test', label = 'p.format')

h <- mtb_bgsub %>% ggplot(aes(x = tissue, y = tnf, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% TNF up') +
  stat_compare_means(method = 't.test', label = 'p.format')

ggarrange (a,b,c,d, nrow = 2, ncol = 2, common.legend = T)
ggarrange (e,f,g,h, nrow = 2, ncol = 2, common.legend = T)




######################
## RESPONSES TO PHA ##
######################


# Generate plots comparing the upregulation of activation markers, cytotoxic proteins, and cytokines
# between infant pbmc and cord blood

a <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% HLA-DR+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

b <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD69+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

c <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD137+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

d <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% CD107a+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

e <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme B+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

f <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% Granzyme K+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% IFN-g+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

h <- pha_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle('% TNF+') +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggarrange(a,b,c,d,e,f,g,h, nrow = 2, ncol = 4, common.legend = T)

ggarrange(c,d,e,f,g,h, nrow = 2, ncol = 3, common.legend = T)

# Generate plots comparing T cell populations
a <- pha_bgsub %>% ggplot(aes(x = tissue, y = hladr, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% HLA-DR up') +
  stat_compare_means(method = 't.test', label = 'p.format')

b <- pha_bgsub %>% ggplot(aes(x = tissue, y = cd69, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD69 up') +
  stat_compare_means(method = 't.test', label = 'p.format')

c <- pha_bgsub %>% ggplot(aes(x = tissue, y = cd137, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD137 up') +
  stat_compare_means(method = 't.test', label = 'p.format')

d <- pha_bgsub %>% ggplot(aes(x = tissue, y = cd107a, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% CD107a up') +
  stat_compare_means(method = 't.test', label = 'p.format')

e <- pha_bgsub %>% ggplot(aes(x = tissue, y = granzymeb, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% Granzyme B up') +
  stat_compare_means(method = 't.test', label = 'p.format')

f <- pha_bgsub %>% ggplot(aes(x = tissue, y = granzymek, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% GranzymeK up') +
  stat_compare_means(method = 't.test', label = 'p.format')

g <- pha_bgsub %>% ggplot(aes(x = tissue, y = ifng, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% IFN-g up') +
  stat_compare_means(method = 't.test', label = 'p.format')

h <- pha_bgsub %>% ggplot(aes(x = tissue, y = tnf, fill = tissue)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  facet_wrap(~subset, ncol = 4) +
  scale_color_grey() +
  stat_summary(fun = 'median',
               size = 0.5,
               geom = 'crossbar',
               color = 'black') +
  ggtitle('% TNF up') +
  stat_compare_means(method = 't.test', label = 'p.format')

ggarrange (a,b,c,d, nrow = 2, ncol = 2, common.legend = T)
ggarrange (e,f,g,h, nrow = 2, ncol = 2, common.legend = T)



########################
## VIEW BATCH EFFECTS ##
########################

# Generate plots comparing the upregulation of activation markers, cytotoxic proteins, and cytokines
# between infant pbmc and cord blood
a <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = hladr, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
  ggtitle('% HLA-DR up') 

b <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd69, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
   ggtitle('% CD69 up')

c <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd137, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
   ggtitle('% CD137 up') 

d <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = cd107a, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
    ggtitle('% CD107a up') 

e <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymeb, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
   ggtitle('% Granzyme B up') 

f <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = granzymek, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
   ggtitle('% Granzyme K up') 

g <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = ifng, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
   ggtitle('% IFN-g up') 

h <- mtb_bgsub %>% filter (subset == 'vg9neg') %>% ggplot(aes(x = tissue, y = tnf, fill = batch)) +
  theme_bw() +
  scale_fill_manual(values = c('grey', 'black')) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               dotsize = 1) +
  scale_color_grey() +
    ggtitle('% TNF up') 

ggarrange(a,b,c,d,e,f,g,h, nrow = 2, ncol = 4, common.legend = T)

