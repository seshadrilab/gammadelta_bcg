# Load the necessary libraries
library(tidyverse)
library(dplyr)
library(ggplot2)
library(cowplot)
library(viridis)
library(ggpubr)
library(ggedit)

# Read in the CSVs
df <- read.csv("infant_cordblood_ics.csv")

# Double-check the included PTIDs (sample donors)
df_temp <- df %>%
  select(c("ptid", "tissue")) %>%
  distinct()
df_temp %>% filter(tissue == "pbmc")

## NOTE: CBMC = Cord Blood Mononuclear Cells

############################
## Background subtraction ##
############################

# Calculate the background-subtracted responses to Mtb lysate and PHA stimulation
# Background subtracted responses = (Stimulated_measurement)-(Unstimulated_measurement)
# Omit total AB T cells from the analyses
# Only batches 2 and 3 were stimulated with PHA
# PHA was used as a control to verify that the cells were responsive to general stimulation

df <- df %>% filter(subset != "tcrab")

df$stim <- factor(df$stim, levels = c("unstim", "mtbl", "pha"))
df$subset <- factor(df$subset, levels = c("vg9neg", "vg9", "cd8", "cd4"))
df$batch <- as.character(df$batch)

df_us <- df %>% filter(stim == "unstim")
df_mtb <- df %>% filter(stim == "mtbl")
df_us23 <- df_us %>% filter(batch != "1") # Table of unstim values from batches 2 & 3 to subtract from PHA stimulated values
df_pha <- df %>% filter(stim == "pha")

metadata <- df_mtb[, 1:15]
mtb_bgsub <- df_mtb[, 16:27] - df_us[, 16:27]
mtb_bgsub <- cbind(metadata, mtb_bgsub)
mtb_bgsub[mtb_bgsub < 0] <- 0

metadata <- df_pha[, 1:15]
pha_bgsub <- df_pha[, 16:27] - df_us23[, 16:27]
pha_bgsub <- cbind(metadata, pha_bgsub)
pha_bgsub[pha_bgsub < 0] <- 0


## Calculate the median frequencies in infants and cord blood ##
medians <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  group_by(tissue) %>%
  summarize(across(19:26, median, na.rm = TRUE))

medians_vg9 <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  group_by(tissue) %>%
  summarize(across(19:26, median, na.rm = TRUE))

medians_us <- df_us %>%
  filter(subset == "vg9neg") %>%
  group_by(tissue) %>%
  summarize(across(19:26, median, na.rm = TRUE))



######################
## REPONSES TO MTBL ##
######################

# Generate plots comparing the upregulation of activation markers, cytotoxic proteins, and cytokines between infant pbmc and cord blood
# Each repeated block visualizes a different marker of interest
a <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% HLA-DR+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.signif") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD69+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.signif") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.signif", label.y = 1.25) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 1.5)

d <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.signif", label.y = 0.8) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 1)

e <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 7) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 8)

f <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.signif", label.y = 4) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 4.5)

g <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 0.7) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 0.8)

h <- mtb_bgsub %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 0.7) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  ylim(0, 0.8)

ggarrange(c, d, e, f, g, h, nrow = 2, ncol = 3, common.legend = T)
ggsave("Fig1E.pdf", width = 6, height = 4, units = "in")


########################
## UNSTIMULATED CELLS ##
########################

# Generate plots comparing the expression of activation markers, cytotoxic proteins, and cytokines between infant pbmc and cord blood
# Each repeated block visualizes a different marker of interest

a <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% HLA-DR+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD69+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 1.1) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

d <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 1.2) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

e <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 70) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

f <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 14) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

g <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 0.3) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

h <- df_us %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format", label.y = 0.7) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggarrange(c, d, e, f, g, h, nrow = 2, ncol = 3, common.legend = T)
ggsave("Fig1F.pdf", width = 6, height = 4, units = "in")



##################################
## REPEAT ANALYSIS IN VG9 CELLS ##
##################################

## MTB STIMULATED VG9 CELLS ##

# Generate plots comparing the upregulation of activation markers, cytotoxic proteins, and cytokines between infant pbmc and cord blood
# Each repeated block visualizes a different marker of interest

a <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% HLA-DR+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD69+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

d <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

e <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

f <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

g <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

h <- mtb_bgsub %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggarrange(c, d, e, f, g, h, nrow = 2, ncol = 3, common.legend = T)
ggsave("SuppFig2A.pdf", width = 6, height = 4, units = "in")



## UNSTIMULATED VG9 CELLS ##
# Generate plots comparing the baseline expression of activation markers, cytotoxic proteins, and cytokines between infant pbmc and cord blood
# Each repeated block visualizes a different marker of interest

a <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = hladr)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% HLA-DR+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd69)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD69+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

d <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

e <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

f <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

g <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

h <- df_us %>%
  filter(subset == "vg9") %>%
  ggplot(aes(x = tissue, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = tissue), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggarrange(c, d, e, f, g, h, nrow = 2, ncol = 3, common.legend = T)
ggsave("SuppFig2B.pdf", width = 6, height = 4, units = "in")



## COMPARE VG9+ AND VG9NEG RESPONSES TO STIMULATION WITH MTB WHOLE CELL LYSATE
# Start with cord blood (CBMC)
# Each repeated block visualizes a different marker of interest
a <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

d <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

e <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


f <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "cbmc") %>%
  ggplot(aes(x = subset, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggarrange(a, b, c, d, e, f, nrow = 2, ncol = 3, common.legend = T)
ggsave("SuppFig2C.pdf", width = 6, height = 4, units = "in")


# Next look in BCG-vaccinated infants (PBMC)
# Each repeated block visualizes a different marker of interest

a <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = cd137)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD137+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

b <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = cd107a)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% CD107a+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

c <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = granzymeb)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme B+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

d <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = granzymek)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Granzyme K+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

e <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = ifng)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )


f <- mtb_bgsub %>%
  filter(subset == "vg9" | subset == "vg9neg") %>%
  filter(tissue == "pbmc") %>%
  ggplot(aes(x = subset, y = tnf)) +
  theme_bw() +
  geom_boxplot(aes(color = subset), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.1, end = 0.4) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Arrange the plots and export the PDF
ggarrange(a, b, c, d, e, f, nrow = 2, ncol = 3, common.legend = T)
ggsave("SuppFig2D.pdf", width = 6, height = 4, units = "in")
