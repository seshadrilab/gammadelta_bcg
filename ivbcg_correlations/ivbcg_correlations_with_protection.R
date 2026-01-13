# Load the necessary libraries
library(tidyverse)
library(cowplot)
library(viridis)
library(ggpubr)
library(zoo)



############
## SET-UP ##
############

# Read in the table
df <- read.csv("vd1_protection_correlations.csv")

# Add a CFU + 1 column to allow logarithmic transformation
df <- df %>% mutate(cfu_nex_p1 = cfu_nex + 1)

# Annotate the protected animals according to the definition in PMID: 37267955 (<100 thoracic CFU)
df$protection <- ifelse(df$cfu_nex_p1 < 101, "protected", "unprotected")
df$protection <- factor(df$protection, levels = c("unprotected", "protected"))
df$workspace <- as.character(df$workspace)

# Change timepoints to class numeric
class(df$timepoint) <- "numeric"

# Filter out irrelevant data (Keep only Vg9neg, overnight stimulation, and assigned protection status)
df <- df %>% filter(subset == "vg9neg" & stim_length == "O/N" & !is.na(protection))
df$ptid %>% unique() # 34 PTIDs in total

# Log-transform CFU
df <- df %>% mutate(cfu_nex_p1_log = log(cfu_nex_p1))

# Generate tables for BG (Unstimulated/background) and PPD
df_us <- df %>% filter(stim == "BG")
df_ppd <- df %>% filter(stim == "PPD")

# Calculate background-subtracted responses to stimulation
df_ppd <- df_us %>%
  subset(select = c("ptid", "timepoint", "icos", "ifng", "il2", "il17", "mr1", "il21", "tnf", "all_cytokine")) %>%
  inner_join(df_ppd, by = c("ptid", "timepoint"), keep = FALSE)

df_ppd <- df_ppd %>% mutate(
  icos_bgs = icos.y - icos.x,
  ifng_bgs = ifng.y - ifng.x,
  il2_bgs = il2.y - il2.x,
  il17_bgs = il17.y - il17.x,
  mr1_bgs = mr1.y - mr1.x,
  il21_bgs = il21.y - il21.x,
  tnf_bgs = tnf.y - tnf.x,
  cytok_bgs = all_cytokine.y - all_cytokine.x
)
df_ppd[df_ppd < 0] <- 0


## EXPORT THE TABLES ##
write.csv(df_ppd, "df_ppd.csv")
write.csv(df_us, "df_us.csv")


# Filter out samples where Vg9neg count <100
# NOTE: No PTIDs were excluded entirely. Some timepoints were excluded for some PTIDS.
ppd_100 <- df_ppd %>% filter(vg9neg_count > 99)
us_100 <- df_us %>% filter(vg9neg_count > 99)


## Calculate the median frequencies by group ##

us_medians <- us_100 %>%
  filter(subset == "vg9neg") %>%
  group_by(timepoint, protection) %>%
  summarize(across(19:30, median, na.rm = TRUE))

ppd_medians <- ppd_100 %>%
  filter(subset == "vg9neg") %>%
  select(c(ptid, timepoint, protection, cytok_bgs)) %>%
  group_by(timepoint, protection) %>%
  summarize(across(cytok_bgs, median, na.rm = TRUE))


#############
## T-tests ##
#############

## Response to PPD'
# The immune measurement changes in each repeated code block
a <- ppd_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, ifng_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- ppd_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, il2_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IL-2+ in response to PPD") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- ppd_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, tnf_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% TNF+ in response to PPD") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

d <- ppd_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, cytok_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Cytokine+ in response to PPD") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a, b, c, d, nrow = 4, ncol = 1, legend = "none")

ggsave("all_cytokine_correlations.pdf", width = 7.5, height = 12, units = "in")
ggsave("Fig6C_right.pdf", plot = d, width = 8.5, height = 3, units = "in")


## Unstimulated cells
# The immune measurement changes in each repeated code block
a <- us_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, non_naive)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Non-naive") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- us_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- us_100 %>%
  filter(subset == "vg9neg") %>%
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell") +
  stat_compare_means(method = "t.test", paired = F, label = "p.format") +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a, b, nrow = 2, ncol = 1, legend = "none")

ggsave("Fig6C_left.pdf", plot = b, width = 8.5, height = 3, units = "in")
