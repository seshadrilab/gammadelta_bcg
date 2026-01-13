library(dplyr)
library(circlize)
library(scales)
library(tidyverse)
library(ggpubr)
library(stringr)

# Read in the TCR count dataframes for each sample donor
dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")

## Define the TCRs that are expanded after IV-BCG ##
#DG0H
dg0h_ct <- dg0h_ct %>% subset(select = -c(X))
dg0h_ct <- dg0h_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
dg0h_ct[is.na(dg0h_ct)] <- 0
dg0h_ct <- dg0h_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct = Week_0_count/sum(dg0h_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(dg0h_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(dg0h_ct$Week_8_count)*100)
dg0h_ct <- dg0h_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)

dg0h_ct <- dg0h_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1,
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)

dg0h_exp <- dg0h_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf") %>% 
  filter(Week_4_count>2 | Week_8_count>2)

dg0h_pre <- dg0h_ct %>% filter(Week_0_count>1) %>%
  anti_join(dg0h_exp, by = 'CTaa')

dg0h_ct_new <- dg0h_ct %>%
  mutate(expanded = case_when(
    CTaa %in% dg0h_exp$CTaa ~ "Resp",
    CTaa %in% dg0h_pre$CTaa ~ "NR_exp",
    .default="NR_nonexp"
  ))



#V113
v113_ct <- v113_ct %>% subset(select = -c(X, Tissue, cdr3_aa2)) %>% unique()
v113_ct <- v113_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v113_ct[is.na(v113_ct)] <- 0
v113_ct <- v113_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
v113_ct <- v113_ct %>% mutate(Week_0_pct = Week_0_count/sum(v113_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v113_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v113_ct$Week_8_count)*100)
v113_ct <- v113_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)

v113_ct <- v113_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1,
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)

v113_exp <- v113_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf") %>% 
  filter(Week_4_count>2 | Week_8_count>2)

v113_pre <- v113_ct %>% filter(Week_0_count>1) %>%
  anti_join(v113_exp, by = 'CTaa')

v113_ct_new <- v113_ct %>%
  mutate(expanded = case_when(
    CTaa %in% v113_exp$CTaa ~ "Resp",
    CTaa %in% v113_pre$CTaa ~ "NR_exp",
    .default="NR_nonexp"
  ))




#V015
v015_ct <- v015_ct %>% subset(select = -c(X))
v015_ct <- v015_ct %>% pivot_wider(names_from = Timepoint, values_from = Freq)
v015_ct[is.na(v015_ct)] <- 0
v015_ct <- v015_ct %>% rename(Week_0_count = Week_0, Week_4_count = Week_4, Week_8_count = Week_8)
v015_ct <- v015_ct %>% mutate(Week_0_pct = Week_0_count/sum(v015_ct$Week_0_count)*100,
                              Week_4_pct = Week_4_count/sum(v015_ct$Week_4_count)*100,
                              Week_8_pct = Week_8_count/sum(v015_ct$Week_8_count)*100)
v015_ct <- v015_ct %>% mutate(ratio_4_to_0 = Week_4_pct/Week_0_pct,
                              ratio_8_to_0 = Week_8_pct/Week_0_pct)

v015_ct <- v015_ct %>% mutate(Week_0_pct_plus1 = Week_0_pct +1,
                              Week_4_pct_plus1 = Week_4_pct +1,
                              Week_8_pct_plus1 = Week_8_pct +1)

v015_exp <- v015_ct %>% filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf") %>% 
  filter(Week_4_count>2 | Week_8_count>2)

v015_pre <- v015_ct %>% filter(Week_0_count>1) %>%
  anti_join(v015_exp, by = 'CTaa')

v015_ct_new <- v015_ct %>%
  mutate(expanded = case_when(
    CTaa %in% v015_exp$CTaa ~ "Resp",
    CTaa %in% v015_pre$CTaa ~ "NR_exp",
    .default="NR_nonexp"
  ))


## Concatenate the tables
concat <- rbind(dg0h_ct_new, v015_ct_new, v113_ct_new)

## Enumerate the expanded TCRs at each timepoint
wk4 <- concat %>% filter(expanded == 'Resp') %>%
  filter(ratio_4_to_0 > 3 | ratio_4_to_0 == "Inf") %>%
  filter(Week_4_count > 3)

wk8 <- concat %>% filter(expanded == 'Resp') %>%
  filter(ratio_8_to_0 > 3 | ratio_8_to_0 == "Inf") %>%
  filter(Week_8_count > 3)


## Visualize a summary of the TCR frequencies over time
# Generate a new dataframe containing log-transformed clonotype frequencies at each time point
temp <- concat %>% subset(select = c(Week_0_pct_plus1, Week_4_pct_plus1, Week_8_pct_plus1, expanded))
temp <- test %>%
  mutate(log_pct_at_Week_0_plus1 = log(Week_0_pct_plus1),
         log_pct_at_Week_4_plus1 = log(Week_4_pct_plus1),
         log_pct_at_Week_8_plus1 = log(Week_8_pct_plus1)) %>%
  group_by(Week_0_pct_plus1, Week_4_pct_plus1, Week_8_pct_plus1) %>%
  mutate(group_size = n())


# Generate a plot of the clonotype frequencies pre-vaccination vs. the clonotype frequencies at week 4
# The log-transformed frequencies are used
# The size of each dot corresponds to the number of clonotypes occupying that space on the plot
p1 <- temp %>% ggplot(aes(x = log_Week_0_count_plus1, y = log_Week_4_count_plus1, color = expanded2)) +
  geom_point(aes(size = group_size)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size(breaks = c(0, 1, 10, 50, 400),
                    transform = 'log',
                    guide = "legend") +
  scale_color_manual(values = c("grey", "black")) + 
  guides(color = guide_legend(nrow = 2))

# Generate a plot of the clonotype frequencies pre-vaccination vs. the clonotype frequencies at week 8
# The log-transformed frequencies are used
# The size of each dot corresponds to the number of clonotypes occupying that space on the plot
p2 <- temp %>% ggplot(aes(x = log_Week_0_count_plus1, y = log_Week_8_count_plus1, color = expanded2)) +
  geom_point(aes(size = group_size)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size(breaks = c(0, 1, 10, 50, 400),
             transform = 'log',
             guide = "legend") +
  scale_color_manual(values = c("grey", "black"))

# Arrange both plots into a single figure
ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T)
ggsave('out/images/fig4a_3x.pdf', width = 6.5, height = 4, units = 'in')


