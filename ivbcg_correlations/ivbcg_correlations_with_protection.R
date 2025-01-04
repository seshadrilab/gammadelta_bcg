
# Load the necessary libraries
library(tidyverse)
library(cowplot)
library(viridis)
library(ggpubr)
library(zoo)

############
## SET-UP ##
############
  
# Read in the tables
df <- read.csv("table_concat_vd1.csv")
df <- df %>% mutate(cfu_nex_p1 = cfu_nex + 1)
df$protection <- ifelse(df$cfu_nex_p1 < 101, 'protected', 'unprotected')
df$protection <- factor(df$protection, levels = c("unprotected", "protected"))
df$workspace <- as.character(df$workspace)

# Cohort metadata
df_md <- df[, 1:16] %>%
  subset(select = -c(timepoint, stim, stim_length, tissue, workspace, date_acquired, tube_id)) %>%
  distinct()
df_md$protection <- ifelse(df_md$cfu_nex < 100, 'protected', 'unprotected')

write.csv(df_md, 'cohort_metadata_table.csv')

df_md$age_at_vax <- as.numeric(df_md$age_at_vax)
df_md$actual_dose <- as.numeric(df_md$actual_dose)

df_md %>%
  group_by(protection) %>%
  summarize(count = n(),
            age = mean(age_at_vax))

df_md %>%
  group_by(protection, vax_cohort) %>%
  summarize(count_cohort = n())

df_md %>%
  group_by(protection, sex) %>%
  summarize(count_sex = n())

df_md %>%
  ggplot(aes(x = actual_dose, fill = vax_cohort, color = vax_cohort)) +
  geom_histogram(alpha = 0.4, binwidth = 1/6, position = 'identity') +
  scale_x_log10() +
  theme_bw()


# Filter out samples where Vg9neg count <100
df <- df %>% filter(vg9neg_count > 99)

# Change timepoints to numerical data
class(df$timepoint) <- 'numeric'

# Filter out irrelevent data
df <- df %>% filter(subset == 'vg9neg' & stim_length == 'O/N' & !is.na(protection))

# Transform CFU
df <- df %>% mutate(cfu_nex_p1_log = log(cfu_nex_p1))

# Generate tables for BG, PPD, and Mtb300
df_us <- df %>% filter(stim == 'BG')
df_ppd <- df %>% filter(stim == 'PPD')
df_mtb <- df %>% filter(stim == 'Mtb300')

# Calculate background-subtracted responses to stimulation
df_ppd <- df_us %>% subset(select = c('ptid', 'timepoint', 'icos', 'ifng', 'il2', 'il17', 'mr1', 'il21', 'tnf', 'all_cytokine')) %>%
  inner_join(df_ppd, by = c('ptid', 'timepoint'), keep = FALSE)

df_ppd <- df_ppd %>% mutate(icos_bgs = icos.y - icos.x,
                            ifng_bgs = ifng.y - ifng.x,
                            il2_bgs = il2.y - il2.x,
                            il17_bgs = il17.y - il17.x,
                            mr1_bgs = mr1.y - mr1.x,
                            il21_bgs = il21.y - il21.x,
                            tnf_bgs = tnf.y - tnf.x,
                            cytok_bgs = all_cytokine.y - all_cytokine.x)
df_ppd[df_ppd < 0] <- 0


df_mtb <- df_us %>% subset(select = c('ptid', 'timepoint', 'icos', 'ifng', 'il2', 'il17', 'mr1', 'il21', 'tnf')) %>%
  inner_join(df_mtb, by = c('ptid', 'timepoint'), keep = FALSE)

df_mtb <- df_mtb %>% mutate(icos_bgs = icos.y - icos.x,
                            ifng_bgs = ifng.y - ifng.x,
                            il2_bgs = il2.y - il2.x,
                            il17_bgs = il17.y - il17.x,
                            mr1_bgs = mr1.y - mr1.x,
                            il21_bgs = il21.y - il21.x,
                            tnf_bgs = tnf.y - tnf.x)
df_mtb[df_mtb < 0] <- 0


## Calculate the median frequencies by group ##

us_medians <- df_us %>%
  filter(subset == 'vg9neg') %>%
  group_by(timepoint, protection) %>%
  summarize(across(19:30, median, na.rm = TRUE))

ppd_medians <- df_ppd %>%
  filter(subset == 'vg9neg') %>%
  select(c(ptid, timepoint, protection, cytok_bgs)) %>%
  group_by(timepoint, protection) %>%
  summarize(across(cytok_bgs, median, na.rm = TRUE))
 


##################################################
## PEARSON CORRELATIONS -- ALL BATCHES COMBINED ##
##################################################

## FUNCTIONAL RESPONSE TO PPD ##

p1 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p2 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p3 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson')  +
  scale_color_viridis(discrete = T, end = 0.8)


p4 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p5 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p6 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p7 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p8 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p9 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p10 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p11 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p12 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p13 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p14 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p15 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p16 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p17 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p18 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)



## FUNCTIONAL RESPONSE TO MTB300 ##

p1 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p2 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p3 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p4 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p5 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p6 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p7 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p8 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p9 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p10 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p11 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p12 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p13 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p14 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p15 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p16 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p17 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p18 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson')

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)



## UNSTIMULATED CELLS ##

p1 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p2 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p3 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p4 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p5 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p6 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p7 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p8 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p9 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p10 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p11 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p12 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p13 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p14 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p15 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p16 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p17 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p18 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p19 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p20 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p21 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p22 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p23 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p24 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p25 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p26 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p27 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p28 <- df_us %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk0") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p29 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p30 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("%Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p31 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("%Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)
ggarrange(p19, p20, p21, p22, p23, p24, p25, p26, p27, ncol = 3, nrow = 3)
ggarrange(p28, p29, p30, p31, ncol = 2, nrow = 2)


## BASELINE RESPONSES ##

p1 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p2 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p3 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p4 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p5 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p6 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)




###########################################
## PEARSON CORRELATIONS -- BATCH EFFECTS ##
###########################################

## FUNCTIONAL RESPONSE TO PPD ##

p1 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p2 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p3 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson')  +
  scale_color_viridis(discrete = T, end = 0.8)


p4 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p5 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p6 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p7 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p8 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p9 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p10 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p11 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p12 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p13 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p14 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p15 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p16 <- df_ppd %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p17 <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p18 <- df_ppd %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)



## FUNCTIONAL RESPONSE TO MTB300 ##

p1 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p2 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p3 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p4 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p5 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p6 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p7 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p8 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p9 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p10 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p11 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p12 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p13 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p14 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p15 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson')


p16 <- df_mtb %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson')

p17 <- df_mtb %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson')

p18 <- df_mtb %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson')

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)



## UNSTIMULATED CELLS ##

p1 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p2 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p3 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p4 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p5 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p6 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p7 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p8 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p9 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p10 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p11 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p12 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p13 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p14 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p15 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p16 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p17 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p18 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p19 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p20 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p21 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, non_naive, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("Memory of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p22 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p23 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p24 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, mr1, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MR1 of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p25 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p26 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p27 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, trm, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TRM of Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)


p28 <- df_us %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk0") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p29 <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p30 <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("%Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

p31 <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, pct_of_tcell, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("%Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
ggarrange(p7, p8, p9, p10, p11, p12, ncol = 3, nrow = 2)
ggarrange(p13, p14, p15, p16, p17, p18, ncol = 3, nrow = 2)
ggarrange(p19, p20, p21, p22, p23, p24, p25, p26, p27, ncol = 3, nrow = 3)
ggarrange(p28, p29, p30, p31, ncol = 2, nrow = 2)


## BASELINE RESPONSES ##

p1 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, ifng_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p2 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il2_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p3 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il17_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-17 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p4 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, il21_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-21 of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p5 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, tnf_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

p6 <- df_ppd %>% filter(timepoint == '0') %>%
  group_by(ptid) %>%
  ggplot(aes(cfu_nex_p1_log, icos_bgs, color = workspace)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("ICOS of Vg9neg, wk0") +
  stat_cor(method = 'pearson')

ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)





#############
## T-tests ##
#############

## Response to PPD
a <- df_ppd %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, ifng_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- df_ppd %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, il2_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IL-2+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- df_ppd %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, tnf_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% TNF+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

d <- df_ppd %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, cytok_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Cytokine+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a,b,c,d, nrow = 4, ncol = 1, legend = 'none')

ggsave('all_cytokine_correlations.pdf', width = 7.5, height = 12, units = 'in')
ggsave('Fig7C.pdf', plot = d, width = 8.5, height = 3, units = 'in')


## Response to Mtb300
a <- df_mtb %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, ifng_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- df_mtb %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, il2_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- df_mtb %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, tnf_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a,b,c, nrow = 3, ncol = 1, legend = 'none')


## Unstimulated cells
a <- df_us %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, non_naive)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Non-naive") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- df_us %>% filter(subset == 'vg9neg') %>% 
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a,b, nrow = 2, ncol = 1, legend = 'none')

ggsave('Fig7B.pdf', plot = b, width = 8.5, height = 3, units = 'in')


## Potential confounders
# Dose of IV-BCG
a <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, ifng_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

b <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, il2_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

c <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, tnf_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(a,b,c, ncol = 3, nrow = 1)


a <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

b <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

c <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(actual_dose, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(a,b,c, ncol = 3, nrow = 1)

# Sex

a <- df_ppd %>% filter(subset == 'vg9neg', timepoint == '4') %>% 
  ggplot(aes(sex, ifng_bgs)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- df_ppd %>% filter(subset == 'vg9neg', timepoint == '4') %>% 
  ggplot(aes(sex, il2_bgs)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% IL-2+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- df_ppd %>% filter(subset == 'vg9neg', timepoint == '4') %>% 
  ggplot(aes(sex, tnf_bgs)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% TNF+ in response to PPD") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a,b,c, legend = 'none', ncol = 3, nrow = 1)

a <- df_us %>% filter(subset == 'vg9neg', timepoint == '2') %>% 
  ggplot(aes(sex, pct_of_tcell)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg, week 2") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

b <- df_us %>% filter(subset == 'vg9neg', timepoint == '4') %>% 
  ggplot(aes(sex, pct_of_tcell)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg, week 4") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

c <- df_us %>% filter(subset == 'vg9neg', timepoint == '8') %>% 
  ggplot(aes(sex, pct_of_tcell)) +
  geom_boxplot(aes(color = sex), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.5, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg, week 8") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  facet_wrap(~timepoint, ncol = 4)

ggarrange(a,b,c, legend = 'none', ncol = 3, nrow = 1)


# Age
df_ppd$age_at_vax <- as.numeric(df_ppd$age_at_vax)
df_us$age_at_vax <- as.numeric(df_us$age_at_vax)

a <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, ifng_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IFNg of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

b <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, il2_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("IL-2 of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

c <- df_ppd %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, tnf_bgs)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("TNF of Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(a,b,c, ncol = 3, nrow = 1)

a <- df_us %>% filter(timepoint == '2') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk2") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

b <- df_us %>% filter(timepoint == '4') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk4") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

c <- df_us %>% filter(timepoint == '8') %>%
  group_by(ptid) %>%
  ggplot(aes(age_at_vax, pct_of_tcell)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("% Vg9neg, wk8") +
  stat_cor(method = 'pearson') +
  scale_color_viridis(discrete = T, end = 0.8)

ggarrange(a,b,c, ncol = 3, nrow = 1)


## VG9NEG FREQUENCY ##

a <- df_us %>% filter(timepoint == '2' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell, Wk2") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

b <- df_us %>% filter(timepoint == '4' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell, Wk4") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

c <- df_us %>% filter(timepoint == '8' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell, Wk8") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

d <- df_us %>% filter(timepoint == '0' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, pct_of_tcell)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% Vg9neg of T cell, Wk0") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


## FUNCTIONAL RESPONSE TO PPD AT WEEK 4 ##

d <- df_ppd %>% filter(timepoint == '4' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, ifng_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IFN-g+, Wk4") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

e <- df_ppd %>% filter(timepoint == '4' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, il2_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% IL-2+, Wk4") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

f <- df_ppd %>% filter(timepoint == '4' & subset == 'vg9neg') %>% 
  ggplot(aes(protection, tnf_bgs)) +
  geom_boxplot(aes(color = protection), width = 0.5, lwd = 1) +
  scale_color_viridis(discrete = T, begin = 0.5, end = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.75, position = position_jitter(0.05)) +
  ggtitle("% TNF+, Wk4") +
  stat_compare_means(method = 't.test', paired = F, label = 'p.format') +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

ggarrange(a,b,c,d,e,f, ncol = 3, nrow = 2, common.legend = T)
