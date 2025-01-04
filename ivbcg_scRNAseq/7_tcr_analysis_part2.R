library(dplyr)
library(Seurat)
library(circlize)
library(scales)
library(tidyverse)
library(ggpubr)
library(stringr)

## Define the TCRs that are expanded after IV-BCG ##
dg0d_ct <- read.csv(file = "out/tcr/dg0d_counts.csv")
dg0h_ct <- read.csv(file = "out/tcr/dg0h_counts.csv")
v113_ct <- read.csv(file = "out/tcr/v113_counts.csv")
v015_ct <- read.csv(file = "out/tcr/v015_counts.csv")

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

dg0h_exp <- dg0h_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf") %>% 
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

v113_exp <- v113_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf") %>% 
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

v015_exp <- v015_ct %>% filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf" | ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf") %>% 
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
  filter(ratio_4_to_0 > 2 | ratio_4_to_0 == "Inf") %>%
  filter(Week_4_count >2)

wk8 <- concat %>% filter(expanded == 'Resp') %>%
  filter(ratio_8_to_0 > 2 | ratio_8_to_0 == "Inf") %>%
  filter(Week_8_count>2)



## PLOT THE TCRs
p1 <- concat %>% ggplot(aes(x = Week_0_pct_plus1, y = Week_4_pct_plus1, color = expanded)) +
  geom_point(size = 1) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

p2 <- concat %>% ggplot(aes(x = Week_0_pct_plus1, y = Week_8_pct_plus1, color = expanded)) +
  geom_point(size = 1) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T)


## Test alternative methods

#  Test -- plot percents
test <- concat %>% subset(select = c(Week_0_pct_plus1, Week_4_pct_plus1, Week_8_pct_plus1, expanded))
test <- test %>%
  mutate(log_pct_at_Week_0_plus1 = log(Week_0_pct_plus1),
         log_pct_at_Week_4_plus1 = log(Week_4_pct_plus1),
         log_pct_at_Week_8_plus1 = log(Week_8_pct_plus1)) %>%
  group_by(Week_0_pct_plus1, Week_4_pct_plus1, Week_8_pct_plus1) %>%
  mutate(group_size = n())

# Scale dots by group size
p1 <- test %>% ggplot(aes(x = log_pct_at_Week_0_plus1, y = log_pct_at_Week_4_plus1, color = expanded)) +
  geom_point(aes(size = group_size)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size(range = c(1, 6),
             breaks = c(0, 1, 3, 10, 100, 1000, 10000),
             guide = "legend")

p2 <- test %>% ggplot(aes(x = Week_0_pct_plus1, y = Week_8_pct_plus1, color = expanded)) +
  geom_point(size = 1) +
  scale_x_continuous(trans = log10_trans()) +
  scale_y_continuous(trans = log10_trans()) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1)

ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T)

#  Test -- plot counts
test <- concat %>% subset(select = c(Week_0_count, Week_4_count, Week_8_count, expanded))
test <- test %>%
  mutate(log_Week_0_count_plus1 = log(Week_0_count + 1),
         log_Week_4_count_plus1 = log(Week_4_count + 1),
         log_Week_8_count_plus1 = log(Week_8_count + 1)) %>%
  group_by(Week_0_count, Week_4_count, Week_8_count) %>%
  mutate(group_size = n())


test <- test %>%
  mutate(expanded2 = case_when(
    expanded == 'Resp' ~ "Responder",
    .default="Non-responder"
  ))


# Scale dots by group size
p1 <- test %>% ggplot(aes(x = log_Week_0_count_plus1, y = log_Week_4_count_plus1, color = expanded2)) +
  geom_point(aes(size = group_size)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size(breaks = c(0, 1, 10, 50, 400),
                    transform = 'log',
                    guide = "legend") +
  scale_color_manual(values = c("grey", "black")) + 
  guides(color = guide_legend(nrow = 2))

p2 <- test %>% ggplot(aes(x = log_Week_0_count_plus1, y = log_Week_8_count_plus1, color = expanded2)) +
  geom_point(aes(size = group_size)) +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1) +
  scale_size(breaks = c(0, 1, 10, 50, 400),
             transform = 'log',
             guide = "legend") +
  scale_color_manual(values = c("grey", "black"))

ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = T)
ggsave('out/images/fig4a.pdf', width = 6.5, height = 4, units = 'in')




## VENN DIAGRAMS
install.packages('ggVennDiagram')
library(ggVennDiagram)

# Read in the combined RDS file and tcr files
ivbcg_gd <- readRDS(file = "out/seurat_objects/ivbcg_vd1_integrated.rds")
dg0d_tcr <- read.csv(file = "out/tcr/dg0d_tcr_sch.csv")
dg0h_tcr <- read.csv(file = "out/tcr/dg0h_tcr_sch.csv")
v113_tcr <- read.csv(file = "out/tcr/v113_tcr_sch.csv")
v015_tcr <- read.csv(file = "out/tcr/v015_tcr_sch.csv")

# Merge the TCR CSVs
ivbcg_tcr <- bind_rows(dg0d_tcr, dg0h_tcr, v113_tcr, v015_tcr)

# Add the metadata to each barcode
barcodes <- FetchData(object = ivbcg_gd, vars = c("seurat_clusters", "Timepoint", "Tissue"))
library(tibble)
barcodes <- tibble::rownames_to_column(barcodes, "barcode")
ivbcg_tcr <- merge(x = barcodes, y = ivbcg_tcr, by = "barcode")

# Create the vectors of TCRs from each condition
pbmc_0 <- ivbcg_tcr %>% filter(Timepoint == 'Week_0' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_4 <- ivbcg_tcr %>% filter(Timepoint == 'Week_4' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_8 <- ivbcg_tcr %>% filter(Timepoint == 'Week_8' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
bal_4 <- ivbcg_tcr %>% filter(Timepoint == 'Week_4' & Tissue == 'BAL') %>% pull(cdr3_aa2)
bal_8 <- ivbcg_tcr %>% filter(Timepoint == 'Week_8' & Tissue == 'BAL') %>% pull(cdr3_aa2)

## Generate the Venn Diagrams 
# Overlap between Tissues and Timepoints
p_04 <- list(pbmc_0, pbmc_4)
ggVennDiagram(p_04)

p_08 <- list(pbmc_0, pbmc_8)
ggVennDiagram(p_08)

p_48 <- list(pbmc_4, pbmc_8)
ggVennDiagram(p_48)

p_all <- list(pbmc_0, pbmc_4, pbmc_8)
ggVennDiagram(p_all, category.names = c('PBMC.Pre', "PBMC.4", "PBMC.8"))
ggsave('out/images/pbmc_all.venn.png')

b_all <- list(bal_4, bal_8)
ggVennDiagram(b_all, category.names = c('BAL.4', "BAL.8"))
ggsave('out/images/bal_all.venn.png')

w4 <- list(pbmc_4, bal_4)
ggVennDiagram(w4, category.names = c('PBMC.4', "BAL.4"))
ggsave('out/images/week4_venn.png')

w8 <- list(pbmc_8, bal_8)
ggVennDiagram(w8, category.names = c('PBMC.8', "BAL.8"))
ggsave('out/images/week8_venn.png')


# Overlap in expanded TCRs and BAL
exp_any <- concat %>% filter(expanded == 'Resp')
exp_any <- exp_any %>% separate(CTaa, c("gamma", "delta")) %>% pull(delta)

exp_b4b8 <- list(exp_any, bal_4, bal_8)
ggVennDiagram(exp_b4b8, category.names = c('Responder', "BAL.4", "BAL.8"))
ggsave('out/images/expanded_venn.png')



## Repeat this analysis for each PTID
ivbcg_tcr_dg0d <- ivbcg_tcr %>% filter(sample == 'dg0d')
ivbcg_tcr_dg0h <- ivbcg_tcr %>% filter(sample == 'dg0h')
ivbcg_tcr_v113 <- ivbcg_tcr %>% filter(sample == 'v113')
ivbcg_tcr_v015 <- ivbcg_tcr %>% filter(sample == 'v015')


# DG0D
pbmc_4 <- ivbcg_tcr_dg0d %>% filter(Timepoint == 'Week_4' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_8 <- ivbcg_tcr_dg0d %>% filter(Timepoint == 'Week_8' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
bal_4 <- ivbcg_tcr_dg0d %>% filter(Timepoint == 'Week_4' & Tissue == 'BAL') %>% pull(cdr3_aa2)
bal_8 <- ivbcg_tcr_dg0d %>% filter(Timepoint == 'Week_8' & Tissue == 'BAL') %>% pull(cdr3_aa2)

w4 <- list(pbmc_4, bal_4)
ggVennDiagram(w4, category.names = c('PBMC.4', "BAL.4"))

w8 <- list(pbmc_8, bal_8)
ggVennDiagram(w8, category.names = c('PBMC.8', "BAL.8"))


# DG0H
pbmc_4 <- ivbcg_tcr_dg0h %>% filter(Timepoint == 'Week_4' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_8 <- ivbcg_tcr_dg0h %>% filter(Timepoint == 'Week_8' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
bal_4 <- ivbcg_tcr_dg0h %>% filter(Timepoint == 'Week_4' & Tissue == 'BAL') %>% pull(cdr3_aa2)
bal_8 <- ivbcg_tcr_dg0h %>% filter(Timepoint == 'Week_8' & Tissue == 'BAL') %>% pull(cdr3_aa2)
dg0h_exp_delta <- dg0h_exp %>% separate(CTaa, c("gamma", "delta")) %>% pull(delta)

w4 <- list(pbmc_4, bal_4)
ggVennDiagram(w4, category.names = c('PBMC.4', "BAL.4"))

w8 <- list(pbmc_8, bal_8)
ggVennDiagram(w8, category.names = c('PBMC.8', "BAL.8"))

exp_b4b8 <- list(dg0h_exp_delta, bal_4, bal_8)
ggVennDiagram(exp_b4b8, category.names = c('Responder', "BAL.4", "BAL.8"))

# V113
pbmc_4 <- ivbcg_tcr_v113 %>% filter(Timepoint == 'Week_4' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_8 <- ivbcg_tcr_v113 %>% filter(Timepoint == 'Week_8' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
bal_4 <- ivbcg_tcr_v113 %>% filter(Timepoint == 'Week_4' & Tissue == 'BAL') %>% pull(cdr3_aa2)
bal_8 <- ivbcg_tcr_v113 %>% filter(Timepoint == 'Week_8' & Tissue == 'BAL') %>% pull(cdr3_aa2)
v113_exp_delta <- v113_exp %>% separate(CTaa, c("gamma", "delta")) %>% pull(delta)


w4 <- list(pbmc_4, bal_4)
ggVennDiagram(w4, category.names = c('PBMC.4', "BAL.4"))

w8 <- list(pbmc_8, bal_8)
ggVennDiagram(w8, category.names = c('PBMC.8', "BAL.8"))

exp_b4b8 <- list(dg0h_exp_delta, bal_4, bal_8)
ggVennDiagram(exp_b4b8, category.names = c('Responder', "BAL.4", "BAL.8"))


# V015
pbmc_4 <- ivbcg_tcr_v015 %>% filter(Timepoint == 'Week_4' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
pbmc_8 <- ivbcg_tcr_v015 %>% filter(Timepoint == 'Week_8' & Tissue == 'PBMC') %>% pull(cdr3_aa2)
bal_4 <- ivbcg_tcr_v015 %>% filter(Timepoint == 'Week_4' & Tissue == 'BAL') %>% pull(cdr3_aa2)
bal_8 <- ivbcg_tcr_v015 %>% filter(Timepoint == 'Week_8' & Tissue == 'BAL') %>% pull(cdr3_aa2)
v015_exp_delta <- v015_exp %>% separate(CTaa, c("gamma", "delta")) %>% pull(delta)


w4 <- list(pbmc_4, bal_4)
ggVennDiagram(w4, category.names = c('PBMC.4', "BAL.4"))

w8 <- list(pbmc_8, bal_8)
ggVennDiagram(w8, category.names = c('PBMC.8', "BAL.8"))

exp_b4b8 <- list(dg0h_exp_delta, bal_4, bal_8)
ggVennDiagram(exp_b4b8, category.names = c('Responder', "BAL.4", "BAL.8"))
