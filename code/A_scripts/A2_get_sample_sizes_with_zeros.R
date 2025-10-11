# R Thornley
# 15/07/2025
# get sample sizes with zeros for appendix

# Project: P1_COMPADRE_DRAGNET
# Script: A1_get_sample_sizes_w_zeros

library(tidyverse)
library(ggpubr)


dat_1 <- 
  dat_1 <- read_csv("results/DRAGNet_T0_T1_overlap.csv")  %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("T0", "T1"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))

dat_1 <- dat_1 %>% select(site_name, New_taxon, taxon_site, trt, year_trt, block, new_max_cover) %>%
  group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  pivot_wider(names_from = "year_trt", values_from = "new_max_cover", values_fn = mean) %>% 
  rowwise() %>%
  filter(!all(c_across(T0:T1) == 0)) %>% ungroup() %>%
  pivot_longer(cols = T0:T1, values_to = "new_max_cover", names_to = "year_trt")

zero_counts <- dat_1 %>%
  group_by(trt, year_trt) %>%
  summarise(zero_counts = sum(new_max_cover == 0, na.rm = TRUE))
all_counts <- dat_1 %>%
  group_by(trt, year_trt) %>%
  tally() %>%
  rename("all_counts" = n)

samples_1 <- zero_counts %>% left_join(all_counts) %>%
  mutate(non_zero_counts = all_counts - zero_counts) %>%
  select(trt, year_trt, zero_counts, non_zero_counts, all_counts)

props_1 <- samples_1 %>% 
  mutate(zero_prop = zero_counts/all_counts) %>%
  mutate(non_zero_prop = 1-zero_prop) %>%
  select(trt, year_trt, all_counts, non_zero_prop, zero_prop) %>%
  pivot_longer(cols = non_zero_prop:zero_prop)

#################################################################################

# dat_2
dat_2 <- read_csv("results/DRAGNet_T0_T2_overlap.csv") %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("T0", "T2"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))

dat_2 <- dat_2 %>% select(site_name, New_taxon, taxon_site, trt, year_trt, block, new_max_cover) %>%
  group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  pivot_wider(names_from = "year_trt", values_from = "new_max_cover", values_fn = mean) %>% 
  rowwise() %>%
  filter(!all(c_across(T0:T2) == 0)) %>% ungroup() %>%
  pivot_longer(cols = T0:T2, values_to = "new_max_cover", names_to = "year_trt")

zero_counts <- dat_2 %>%
  group_by(trt, year_trt) %>%
  summarise(zero_counts = sum(new_max_cover == 0, na.rm = TRUE))
all_counts <- dat_2 %>%
  group_by(trt, year_trt) %>%
  tally() %>%
  rename("all_counts" = n)

samples_2 <- zero_counts %>% left_join(all_counts) %>%
  mutate(non_zero_counts = all_counts - zero_counts) %>%
  select(trt, year_trt, zero_counts, non_zero_counts, all_counts)

props_2 <- samples_2 %>% 
  mutate(zero_prop = zero_counts/all_counts) %>%
  mutate(non_zero_prop = 1-zero_prop) %>%
  select(trt, year_trt, all_counts, non_zero_prop, zero_prop) %>%
  pivot_longer(cols = non_zero_prop:zero_prop)

################################################################################

# dat_3
dat_3 <- read_csv("results/DRAGNet_T0_T3_overlap.csv") %>% 
  mutate(trt = factor(trt, levels = c("Control", "Disturbance", "NPK+Disturbance", "NPK", "NPK_Cessation"))) %>%
  filter(!trt == "NPK_Cessation") %>%
  mutate(year_trt = factor(year_trt, levels = c("T0", "T3"))) %>% 
  mutate(site_name = factor(site_name))%>%
  mutate(New_taxon = factor(New_taxon)) %>%
  mutate(taxon_site = paste0(New_taxon, "_", site_name)) %>% 
  mutate(taxon_site = factor(taxon_site)) %>%
  dplyr::select(site_name, New_taxon, taxon_site, year_trt, trt, max_cover, block) %>%
  mutate(new_max_cover = case_when(max_cover == 0 ~ 0,
                                   max_cover > 0 ~ max_cover/100,
                                   TRUE ~ max_cover)) %>%
  mutate(new_max_cover = if_else(new_max_cover == 1, 0.99, new_max_cover))

dat_3 <- dat_3 %>% select(site_name, New_taxon, taxon_site, trt, year_trt, block, new_max_cover) %>%
  group_by(site_name, New_taxon, taxon_site, trt, block) %>%
  pivot_wider(names_from = "year_trt", values_from = "new_max_cover", values_fn = mean) %>% 
  rowwise() %>%
  filter(!all(c_across(T0:T3) == 0)) %>% ungroup() %>%
  pivot_longer(cols = T0:T3, values_to = "new_max_cover", names_to = "year_trt")

zero_counts <- dat_3 %>%
  group_by(trt, year_trt) %>%
  summarise(zero_counts = sum(new_max_cover == 0, na.rm = TRUE))
all_counts <- dat_3 %>%
  group_by(trt, year_trt) %>%
  tally() %>%
  rename("all_counts" = n)

samples_3 <- zero_counts %>% left_join(all_counts) %>%
  mutate(non_zero_counts = all_counts - zero_counts) %>%
  select(trt, year_trt, zero_counts, non_zero_counts, all_counts)

props_3 <- samples_3 %>% 
  mutate(zero_prop = zero_counts/all_counts) %>%
  mutate(non_zero_prop = 1-zero_prop) %>%
  select(trt, year_trt, all_counts, non_zero_prop, zero_prop) %>%
  pivot_longer(cols = non_zero_prop:zero_prop)


################################################################################

p1 <- 
  ggplot(props_1, aes(year_trt, value, fill = name)) + 
  geom_col(position = "stack") + facet_wrap(~trt, ncol = 5) + theme_bw()+
  scale_fill_manual(values = c("non_zero_prop" = "grey40", "zero_prop" = "grey80"),
                    labels = c("non_zero_prop" = "Non Zeros", "zero_prop" = "Zeros"))+
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill = "white", color = "black"),  
        strip.text = element_text(color = "black"))+
  ylab("Proportion of data")+
  xlab("Year treatment")+
  ggtitle("Proportion of zeros for T0-T1 Analysis")
p1

p2 <- 
  ggplot(props_2, aes(year_trt, value, fill = name)) + 
  geom_col(position = "stack") + facet_wrap(~trt, ncol = 5) + theme_bw()+
  scale_fill_manual(values = c("non_zero_prop" = "grey40", "zero_prop" = "grey80"),
                    labels = c("non_zero_prop" = "Non Zeros", "zero_prop" = "Zeros"))+
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill = "white", color = "black"),  
         strip.text = element_text(color = "black"))+             
  ylab("Proportion of data")+
  xlab("Year treatment")+
  ggtitle("Proportion of zeros for T0-T2 Analysis")
p2

p3 <- 
  ggplot(props_3, aes(year_trt, value, fill = name)) + 
  geom_col(position = "stack") + facet_wrap(~trt, ncol = 5) + theme_bw()+
  scale_fill_manual(values = c("non_zero_prop" = "grey40", "zero_prop" = "grey80"),
                    labels = c("non_zero_prop" = "Non Zeros", "zero_prop" = "Zeros"))+
  theme(legend.title = element_blank())+
  theme(strip.background = element_rect(fill = "white", color = "black"),  
        strip.text = element_text(color = "black"))+             
  ylab("Proportion of data")+
  xlab("Year treatment")+
  ggtitle("Proportion of zeros for T0-T3 Analysis")
p3

# now create some nice tables to go with the barcharts
# tidy up names

names(samples_1) <- c("Treatment", "Year", "Number of \nzero samples", "Number of \nnon-zero samples", "Total \nsample number")
names(samples_2) <- c("Treatment", "Year", "Number of \nzero samples", "Number of \nnon-zero samples", "Total \nsample number")
names(samples_3) <- c("Treatment", "Year", "Number of \nzero samples", "Number of \nnon-zero samples", "Total \nsample number")
t1 <- samples_1 %>% ggtexttable(theme = ttheme("minimal"), rows = NULL)
t2 <- samples_2 %>% ggtexttable(theme = ttheme("minimal"), rows = NULL)
t3 <- samples_3 %>% ggtexttable(theme = ttheme("minimal"), rows = NULL)

all <- ggarrange(p1, p2, p3, t1, t2, t3, common.legend = TRUE)
all

ggsave("figures/A1_sample_size_zeros.jpeg", all, height = 7.5, width = 18)

all <- ggarrange(p1, t1, p2, t2, p3, t3, common.legend = TRUE, ncol = 2, nrow = 3)
all

ggsave("figures/A1_sample_size_zeros_again.jpeg", all, height = 14, width = 12)

################################################################################

