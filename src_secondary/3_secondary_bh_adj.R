#### Adjust all pvalues with BH procedure ####
rm(list=ls())

#source(here::here("0-config.R"))

# load all results
H1_res <- readRDS(here('results/unadjusted/H1_whz_res.RDS'))
H2_res <- readRDS(here('results/unadjusted/H2_whz_res.RDS'))
H3_res <- readRDS(here('results/unadjusted/H3_whz_res.RDS'))
H4_res <- readRDS(here('results/unadjusted/H4_whz_res.RDS'))

H1_adj_res <- readRDS(here('results/adjusted/H1_whz_adj_res.RDS'))
H2_adj_res <- readRDS(here('results/adjusted/H2_whz_adj_res.RDS'))
H3_adj_res <- readRDS(here('results/adjusted/H3_whz_adj_res.RDS'))
H4_adj_res <- readRDS(here('results/adjusted/H4_whz_adj_res.RDS'))

H1_res$H = 1
H2_res$H = 2
H3_res$H = 3
H4_res$H = 4

H1_adj_res$H = 1
H2_adj_res$H = 2
H3_adj_res$H = 3
H4_adj_res$H = 4

full_res <- rbind(H1_res, H2_res, H3_res, H4_res)
full_adj_res <- rbind(H1_adj_res, H2_adj_res, H3_adj_res, H4_adj_res)

full_res <- full_res %>% group_by(H) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

full_adj_res <- full_adj_res %>% group_by(H) %>% 
  mutate(BH.Pval=p.adjust(Pval, method="BH")) %>%
  ungroup() %>%
  as.data.frame()

saveRDS(full_res %>% filter(H==1) %>% select(-H), here("results/unadjusted/H1_whz_res.RDS"))
saveRDS(full_res %>% filter(H==2) %>% select(-H), here("results/unadjusted/H2_whz_res.RDS"))
saveRDS(full_res %>% filter(H==3) %>% select(-H), here("results/unadjusted/H3_whz_res.RDS"))
saveRDS(full_res %>% filter(H==4) %>% select(-H), here("results/unadjusted/H4_whz_res.RDS"))

saveRDS(full_adj_res %>% filter(H==1) %>% select(-H), here("results/adjusted/H1_whz_adj_res.RDS"))
saveRDS(full_adj_res %>% filter(H==2) %>% select(-H), here("results/adjusted/H2_whz_adj_res.RDS"))
saveRDS(full_adj_res %>% filter(H==3) %>% select(-H), here("results/adjusted/H3_whz_adj_res.RDS"))
saveRDS(full_adj_res %>% filter(H==4) %>% select(-H), here("results/adjusted/H4_whz)adj_res.RDS"))