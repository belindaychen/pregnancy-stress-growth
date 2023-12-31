rm(list=ls())

source(here::here("0-config.R"))

d<-readRDS(paste0(dropboxDir, "Data/Cleaned/Audrie/pregnancy_child_immune_covariates_data.RDS"))

Xvars <- c("vitD_nmol_per_L", "logFERR_inf", "logSTFR_inf", "logRBP_inf", 
           "vit_A_def", "iron_def", "vit_D_def", "ln_preg_cort", "ln_preg_estri",
           "logCRP", "logAGP", "mom_t0_ln_ifn", "sumscore_t0_mom_Z") 

filtering <- function(row){
  any(!is.na(row))
}
has_exp<-d[apply(select(d, Xvars), 1, filtering),]
has_outy1<-d[apply(select(d, grep("t2_ln", names(d)), sumscore_t2_Z), 1, filtering),]
has_outy2<-d[apply(select(d, grep("t3_ln", names(d)), sumscore_t3_Z), 1, filtering),]

has_exp$dataid%>%unique()%>%length()
has_outy1$childid%>%unique()%>%length()
has_outy2$childid%>%unique()%>%length()

y1 <- has_exp[apply(select(has_exp, grep("t2_ln", names(d)), sumscore_t2_Z), 1, filtering),] 
y2 <- has_exp[apply(select(has_exp, c("t3_ln_crp", "t3_ln_agp", "t3_ln_ifn","sumscore_t3_Z"), sumscore_t3_Z), 1, filtering),]

nrow(y1)
y1$dataid %>% unique() %>% length()
y1$childid %>% unique() %>% length()

nrow(y2)
y2$dataid %>% unique() %>% length()
y2$childid %>% unique() %>% length()