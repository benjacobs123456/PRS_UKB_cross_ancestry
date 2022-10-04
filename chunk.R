library(tidyverse)
setwd("/data/home/hmy117/PRS_UKB_cross_ancestry/")
args = commandArgs(trailingOnly=TRUE)
pheno = read_tsv("ms_pheno.tsv")

ms = pheno %>% filter(MS_status==1) %>% sample_n(size=42,replace = F)
cont = pheno %>% filter(MS_status==0) %>% sample_n(size=40490,replace = F)
# sample
pheno_sample = bind_rows(ms,cont)
write_tsv(pheno_sample,path = paste0("prs_subset_iter_",args[1],".tsv"))

