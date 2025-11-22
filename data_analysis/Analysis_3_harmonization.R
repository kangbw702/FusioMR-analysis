# using R/4.2.1

library(data.table)
library(tidyverse)

source("functions/utilities.R")

# LDL on ischemic stroke, EUR SAS
# harmonize outcomes
# load data
out1 = fread("sum_stat/GCST90104540_ais_eur.h.tsv.gz", select = c("chromosome", "rsid", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")) %>% 
dplyr::rename(chr = chromosome, REF = other_allele, ALT = effect_allele, out1.af = effect_allele_frequency, out1.beta = beta, out1.se = standard_error, out1.p = p_value) %>% 
drop_na() %>% distinct(rsid, .keep_all = TRUE) %>% drop_na() 

out2 = fread("sum_stat/GCST90104560_ais_sas.h.tsv.gz", select = c("chromosome", "rsid", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value")) %>% 
dplyr::rename(chr = chromosome, REF = other_allele, ALT = effect_allele, out2.af = effect_allele_frequency, out2.beta = beta, out2.se = standard_error, out2.p = p_value) %>% 
drop_na() %>% distinct(rsid, .keep_all = TRUE) %>% drop_na() 

# harmonize two outcome gwas
df_out = out1 %>% left_join(out2, by = c('rsid', 'chr')) %>% 
mutate(align = ifelse(ALT.x==ALT.y & REF.x==REF.y, 1, ifelse(ALT.x==REF.y & REF.x==ALT.y, -1, NA))) %>% 
dplyr::filter(!(is.na(align))) %>% mutate(out2.beta = out2.beta*align, out2.af = ifelse(align==1, out2.af, 1-out2.af)) %>% 
dplyr::select(-REF.y, -ALT.y, -align) %>% dplyr::rename(EFF=ALT.x, REF=REF.x) 

# load exposure data
exp1 = fread("sum_stat/GCST90239658_ldl_eur.h.tsv.gz", select = c("chromosome", "rsid", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "effect_allele_frequency")) %>% 
dplyr::rename(chr = chromosome, REF = other_allele, ALT = effect_allele, exp1.beta = beta, exp1.se = standard_error, exp1.p = p_value, exp1.af = effect_allele_frequency) %>% 
drop_na() %>% dplyr::filter(nchar(REF) == 1, nchar(ALT) == 1) %>% distinct(rsid, .keep_all = TRUE) 

exp2 = fread("sum_stat/GCST90239660_ldl_sas.h.tsv.gz", select = c("chromosome", "rsid", "effect_allele", "other_allele", "beta", "standard_error", "p_value", "effect_allele_frequency")) %>% 
dplyr::rename(chr = chromosome, REF = other_allele, ALT = effect_allele, exp2.beta = beta, exp2.se = standard_error, exp2.p = p_value, exp2.af = effect_allele_frequency) %>% 
drop_na() %>% dplyr::filter(nchar(REF) == 1, nchar(ALT) == 1) %>% distinct(rsid, .keep_all = TRUE) 

# harmonize two exposure gwas
df_exp = exp1 %>% left_join(exp2, by = c('rsid', 'chr')) %>% mutate(align = ifelse(ALT.x==ALT.y & REF.x==REF.y, 1, ifelse(ALT.x==REF.y & REF.x==ALT.y, -1, NA))) %>% 
dplyr::filter(!(is.na(align))) %>% mutate(exp2.beta = exp2.beta*align, exp2.af = ifelse(align==1, exp2.af, 1-exp2.af)) %>% 
dplyr::select(-REF.y, -ALT.y, -align) %>% dplyr::rename(EFF=ALT.x, REF=REF.x) 

# harmonize exposures + outcomes, df_exp + df_out
df = df_exp %>% left_join(df_out, by = c('rsid', 'chr')) %>% mutate(align = ifelse(EFF.x==EFF.y & REF.x==REF.y, 1, ifelse(EFF.x==REF.y & REF.x==EFF.y, -1, NA))) %>% 
dplyr::filter(!(is.na(align))) %>% mutate(out1.beta = out1.beta*align, out2.beta = out2.beta*align, out1.af = ifelse(align==1, out1.af, 1-out1.af), out2.af = ifelse(align==1, out2.af, 1-out2.af)) %>% dplyr::select(-REF.y, -EFF.y, -align) %>% 
dplyr::rename(EFF=EFF.x, REF=REF.x) 

# save intermediate table
fwrite(df, "harmo_tbl_ldl_ais_eur_sas.csv")







