# clumping using R/4.2.1

library(data.table)
library(tidyverse)

source("functions/utilities.R")

# EUR+SAS
pcut1 = 1e-4; pcut2 = 1e-3; clump_kb = 100; clump_r2 = 0.01
df = fread("harmo_tbl_ldl_ais_eur_sas.csv") 
df = df %>% mutate(across(exp1.beta:out2.p, ~ as.numeric(.)))
df1 = df[df$exp1.p < pcut1 & df$exp2.p < pcut2, ] 
# clumping 
clumped_df = clump(df1, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "1kg/SAS", pval_col = "exp1.p", plink_bin = genetics.binaRies::get_plink_binary())
is.data.table(clumped_df)
nrow(clumped_df)
fwrite(clumped_df, "dfclump_ldrefsas_ldl_ais_eur_sas.csv")

# SAS
pcut = 1e-3; clump_kb = 100; clump_r2 = 0.01
df = fread("harmo_tbl_ldl_ais_eur_sas.csv") 
df = df %>% mutate(across(exp1.beta:out2.p, ~ as.numeric(.)))
df1 = df[df$exp2.p < pcut, ] 
# clumping
clumped_df = clump(df1, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "1kg/SAS", pval_col = "exp2.p", plink_bin = genetics.binaRies::get_plink_binary())
is.data.table(clumped_df)
nrow(clumped_df)
fwrite(clumped_df, paste0("dfclump_ldl_ais_sas_", pcut, ".csv"))









