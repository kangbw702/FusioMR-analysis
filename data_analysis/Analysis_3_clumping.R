# clumping using R/4.2.1

library(data.table)
library(tidyverse)

source("/gpfs/data/linchen-lab/Bowei/fusiomr_sc_bk/r1/code/functions.R")

pcut1 = 1e-4; pcut2 = 1e-3; clump_kb = 100; clump_r2 = 0.01

# LDL on ischemic stroke, EUR SAS
df = fread("harmo_tbl_ldl_ais_eur_sas.csv") # # 2,809,918 rows
df = df %>% mutate(across(exp1.beta:out2.p, ~ as.numeric(.)))
df1 = df[df$exp1.p < pcut1 & df$exp2.p < pcut2, ] # take intersection
nrow(df1) # 2163

# clumping
bfile_eur = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR"
bfile_sas = "/gpfs/data/linchen-lab/Bowei/analysis_3/1kg/SAS"
clumped_df = clump(df1, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = bfile_sas, pval_col = "exp1.p", plink_bin = genetics.binaRies::get_plink_binary())
is.data.table(clumped_df)
nrow(clumped_df) # 131

fwrite(clumped_df, "dfclump_ldrefsas_ldl_ais_eur_sas.csv")







