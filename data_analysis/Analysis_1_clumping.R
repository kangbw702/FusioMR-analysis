# clumping using R/4.2.1

library(data.table)
library(tidyverse)
library(pbapply)

source("code/functions.R")

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

pcut = 0.001
clump_kb = 50; clump_r2 = 0.1
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
diseases_short = "ad"

for (i in 1:length(ctypes_short)) {
print(paste0("ctype-", i, ", ", ctypes_short[i]))
dat = fread(paste0('../harmo_out/eqtl_', diseases_short, '_', ctypes_short[i], "_", chr, '.txt.gz')) 
dat1 = dat %>% dplyr::select(gene_symbol:se, b_out, se_out) %>% dplyr::filter(se != Inf, se_out != Inf, se > 0, se_out > 0) %>% drop_na() 
gene_list1 = unique(dat1$gene_symbol)
cat(length(gene_list1), "genes are available in", ctypes_short[i], "chr", chr, "\n")
df_clump_list = pblapply(seq_along(gene_list1), function(j) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
dat1 = dat1 %>% dplyr::filter(gene_symbol == gene_list1[j]) 
tmp = dat1 %>% dplyr::filter(p <= pcut)
if (nrow(tmp) > 0) {
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
if (all(class(clumped_df) == c("data.table", "data.frame")) && nrow(clumped_df) >= 1 && nrow(clumped_df) < 1000) {
return(clumped_df)
}}
NULL
})

names(df_clump_list) = gene_list1
df_clump_list = base::Filter(Negate(is.null), df_clump_list)
cat(length(df_clump_list), "genes are left after clumping, nIVs in [1, 1000).\n")

saveRDS(df_clump_list, paste0("df_clump/dfclump_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], "_", chr, '.rds'))
} # end loop cell type

cat("step1 clumping for testing finished!\n")






