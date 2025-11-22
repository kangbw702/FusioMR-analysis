# clumping using R/4.2.1

library(data.table)
library(tidyverse)
library(pbapply)

# clumping function
clump <- function(dat, SNP_col = "eQTL_variant_id", pval_col = "rowmeta", clump_kb = 250, clump_r2 = 0.1, clump_p = 0.999, bfile = "/scratch/t.phs.yihaolu/GWAS_summary_statistics/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01", plink_bin = "plink", pop="EUR") {
df <- data.frame(rsid = dat[, ..SNP_col], pval = dat[,..pval_col])
colnames(df) = c("rsid", "pval")
out <- tryCatch({
ieugwasr::ld_clump(df, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin = plink_bin, pop = pop)
}, silent = TRUE, error = function(x) return(NA)
)
if(length(out)==1) return(NA)
MRdat <- dat[which(unlist(dat[,..SNP_col]) %in% out$rsid),]
return(MRdat)
}

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

pcut = 0.001; clump_kb = 50; clump_r2 = 0.01
traits = c("stroke", "af") 
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")

for (i in 1:length(tissue_short)) {
tissue = tissue_short[i]
		
res = NULL
cat("tissue", i, tissue, "chr", chr, "\n")
dat = fread(paste0('harmo/apaqtl_', traits[1], '_', traits[2], '_', tissue, '_', chr, '.txt.gz')) 
dat1 = dat %>% dplyr::filter(se != Inf, se_out_1 != Inf, se_out_2 != Inf, se > 0, se_out_1 > 0, se_out_2 > 0) %>% drop_na()
dat1 = dat1 %>% group_by(gene_symbol) %>% distinct(rsid, .keep_all = TRUE) %>% ungroup()
gene_list1 = unique(dat1$gene_symbol)
cat(length(gene_list1), "genes are available in", tissue, "chr", chr, "\n")

df_clump_list = pblapply(seq_along(gene_list1), function(j) {
message(paste0("Processing gene-", j, ": ", gene_list1[j]))
dat1 = dat1 %>% dplyr::filter(gene_symbol == gene_list1[j]) 
tmp = dat1 %>% dplyr::filter(p <= pcut)
class(tmp) = 'data.table'
if (nrow(tmp) > 0) {
clumped_df = clump(tmp, SNP_col = "rsid", clump_kb = clump_kb, clump_r2 = clump_r2, bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", pval_col = "p", plink_bin = genetics.binaRies::get_plink_binary())
if (class(clumped_df) %in% c("data.table", "data.frame") && nrow(clumped_df) >= 1 && nrow(clumped_df) < 1000) {
return(clumped_df)
}}
NULL
})

names(df_clump_list) = gene_list1
df_clump_list = base::Filter(Negate(is.null), df_clump_list)
cat(length(df_clump_list), "genes are left after clumping, nIVs in [1, 1000).\n")

saveRDS(df_clump_list, paste0("df_clump/dfclump_", pcut, "_", clump_kb, "_", clump_r2, "_", tissue, "_", chr, '.rds'))
} # end loop cell type

cat(length(gene_list1), "genes are available in", tissue, "chr", chr, "\n")
cat(length(df_clump_list), "genes are left after clumping, nIVs in [1, 1000).\n")
cat("step1 clumping finished!\n")







