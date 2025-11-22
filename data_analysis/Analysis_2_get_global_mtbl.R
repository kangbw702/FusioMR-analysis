# calculate global mom R/4.2.1

library(data.table)
library(tidyverse)
library(pbapply)

source("functions/set_var_prior.R")

#=====================================================================================================================

# Build a global EB center from many genes' local MoMs (vector),
# using a trimmed median of positive values (robust, no recomputation inside priors).
build_global_mom <- function(local_mom_vec, q = 0.5, trim = 0.05, min_pos = 20) {
x = local_mom_vec[is.finite(local_mom_vec) & local_mom_vec > 0]
if (length(x) == 0) return(0)
if (length(x) < min_pos) return(stats::median(x, na.rm = TRUE))
lo = stats::quantile(x, trim, na.rm = TRUE)
hi = stats::quantile(x, 1 - trim, na.rm = TRUE)
x = x[x >= lo & x <= hi]
if (!length(x)) return(0)
stats::quantile(x, q, na.rm = TRUE)
}

#=====================================================================================================================

# aggregate scale (variances) and dependence (correlation) separately using robust summaries; then recompose. This is very stable for small K_g and outliers.

build_global_mom_covmat <- function(mat, q = 0.5, eps = 1e-8) {
# mat: each row is a gene, col1,col2,col3 = s11,s22,s12 in 2x2 PSD matrices
A = mat[,1]
B = mat[,2]
C = mat[,3]
# guard small/negative due to numerical noise
A = pmax(A, eps); B = pmax(B, eps)
r = C / sqrt(A*B); r = pmin(pmax(r, -1+1e-10), 1-1e-10)
v1_bar = stats::quantile(log(A), q, na.rm = TRUE)
v2_bar = stats::quantile(log(B), q, na.rm = TRUE)
z_bar = stats::quantile(atanh(r), 0.5, na.rm = TRUE) 
r_bar = tanh(z_bar)
s1 = exp(v1_bar); s2 = exp(v2_bar)
matrix(c(s1, r_bar*sqrt(s1*s2), r_bar*sqrt(s1*s2), s2), 2, 2)
}

#=====================================================================================================================


pcut = 0.001; clump_kb = 50; clump_r2 = 0.01
glob_q = 0.9; kappa_gamma = 1; kappa_theta = 2
K_cut = 3

traits = c("stroke", "af") 
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")

#for (i in 1:length(tissue_short)) {
#for (i in 13:length(tissue_short)) {
i=12
tissue = tissue_short[i]
		
# calculate global moms
res = NULL
#for (chr in 1:22) {
for (chr in c(1:9,11:22)) {
cat("tissue", i, tissue, "chr", chr, "\n")

df_clump_list = readRDS(paste0('df_clump/dfclump_', pcut, '_', clump_kb, '_', clump_r2, '_', tissue, '_', chr, '.rds')) 

tbl = do.call(rbind, pblapply(df_clump_list, function(df) {
ghat = df$beta; gse = df$se
Ghat_mat = matrix(c(df$b_out_1, df$b_out_2), ncol=2)
Gse_mat = matrix(c(df$se_out_1, df$se_out_2), ncol=2)
moms = set_variance_priors_m2(ghat, gse, Ghat_mat, Gse_mat, beta0 = NULL, K = NULL, Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0, c_gamma = 0.5, c_theta = 0.8, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.10, kappa_gamma = kappa_gamma, kappa_theta = kappa_theta)
c(mom_gamma = moms$gamma$prior_mean, mom_theta11 = moms$theta$prior_mean[1,1], mom_theta22 = moms$theta$prior_mean[2,2], mom_theta12 = moms$theta$prior_mean[1,2], beta0 = moms$beta0, K = moms$K)
}))

tbl = as.data.frame(tbl)
tbl$gene_symbol = names(df_clump_list)
tbl$tissue = tissue
tbl$chr = chr
rownames(tbl) = NULL
res = bind_rows(res, tbl) 
} # end loop chr

# Globals per cell type using the same aux pipeline:
# trimmed median of positive locals
res_g = res[res$K >= K_cut,]
global_mean_gamma = build_global_mom(res_g$mom_gamma, q = glob_q, trim = 0.05, min_pos = 20)  
global_mean_theta = build_global_mom_covmat(as.matrix(res_g[,2:4]), q = glob_q, eps = 1e-8)
res$global_K = nrow(res_g)
res$global_m_gamma = global_mean_gamma
res$global_m_theta11 = global_mean_theta[1,1]
res$global_m_theta22 = global_mean_theta[2,2]
res$global_m_theta12 = global_mean_theta[1,2]

fwrite(res, paste0("mtbls/mtbls_apa_", glob_q, "_", pcut, "_", K_cut, "_", kappa_gamma, "_", kappa_theta, "_", clump_kb, "_", clump_r2, "_", tissue, ".txt"))
#} # end loop tissue

cat("step2, calculating local and global mom from clump df, finished!\n")



#=====================================================================================================================
# UPDATE global moms. mtbls ALREADY exists

traits = c("stroke", "af") 
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Tibial","Left_Ventricle","Atrial_Appendage")

#for (i in 1:length(tissue_short)) {
i=6
tissue = tissue_short[i]
cat("tissue\n")

# read in existing mtbl
pcut = 0.001; clump_kb = 50; clump_r2 = 0.01; glob_q = 0.75; kappa_gamma = 1; kappa_theta = 1; K_cut = 3
res = fread(paste0("mtbls/mtbls_apa_", glob_q, "_", pcut, "_", K_cut, "_", kappa_gamma, "_", kappa_theta, "_", clump_kb, "_", clump_r2, "_", tissue, ".txt"))

# Globals per cell type using the same aux pipeline:
# trimmed median of positive locals
res_g = res[res$K >= K_cut,]
global_mean_gamma = build_global_mom(res_g$mom_gamma, q = glob_q, trim = 0.05, min_pos = 20)  
global_mean_theta = build_global_mom_covmat(as.matrix(res_g[,2:4]), q = glob_q)
res$global_K = nrow(res_g)
res$global_m_gamma = global_mean_gamma
res$global_m_theta11 = global_mean_theta[1,1]
res$global_m_theta22 = global_mean_theta[2,2]
res$global_m_theta12 = global_mean_theta[1,2]

fwrite(res, paste0("mtbls/mtbls_apa_", glob_q, "_", pcut, "_", K_cut, "_", kappa_gamma, "_", kappa_theta, "_", clump_kb, "_", clump_r2, "_", tissue, ".txt"))
#} # end loop tissue

cat("step2-u, updating global mom, finished!\n")






