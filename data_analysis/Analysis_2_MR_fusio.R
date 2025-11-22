library(data.table)
library(tidyverse)
library(pbapply)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)

source("functions/utilites.R")
source("functions/set_var_prior.R")
source("functions/set_init_seso.R")
source("functions/set_init_semo.R")
Rcpp::sourceCpp("fusiomr/fusiomr_s_uhp_only.cpp")
Rcpp::sourceCpp("fusiomr/fusiomr_m_shared_exp.cpp")

array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
chr = array_id
print(chr)

pcut = 0.001; clump_kb = 50; clump_r2 = 0.01
glob_q = 0.9; kappa_hybrid = 50; c_gamma = 4; c_theta = c_gamma*1.5; kappa_gamma_s = 1; kappa_theta_s = 1.5; kappa_gamma_m = 1; kappa_theta_m = 1.25; K_cut = 3
hybrid = TRUE
niv_cut = 5; out_dir = "res_r1"

traits = c("stroke", "af") 
tissue_short = c("Amygdala","BA24","Caudate_basal_ganglia","Cerebellar_Hemi","Cerebellum","Cortex","BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord","Substantia_nigra","Blood","Aorta","Coronary","Left_Ventricle","Atrial_Appendage")

for (i in 1:length(tissue_short)) {
tissue = tissue_short[i]

df_clump_list = readRDS(paste0('df_clump/dfclump_', pcut, '_', clump_kb, '_', clump_r2, '_', tissue, '_', chr, '.rds')) 
gene_list = names(df_clump_list)
# get loc and glob mom tables from step 2
mtbls = fread(paste0("mtbls/mtbls_apa_", glob_q, "_", pcut, "_", K_cut, "_1_2_", clump_kb, "_", clump_r2, "_", tissue, ".txt"))
gene_list = intersect(gene_list, unique(mtbls$gene_symbol))
cat(tissue, "chr", chr, ",", length(gene_list), "genes are tested.\n")

for (j in 1:length(gene_list)) {
gg = gene_list[j]
df = df_clump_list[[gg]]
b_exp = df$beta; se_exp = df$se
b_out_1 = df$b_out_1; se_out_1 = df$se_out_1
b_out_2 = df$b_out_2; se_out_2 = df$se_out_2
K_post = nrow(df)
cat("Processing gene", j, gg, K_post, "IVs to process. Total genes:", length(gene_list), tissue, "\n")

# global
global_m_gamma = mtbls$global_m_gamma[mtbls$gene_symbol == gg]
global_m_theta11 = mtbls$global_m_theta11[mtbls$gene_symbol == gg]
global_m_theta22 = mtbls$global_m_theta22[mtbls$gene_symbol == gg]
global_m_theta12 = mtbls$global_m_theta12[mtbls$gene_symbol == gg]
global_m_theta_mat = matrix(c(global_m_theta11, global_m_theta12, global_m_theta12, global_m_theta22), 2, 2)

if (K_post >= niv_cut) {
# run MR
# FusioMRm semo, 1-stroke, 2-af
niter = 20000
# hyper parameters
# parameters for variance priors
vp = set_variance_priors_m2(ghat = b_exp, gse = se_exp, Ghat_mat = matrix(c(b_out_1, b_out_2), ncol=2), Gse_mat = matrix(c(se_out_1, se_out_2), ncol=2), beta0 = NULL, K = K_post, Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = global_m_gamma, global_mean_theta = global_m_theta_mat, hybrid = hybrid, kappa_hybrid = kappa_hybrid, z_thresh = NULL, trim = 0.0, kappa_gamma = kappa_gamma_m, kappa_theta = kappa_theta_m)
a_gamma_prior = vp$gamma$a
b_gamma_prior = vp$gamma$b
nu_theta_prior = vp$theta$nu
Phi_theta_prior = vp$theta$Phi
# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup_semo_uhp_only(niter, K_post, beta_1_init = vp$beta0[1], beta_2_init = vp$beta0[2], sigma_gamma_init = sqrt(vp$gamma$prior_mean))
# rcpp
res_semo = gibbs_semo_uhp_only_rcpp(niter, K_post, start_val$beta_1_tk, start_val$beta_2_tk, start_val$theta_1_tk, start_val$theta_2_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, b_out_1, b_out_2, (se_out_1)^2, (se_out_2)^2, b_exp, (se_exp)^2, a_gamma_prior, b_gamma_prior, vp$theta$prior_mean, nu_theta_prior, Phi_theta_prior)
ids = (niter/2 + 1):niter
b1_semo = mean(res_semo$beta_1_tk[ids], na.rm = T)
se1_semo = sd(res_semo$beta_1_tk[ids], na.rm = T)
p1_semo = 2*exp(pnorm(abs(b1_semo)/se1_semo, lower.tail = FALSE, log.p = TRUE))  

b2_semo = mean(res_semo$beta_2_tk[ids], na.rm = T)
se2_semo = sd(res_semo$beta_2_tk[ids], na.rm = T)
p2_semo = 2*exp(pnorm(abs(b2_semo)/se2_semo, lower.tail = FALSE, log.p = TRUE)) 

#bci_1_semo = quantile(res_semo$beta_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
#bci_2_semo = quantile(res_semo$beta_2_tk[ids], probs=c(0.025,0.975), na.rm=T)

# Fusio-s, outcome2 AF
niter = 20000
# hyper parameters
# parameters for variance priors
vp = set_variance_priors(ghat = b_exp, gse = se_exp, Ghat = b_out_2, Gse = se_out_2, beta0 = NULL, K = K_post, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = global_m_gamma, global_mean_theta = global_m_theta22, hybrid = hybrid, kappa_hybrid = kappa_hybrid, z_thresh = NULL, trim = 0.0, kappa_gamma = kappa_gamma_s, kappa_theta = kappa_theta_s)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
# starting values. Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K_post, alpha_init = 1, beta_init = vp$beta0, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean))
# rcpp
res_s2 = gibbs_seso_uhp_only_cpp(niter, K_post, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out_2, b_exp, (se_out_2)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
ids = (niter/2 + 1):niter
bhat_s2 = mean(res_s2$beta_tk[ids], na.rm = T)
se_bhat_s2 = sd(res_s2$beta_tk[ids], na.rm = T)
p_s2 = 2*exp(pnorm(abs(bhat_s2)/se_bhat_s2, lower.tail = FALSE, log.p = TRUE)) 
# 1-coverage, bCI
#bci_s2 = quantile(res_s2$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)

# write out
tmp = c(tissue, chr, gg, K_post, b1_semo, se1_semo, p1_semo, b2_semo, se2_semo, p2_semo, bhat_s2, se_bhat_s2, p_s2)
write.table(t(tmp), paste0(out_dir, '/fusio_res_', tissue, '_s7.txt'), quote=F, row.names=F, col.names=F, append=T)
} # end if
} # end loop j
} # end loop i

cat("step3, run FusioMR methods, finished! Chr", chr, "\n")






