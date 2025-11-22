# run MR using R/4.5.1

library(data.table)
library(tidyverse)
library(MendelianRandomization)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)

source("code/functions.R")
source("functions/set_init_seso.R")
source("functions/set_var_prior.R")
Rcpp::sourceCpp("fusiomr/fusiomr_s_uhp_only.cpp")

# function to run MR, seso
run_mr = function(gene, b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior, sigma_gamma_init, sigma_theta_init) {
# FusioMRs, uhp only
niter = 20000
K = length(b_exp)
start_val = init_setup(niter, K, alpha_init = 1, beta_init = 0, sigma_gamma_init = sigma_gamma_init, sigma_theta_init = sigma_theta_init)
res = gibbs_seso_uhp_only_cpp(
niter, K, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk,
start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk,
b_out, b_exp, (se_out)^2, (se_exp)^2,
a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior
)
if(all(!is.na(res))) {
ids = (niter/2 + 1):niter
bhat = mean(res$beta_tk[ids], na.rm = T)
se_bhat = sd(res$beta_tk[ids], na.rm = T)
pval = 2*(1-pnorm(abs(bhat)/se_bhat))
} else { 
bhat = se_bhat = pval = NA 
}

# competing methods
mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
b_ivw_fixed = IVW_f$Estimate; se_ivw_fixed = IVW_f$StdError; p_ivw_fixed = IVW_f$Pvalue
Egger = tryCatch({MendelianRandomization::mr_egger(mr.obj)}, error = function(e) {NA})
if(!is.null(Egger) & class(Egger)=='Egger') { 
b_egger = Egger$Estimate; se_egger = Egger$StdError.Est; p_egger = Egger$Pvalue.Est
} else b_egger = se_egger = p_egger = NA
cml = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = FALSE, n = 190)}, error = function(e) {NA})
if(!is.null(cml)) {  
b_cml = cml$Estimate; se_cml = cml$StdError; p_cml = cml$Pvalue
} else b_cml = se_cml = p_cml = NA
cml_dp = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = TRUE, num_pert = 200, n = 190)}, error = function(e) {NA})
if(!is.null(cml_dp)) {  
b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError; p_cml_dp = cml_dp$Pvalue
} else b_cml_dp = se_cml_dp = p_cml_dp = NA

# out
out_other_methods = c(length(b_exp), b_ivw_fixed, se_ivw_fixed, p_ivw_fixed, b_egger, se_egger, p_egger, b_cml, se_cml, p_cml, b_cml_dp, se_cml_dp, p_cml_dp)
output = c(gene, out_other_methods, c(bhat, se_bhat, pval), sigma_gamma_init^2, sigma_theta_init^2)
names(output) = c("gene", "niv", "b_ivw", "se_ivw", "p_ivw", "b_egger", "se_egger", "p_egger", "b_cml", "se_cml", "p_cml", "b_cml_dp", "se_cml_dp", "p_cml_dp", "b_fusio", "se_fusio", "p_fusio", "sigma2_gamma_prior_mean", "sigma2_theta_prior_mean")
return(output)
}

args <- commandArgs(trailingOnly = T)
chr = as.numeric(args[1])
print(chr)

pcut = 0.001; clump_kb = 50; clump_r2 = 0.1; glob_q = 0.9
floor_frac = 0.6; kappa = 50; c_gamma = 1; c_theta = c_gamma*1.5
ctypes_short = c("end","ast","exc","inh","mic","opc","oli","per")
diseases_short = "ad"

for (i in 1:length(ctypes_short)) {
print(paste0("ctype-", i, ", ", ctypes_short[i]))
df_clump_list = readRDS(paste0("df_clump/dfclump_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], "_", chr, '.rds'))
gene_list = names(df_clump_list)

# get loc and glob mom tables 
mtbls = fread(paste0("mtbls/mtbls_testclump_", diseases_short, "_", glob_q, "_", pcut, "_", clump_kb, "_", clump_r2, "_", ctypes_short[i], ".txt")) 
gene_list = intersect(gene_list, unique(mtbls$gene_symbol))
cat(length(gene_list), "genes are tested.\n")

sumtbl1 = NULL
sumtbl1 = do.call(rbind, lapply(seq_along(gene_list), function(j) {
gg = gene_list[j]
b_exp = df_clump_list[[gg]]$beta; se_exp = df_clump_list[[gg]]$se; b_out = df_clump_list[[gg]]$b_out; se_out = df_clump_list[[gg]]$se_out
K_post = nrow(df_clump_list[[gg]])
cat("Processing gene", j, gg, K_post, "IVs to process\n")
# set variance priors
# Globals per cell type (compute once) using the same aux pipeline:
local_m_gamma  = mtbls$mom_gamma[mtbls$gene_symbol == gg]
local_m_theta  = mtbls$mom_theta[mtbls$gene_symbol == gg]
global_m_gamma = mtbls$global_m_gamma[mtbls$gene_symbol == gg]
global_m_theta = mtbls$global_m_theta[mtbls$gene_symbol == gg]
# flooring and blending
lgamma = max(local_m_gamma, floor_frac * global_m_gamma)
ltheta = max(local_m_theta, floor_frac * global_m_theta)
eta = K_post / (K_post + kappa)
m_gamma = eta * lgamma + (1 - eta) * global_m_gamma
m_theta = eta * ltheta + (1 - eta) * global_m_theta
# get a,b for Inv-gamma
a_gamma_prior = max(1.1, c_gamma * pmax(5, pmin(K_post, 20)) / 2)
a_theta_prior = max(1.1, c_theta * pmax(5, pmin(K_post, 20)) / 2)
b_gamma_prior = (a_gamma_prior - 1) * m_gamma
b_theta_prior = (a_theta_prior - 1) * m_theta
# initial value for variances
sigma_gamma_init = sqrt(m_gamma)
sigma_theta_init = sqrt(m_theta)
res = run_mr(gg, b_exp, se_exp, b_out, se_out, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior, sigma_gamma_init, sigma_theta_init)
return(res)
}
))

fwrite(sumtbl1, paste0("sumtbl_bryois/sumtbl_", diseases_short, "_", pcut_aux, "_", pcut, "_", clump_kb, "_", clump_r2, "_", glob_q, "_", floor_frac, "_", kappa, "_", c_gamma, "_", ctypes_short[i], "_", chr, '.txt'))
} # end loop cell type

cat("step4, run MR finished!\n")



