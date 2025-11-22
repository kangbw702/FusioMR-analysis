# R/4.5.1

library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(Rcpp)
library(MendelianRandomization)

source("functions/set_var_prior.R")
source("functions/set_init_seso.R")
source("functions/utilities.R")
Rcpp::sourceCpp("fusiomr/fusiomr_s_with_chp.cpp")


# LDL on ischemic stroke, SAS
pcut = 1e-3
dat = fread(paste0("dfclump_ldl_ais_sas_", pcut, ".csv")) 

# get summary statistics 
b_exp_2 = dat$exp2.beta; se_exp_2 = dat$exp2.se; 
b_out_2 = dat$out2.beta; se_out_2 = dat$out2.se;
K = length(b_exp_2)

# FusioMRs empirical (outcome 2, SAS)
niter = 20000
# parameters for variance priors
vp = set_variance_priors(ghat = b_exp_2, gse = se_exp_2, Ghat = b_out_2, Gse = se_out_2, beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = 1, c_theta = 1.5, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = 1, kappa_theta = 1)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
a_q_prior = b_q_prior = 1
# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K, alpha_init = 1, beta_init = vp$beta0, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean), q_init = 0)
# rcpp
res2 = gibbs_rcpp_nopr(niter, K, beta_tk = start_val$beta_tk, alpha_tk = start_val$alpha_tk, eta_tk = start_val$eta_tk, theta_tk = start_val$theta_tk, gamma_tk = start_val$gamma_tk, q_tk = start_val$q_tk, b_out_2, b_exp_2, (se_out_2)^2, (se_exp_2)^2, sigma2_gamma_tk = start_val$sigma2_gamma_tk, sigma2_theta_tk = start_val$sigma2_theta_tk, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior, a_q_prior, b_q_prior) 
ids = (niter/2 + 1):niter
post_flip_s2 = label_flip(niter, res2, rep(0, K))
b_s2 = post_flip_s2$b_mean; se_s2 = post_flip_s2$b_sd; p_s2 = 2*exp(pnorm(abs(b_s2/se_s2), lower.tail=FALSE, log.p=TRUE))
bci_s2 = post_flip_s2$bci

set.seed(1)
# create mr object
mr.obj = MendelianRandomization::mr_input(bx = b_exp_2, bxse = se_exp_2, by = b_out_2, byse = se_out_2)
# IVW fixed
IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
b_ivw = IVW_f$Estimate; se_ivw = IVW_f$StdError; p_ivw = IVW_f$Pvalue
# MR-Egger
Egger = tryCatch({MendelianRandomization::mr_egger(mr.obj)}, error = function(e) {NA})
b_egger = Egger$Estimate; se_egger = Egger$StdError.Est; p_egger = Egger$Pvalue.Est
# cML
cml = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = FALSE, n = 11000)}, error = function(e) {NA})
b_cml = cml$Estimate; se_cml = cml$StdError; p_cml = cml$Pvalue
# n_iiv_c = sum(cml_c$BIC_invalid)
# cML-DP
cml_dp = tryCatch({MendelianRandomization::mr_cML(mr.obj, MA = TRUE, DP = TRUE, n = 11000, num_pert = 200)}, error = function(e) {NA})
b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError; p_cml_dp = cml_dp$Pvalue

# summarize
res = c(K = K, b_fusios = b_s2, se_fusios = se_s2, p_fusios = p_s2, b_ivw, se_ivw, p_ivw, b_egger, se_egger, p_egger, b_cml, se_cml, p_cml, b_cml_dp, se_cml_dp, p_cml_dp)
res







