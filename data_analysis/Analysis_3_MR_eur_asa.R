# R/4.5.1

library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(Rcpp)
library(MR2)

source("functions/set_var_prior.R")
source("functions/utilities.R")
source("functions/set_init_memo.R")
Rcpp::sourceCpp("fusiomr/fusiomr_m_multi_exp.cpp")

# LDL on ischemic stroke, EUR SAS
# pcut1 = 1e-4; pcut2 = 1e-3; clump_kb = 100; clump_r2 = 0.01
dat = fread("dfclump_ldrefsas_ldl_ais_eur_sas.csv") 

# get summary statistics 
b_exp_1 = dat$exp1.beta; se_exp_1 = dat$exp1.se;
b_exp_2 = dat$exp2.beta; se_exp_2 = dat$exp2.se; 
b_out_1 = dat$out1.beta; se_out_1 = dat$out1.se;
b_out_2 = dat$out2.beta; se_out_2 = dat$out2.se;
K = length(b_exp_2)

# FusioMRm
niter = 20000
# hyper parameters
# parameters for variance priors
vp = set_variance_priors_m2x2_diag(ghat_mat = matrix(c(b_exp_1, b_exp_2), ncol=2), gse_mat = matrix(c(se_exp_1, se_exp_2), ncol=2), Ghat_mat = matrix(c(b_out_1, b_out_2), ncol=2), Gse_mat = matrix(c(se_out_1, se_out_2), ncol=2), B0 = NULL, K = K, Kmin = 5, Kmax = 20, rho12 = 0, rho_gg = 0, rho_gj = list(c(0,0), c(0,0)), c_gamma = 1, c_theta = 1.5, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, kappa_gamma = 1, kappa_theta = 1)
nu_gamma_prior = vp$gamma$nu; Phi_gamma_prior = vp$gamma$Phi
nu_theta_prior = vp$theta$nu; Phi_theta_prior = vp$theta$Phi
pst_init = bern_prob(0.5, 0.1, 0.1)
cc = rep(1, 4) # Dirichlet parameters
# starting values
# Use prior mean as init for Sigma_gamma and Sigma_theta
start_val = init_setup_memo(niter, K, alpha_1_init = 0.1, alpha_2_init = 0.1, beta_1_init = vp$B0[1], beta_2_init =vp$B0[2], eta_1_init = rep(0,K), eta_2_init = rep(0,K), pst_init = pst_init) 
res_memo = gibbs_joint_rcpp_nopr(niter, K, beta_1_tk = start_val$beta_1_tk, beta_2_tk = start_val$beta_2_tk, alpha_1_tk = start_val$alpha_1_tk, alpha_2_tk = start_val$alpha_2_tk, eta_1_tk = start_val$eta_1_tk, eta_2_tk = start_val$eta_2_tk, theta_1_tk = start_val$theta_1_tk, theta_2_tk = start_val$theta_2_tk, gamma_1_tk = start_val$gamma_1_tk, gamma_2_tk = start_val$gamma_2_tk, pst_tk = start_val$pst_tk, b_out_1, b_exp_1, (se_out_1)^2, (se_exp_1)^2, b_out_2, b_exp_2, (se_out_2)^2, (se_exp_2)^2, vp$gamma$prior_mean, vp$theta$prior_mean, nu_gamma_prior, Phi_gamma_prior, nu_theta_prior, Phi_theta_prior, cc) 
ids = (niter/2 + 1):niter
post_flip = label_flip_joint(niter, res_memo, rep(0,K), rep(0,K))
b1_memo = post_flip$b1_mean; se1_memo = post_flip$b1_sd; p1_memo = 2*exp(pnorm(abs(b1_memo/se1_memo), lower.tail=FALSE, log.p=TRUE))
b2_memo = post_flip$b2_mean; se2_memo = post_flip$b2_sd; p2_memo = 2*exp(pnorm(abs(b2_memo/se2_memo), lower.tail=FALSE, log.p=TRUE))
bci_1_memo = post_flip$bci1
bci_2_memo = post_flip$bci2

# MR2, multi-response
betaHat_Y = as.matrix(cbind(b_out_1, b_out_2))
betaHat_X = as.matrix(cbind(b_exp_1, b_exp_2))
colnames(betaHat_Y) = colnames(betaHat_X) = NULL
MR2_output = MR2(betaHat_Y, betaHat_X, EVgamma = c(1,2), niter = 7500, burnin = 2500, thin = 5, monitor = 1000)
# head(MR2_output$postMean$theta)
p_MR2_1 = MR2_output$samplerPar$pval[1,1]
p_MR2_2 = MR2_output$samplerPar$pval[2,2]
PostProc_output = PostProc(MR2_output, betaHat_Y, betaHat_X)
b_MR2_1 = PostProc_output$thetaPost[1,1]
b_MR2_2 = PostProc_output$thetaPost[2,2]
bci_MR2_1 = PostProc_output$thetaPost_CI[2,1,]
bci_MR2_2 = PostProc_output$thetaPost_CI[2,2,]


# summarize
res = c(K = K, b_fusiom = b2_memo, se_fusiom = se2_memo, p_fusiom = p2_memo, b_MR2_2, p_MR2_2)
res







