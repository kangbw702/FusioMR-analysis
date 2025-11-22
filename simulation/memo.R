library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(MendelianRandomization)
library(Rcpp)
library(MR2)

source("dgf/dgm6.R")
source("functions/set_var_prior.R")
source("functions/set_init_seso.R")
source("functions/set_init_memo.R")
source("functions/utilities.R")
Rcpp::sourceCpp("fusiomr/gibbs_seso_uhpchp.cpp")
Rcpp::sourceCpp("fusiomr/gibbs_joint_rcpp_nopr_pp_bound.cpp")
Rcpp::sourceCpp("dgf/fastlm.cpp")

bern_prob <- function(rho, p, q) {
p00 = rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
p10 = 1 - q - p00
p01 = 1 - p - p00
p11 = p + q + p00 -  1
p00 = max(0, p00); p10 = max(0, p10); p01 = max(0, p01); p11 = max(0, p11)
return(c(p00, p01, p10, p11))
}

args <- commandArgs(trailingOnly = T)
param = read.csv(args[1], header = T)
offset = as.numeric(args[2])
args = as.numeric(unlist(param[offset,]))
print(args)

#param = read.csv("param_sim5.csv", header = T)
#args = as.numeric(unlist(param[9,]))

m = args[1]
nx1 = args[2]
nx2 = args[3]
ny1 = args[4]
ny2 = args[5]
a_gamma1 = args[6]*(-1)
b_gamma1 = args[6]
a_gamma2 = args[7]*(-1)
b_gamma2 = args[7]
rho_gamma = args[8]
a_f = args[9]
b_f = args[10]
q_uhp1 = args[11]
q_uhp2 = args[12]
a_alpha1 = args[13]*(-1)
b_alpha1 = args[13]
a_alpha2 = args[14]*(-1)
b_alpha2 = args[14]
rho_alpha = args[15]
a_phi1 = args[16]*(-1)
b_phi1 = args[16]
a_phi2 = args[17]*(-1)
b_phi2 = args[17]
rho_eta = args[18]
theta1 = args[19]
theta2 = args[20]
q_chp1 = args[21]
q_chp2 = args[22]
p_cutoff = args[23]
seed_k = args[24]
kappa_gamma_m = args[25]
kappa_theta_m = args[26]
kappa_gamma_s = args[27]
kappa_theta_s = args[28]
c_gamma = args[29]
c_theta = args[30]

t1 = Sys.time() 
out = NULL
nrep = 200
for (ii in 1:nrep) {
set.seed(ii+200*(seed_k-1))
if (ii%%1 == 0) print(ii)
# generate summary statistics
GWAS_summary = dgm6(m, nx1, nx2, ny1, ny2, a_gamma1, b_gamma1, a_gamma2, b_gamma2, rho_gamma, a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, a_phi1, b_phi1, a_phi2, b_phi2, rho_eta, theta1, theta2, q_uhp1, q_uhp2, q_chp1, q_chp2, beta_XU=1, beta_YU=1)
# select IV
b_exp_1_raw = c(GWAS_summary$b_exp_1)
b_exp_2_raw = c(GWAS_summary$b_exp_2)
b_out_1_raw = c(GWAS_summary$b_out_1)
b_out_2_raw = c(GWAS_summary$b_out_2)
se_exp_1_raw = c(GWAS_summary$se_exp_1)
se_exp_2_raw = c(GWAS_summary$se_exp_2)
se_out_1_raw = c(GWAS_summary$se_out_1)
se_out_2_raw = c(GWAS_summary$se_out_2)
z_gamma_1 = abs(b_exp_1_raw) / se_exp_1_raw
z_gamma_2 = abs(b_exp_2_raw) / se_exp_2_raw
p_gamma_1 = 2*(1 - pnorm(z_gamma_1))
p_gamma_2 = 2*(1 - pnorm(z_gamma_2))
iv_keep_1 = p_gamma_1 < p_cutoff
iv_keep_2 = p_gamma_2 < p_cutoff
K1 = sum(iv_keep_1)
K2 = sum(iv_keep_2)
iv_keep = iv_keep_1 | iv_keep_2
K = sum(iv_keep)

if (K>=5) {
# get summary statistics after selection
b_out_1s = b_out_1_raw[iv_keep_1]; se_out_1s = se_out_1_raw[iv_keep_1]
b_exp_1s = b_exp_1_raw[iv_keep_1]; se_exp_1s = se_exp_1_raw[iv_keep_1]
b_out_2s = b_out_2_raw[iv_keep_2]; se_out_2s = se_out_2_raw[iv_keep_2]
b_exp_2s = b_exp_2_raw[iv_keep_2]; se_exp_2s = se_exp_2_raw[iv_keep_2]
b_out_1 = b_out_1_raw[iv_keep]; se_out_1 = se_out_1_raw[iv_keep]
b_exp_1 = b_exp_1_raw[iv_keep]; se_exp_1 = se_exp_1_raw[iv_keep]
b_out_2 = b_out_2_raw[iv_keep]; se_out_2 = se_out_2_raw[iv_keep]
b_exp_2 = b_exp_2_raw[iv_keep]; se_exp_2 = se_exp_2_raw[iv_keep]

# run MR
# FusioMRm
niter = 20000
# hyper parameters
# parameters for variance priors
# use all SNPs before pval thresholding
vp = set_variance_priors_m2x2_diag(ghat_mat = matrix(c(b_exp_1_raw, b_exp_2_raw), ncol=2), gse_mat = matrix(c(se_exp_1_raw, se_exp_2_raw), ncol=2), Ghat_mat = matrix(c(b_out_1_raw, b_out_2_raw), ncol=2), Gse_mat = matrix(c(se_out_1_raw, se_out_2_raw), ncol=2), B0 = NULL, K = K, Kmin = 5, Kmax = 20, rho12 = 0, rho_gg = 0, rho_gj = list(c(0,0), c(0,0)), c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, kappa_gamma = kappa_gamma_m, kappa_theta = kappa_theta_m)

nu_gamma_prior = vp$gamma$nu; Phi_gamma_prior = vp$gamma$Phi
nu_theta_prior = vp$theta$nu; Phi_theta_prior = vp$theta$Phi
pst_init = bern_prob(rho_eta, q_chp1, q_chp2)
cc = rep(1, 4) # Dirichlet parameters

# starting values
# Use prior mean as init for Sigma_gamma and Sigma_theta
start_val = init_setup_memo(niter, K, alpha_1_init = 0.1, alpha_2_init = 0.1, beta_1_init = theta1, beta_2_init = theta2, eta_1_init = rep(0,K), eta_2_init = rep(0,K), pst_init = pst_init) 

# rcpp
res_memo = gibbs_joint_rcpp_nopr(niter, K, beta_1_tk = start_val$beta_1_tk, beta_2_tk = start_val$beta_2_tk, alpha_1_tk = start_val$alpha_1_tk, alpha_2_tk = start_val$alpha_2_tk, eta_1_tk = start_val$eta_1_tk, eta_2_tk = start_val$eta_2_tk, theta_1_tk = start_val$theta_1_tk, theta_2_tk = start_val$theta_2_tk, gamma_1_tk = start_val$gamma_1_tk, gamma_2_tk = start_val$gamma_2_tk, pst_tk = start_val$pst_tk, b_out_1, b_exp_1, (se_out_1)^2, (se_exp_1)^2, b_out_2, b_exp_2, (se_out_2)^2, (se_exp_2)^2, vp$gamma$prior_mean, vp$theta$prior_mean, nu_gamma_prior, Phi_gamma_prior, nu_theta_prior, Phi_theta_prior, cc) 

ids = (niter/2 + 1):niter
post_flip = label_flip_joint(niter, res_memo, rep(0,K), rep(0,K))
b1_memo = post_flip$b1_mean; se1_memo = post_flip$b1_sd
b2_memo = post_flip$b2_mean; se2_memo = post_flip$b2_sd
bci_1_memo = post_flip$bci1
bci_2_memo = post_flip$bci2
# sigma2_gamma1 = mean(res_memo$sigma2_gamma1_tk[ids], na.rm = T)
# sigma2_gamma2 = mean(res_memo$sigma2_gamma2_tk[ids], na.rm = T)
# sigma2_theta1 = mean(res_memo$sigma2_theta1_tk[ids], na.rm = T)
# sigma2_theta2 = mean(res_memo$sigma2_theta2_tk[ids], na.rm = T)


# FusioMRs empirical (outcome 2, smaller population)
niter = 20000
# parameters for variance priors
# use all SNPs before pval thresholding
vp = set_variance_priors(ghat = b_exp_2_raw, gse = se_exp_2_raw, Ghat = b_out_2_raw, Gse = se_out_2_raw, beta0 = NULL, K = K2, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = kappa_gamma_s, kappa_theta = kappa_theta_s)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
a_q_prior = b_q_prior = 1

# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K2, alpha_init = 1, beta_init = theta2, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean), q_init = 0)

# rcpp
res2 = gibbs_rcpp_nopr(niter, K2, beta_tk = start_val$beta_tk, alpha_tk = start_val$alpha_tk, eta_tk = start_val$eta_tk, theta_tk = start_val$theta_tk, gamma_tk = start_val$gamma_tk, q_tk = start_val$q_tk, b_out_2s, b_exp_2s, (se_out_2s)^2, (se_exp_2s)^2, sigma2_gamma_tk = start_val$sigma2_gamma_tk, sigma2_theta_tk = start_val$sigma2_theta_tk, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior, a_q_prior, b_q_prior) 

ids = (niter/2 + 1):niter
post_flip = label_flip(niter, res2, rep(0, K2))
b_s2 = post_flip$b_mean; se_s2 = post_flip$b_sd
bci_s2 = post_flip$bci


# output
tmp = c(seed = ii, K = K, K1 = K1, K2 = K2, nx1=nx1, nx2=nx2, ny1=ny1, ny2=ny2,
s2_gamma1_theo = (b_gamma1-a_gamma1)^2/12, s2_theta1_theo = (b_alpha1-a_alpha1)^2/12, 
s2_gamma1_emp = max(1e-3, (var(b_exp_1_raw) - mean(se_exp_1_raw)^2)), 
s2_theta1_emp = max(1e-6, (var(b_out_1_raw)-mean(se_out_1_raw^2))), 
b_fusiom = b2_memo, b_fusios = b_s2, se_fusiom = se2_memo, se_fusios = se_s2, cover_fusiom = bci_2_memo[2] > theta2 & bci_2_memo[1] < theta2, cover_fusios = bci_s2[2] > theta2 & bci_s2[1] < theta2, fusiom = bci_2_memo[2] < 0 | bci_2_memo[1] > 0, fusios = bci_s2[2] < 0 | bci_s2[1] > 0)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
colnames(out) = c('seed', 'K', 'K1', 'K2', 'nx1', 'nx2', 'ny1', 'ny2', 's2_gamma_theo', 's2_theta_theo', 's2_gamma_emp', 's2_theta_emp', 'b_fusiom', 'b_fusios', 'se_fusiom', 'se_fusios', 'cover_fusiom', 'cover_fusios', 'fusiom', 'fusios')

out_file_name = paste0("out/sim5d_m_", m, "_bgamma1_", b_gamma1, "_bgamma2_", b_gamma2, "_rhogamma_", rho_gamma, "_quhp1_", q_uhp1, "_quhp2_", q_uhp2, "_balpha1_", b_alpha1, "_balpha2_", b_alpha2, "_rhoalpha_", rho_alpha, "_bphi1_", b_phi1, "_bphi2_", b_phi2, "_rhoeta_", rho_eta, "_theta1_", theta1, "_theta2_", theta2, "_qchp1_", q_chp1, "_qchp2_", q_chp2, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)

cat("File saved!\n", out_file_name, "\n")

t2 = Sys.time()
cat(difftime(t2, t1, units = "mins"), "\n")





