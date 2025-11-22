library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(MendelianRandomization)
library(Rcpp)
library(MR2)

source("dgf/dgm5.R")
source("functions/set_var_prior.R")
source("functions/set_init_seso.R")
source("functions/set_init_semo.R")
Rcpp::sourceCpp("fusiomr/fusiomr_s.cpp")
Rcpp::sourceCpp("fusiomr/fusiomr_m_shared_exp.cpp")
Rcpp::sourceCpp("dgf/fastlm.cpp")

args <- commandArgs(trailingOnly = T)

m = args[1]
nx = args[2]
ny1 = args[3]
ny2 = args[4]
a_gamma = args[5]*(-1)
b_gamma = args[5]
a_f = args[6]
b_f = args[7]
q_uhp1 = args[8]
q_uhp2 = args[9]
a_alpha1 = args[10]*(-1)
b_alpha1 = args[10]
a_alpha2 = args[11]*(-1)
b_alpha2 = args[11]
rho_alpha = args[12]
theta1 = args[13]
theta2 = args[14]
p_cutoff = args[15]
seed_k = args[16]
kappa_gamma_m = args[17]
kappa_theta_m = args[18]
kappa_gamma_s = args[19]
kappa_theta_s = args[20]
c_gamma = args[21]
c_theta = args[22]

out = NULL
nrep = 1000
for (ii in 1:nrep) {
set.seed(ii+1000*(seed_k-1))
if (ii%%100 == 0) print(ii)
# generate summary statistics
GWAS_summary = dgm5(m, nx, ny1, ny2, a_gamma, b_gamma, a_f, b_f, a_alpha1, b_alpha1, a_alpha2, b_alpha2, rho_alpha, theta1, theta2, q_uhp1, q_uhp2)
# select IV
b_exp_raw = c(GWAS_summary$b_exp); se_exp_raw = c(GWAS_summary$se_exp)
b_out_1_raw = c(GWAS_summary$b_out_1); se_out_1_raw = c(GWAS_summary$se_out_1)
b_out_2_raw = c(GWAS_summary$b_out_2); se_out_2_raw = c(GWAS_summary$se_out_2)
z_gamma = abs(b_exp_raw) / se_exp_raw
p_gamma = 2*(1 - pnorm(z_gamma))
iv_keep = p_gamma < p_cutoff
K = sum(iv_keep)

if (K>=5) {
# get summary statistics
b_out_1 = b_out_1_raw[iv_keep]; se_out_1 = se_out_1_raw[iv_keep]
b_out_2 = b_out_2_raw[iv_keep]; se_out_2 = se_out_2_raw[iv_keep]
b_exp = b_exp_raw[iv_keep]; se_exp = se_exp_raw[iv_keep]
    
# run MR
# FusioMRm semo
niter = 20000
# hyper parameters
# parameters for variance priors
# use all SNPs before pval thresholding
vp = set_variance_priors_m2(ghat = b_exp_raw, gse = se_exp_raw, Ghat_mat = matrix(c(b_out_1_raw, b_out_2_raw), ncol=2), Gse_mat = matrix(c(se_out_1_raw, se_out_2_raw), ncol=2), beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = kappa_gamma_m, kappa_theta = kappa_theta_m)
a_gamma_prior = vp$gamma$a
b_gamma_prior = vp$gamma$b
nu_theta_prior = vp$theta$nu
Phi_theta_prior = vp$theta$Phi
# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup_semo_uhp_only(niter, K, beta_1_init = theta1, beta_2_init = theta2, sigma_gamma_init = sqrt(vp$gamma$prior_mean))
# rcpp
res_semo = gibbs_semo_uhp_only_rcpp(niter, K, start_val$beta_1_tk, start_val$beta_2_tk, start_val$theta_1_tk, start_val$theta_2_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, b_out_1, b_out_2, (se_out_1)^2, (se_out_2)^2, b_exp, (se_exp)^2, a_gamma_prior, b_gamma_prior, vp$theta$prior_mean, nu_theta_prior, Phi_theta_prior)
ids = (niter/2 + 1):niter
b1_semo = mean(res_semo$beta_1_tk[ids], na.rm = T)
se1_semo = sd(res_semo$beta_1_tk[ids], na.rm = T)
b2_semo = mean(res_semo$beta_2_tk[ids], na.rm = T)
se2_semo = sd(res_semo$beta_2_tk[ids], na.rm = T)
bci_1_semo = quantile(res_semo$beta_1_tk[ids], probs=c(0.025,0.975), na.rm=T)
bci_2_semo = quantile(res_semo$beta_2_tk[ids], probs=c(0.025,0.975), na.rm=T)
# sigma2_gamma = mean(res_semo$sigma2_gamma_tk[ids], na.rm = T)
# sigma2_theta1 = mean(res_semo$sigma2_theta1_tk[ids], na.rm = T)
# sigma2_theta2 = mean(res_semo$sigma2_theta2_tk[ids], na.rm = T)

# FusioMRs empirical (outcome 2)
niter = 20000
# hyper parameters
# parameters for variance priors
# use all SNPs before pval thresholding
vp = set_variance_priors(ghat = b_exp_raw, gse = se_exp_raw, Ghat = b_out_2_raw, Gse = se_out_2_raw, beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = kappa_gamma_s, kappa_theta = kappa_theta_s)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
# starting values. Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K, alpha_init = 1, beta_init = theta2, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean))
# rcpp
res1 = gibbs_seso_uhp_only_cpp(niter, K, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk, start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk, b_out_2, b_exp, (se_out_2)^2, (se_exp)^2, a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior)
ids = (niter/2 + 1):niter
bhat_s1 = mean(res1$beta_tk[ids], na.rm = T)
se_bhat_s1 = sd(res1$beta_tk[ids], na.rm = T)
# sigma2_gamma2 = mean(res2$sigma2_gamma_tk[ids])
# sigma2_theta2 = mean(res2$sigma2_theta_tk[ids])
# 1-coverage, bCI
bci_s1 = quantile(res1$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)
# prop_neg = mean(res2$beta_tk[ids] < 0) # F_hat(0)
# pseudo_p = 2*min(prop_neg, 1-prop_neg)
# prop_small = mean(abs(res2$beta_tk[ids]) < se_bhat2*0.1)

# competing methods (single-outcome MR, outcome 2)
mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out_2, byse = se_out_2)
# IVW fixed
IVW_f = MendelianRandomization::mr_ivw(mr.obj, model = 'fixed')
b_ivw = IVW_f$Estimate; se_ivw = IVW_f$StdError
# MR-Egger
Egger = try(MendelianRandomization::mr_egger(mr.obj))
if(class(Egger) != 'try-error' & !is.null(Egger)) {
b_egger = Egger$Estimate; se_egger = Egger$StdError.Est
} else {b_egger = se_egger = NA}
# cML
cml = MendelianRandomization::mr_cML(mr.obj, MA = T, DP = F, num_pert = 200, n = nx)
b_cml = cml$Estimate; se_cml = cml$StdError
# cML-DP
cml_dp = MendelianRandomization::mr_cML(mr.obj, MA = T, DP = T, num_pert = 200, n = nx)
b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError

# competing methods (multi-response)
# MR2
betaHat_Y = as.matrix(cbind(b_out_1, b_out_2))
betaHat_X = as.matrix(cbind(b_exp))
colnames(betaHat_Y) = colnames(betaHat_X) = NULL
MR2_output = MR2(betaHat_Y, betaHat_X, EVgamma = 0.5, niter = 7500, burnin = 2500, thin = 5, monitor = 1000)
# head(MR2_output$postMean$theta)
#p_MR2_1 = MR2_output$samplerPar$pval[1]
p_MR2_2 = MR2_output$samplerPar$pval[2]
PostProc_output = PostProc(MR2_output, betaHat_Y, betaHat_X)
#b_MR2_1 = PostProc_output$thetaPost[1]
b_MR2_2 = PostProc_output$thetaPost[2]
#bci_MR2_1 = PostProc_output$thetaPost_CI[1,1,]
bci_MR2_2 = PostProc_output$thetaPost_CI[1,2,]

# output
tmp = c(seed = ii, K = K, 
b_fusiom = b2_semo, b_fusios = bhat_s1, b_ivw = b_ivw, b_egger = b_egger, b_cml = b_cml, b_cml_dp = b_cml_dp, b_MR2 = b_MR2_2, 
se_fusiom = se2_semo, se_fusios = se_bhat_s1, se_ivw = se_ivw, se_egger = se_egger, se_cml = se_cml, se_cml_dp = se_cml_dp, se_MR2 = NA, 
cover_fusiom = bci_2_semo[2] > theta2 & bci_2_semo[1] < theta2, cover_fusios = bci_s1[2] > theta2 & bci_s1[1] < theta2, cover_ivw = b_ivw + 1.96*se_ivw > theta2 & b_ivw - 1.96*se_ivw < theta2, cover_egger = b_egger + 1.96*se_egger > theta2 & b_egger - 1.96*se_egger < theta2, cover_cml = b_cml + 1.96*se_cml > theta2 & b_cml - 1.96*se_cml < theta2, cover_cml_dp = b_cml_dp + 1.96*se_cml_dp > theta2 & b_cml_dp - 1.96*se_cml_dp < theta2, cover_MR2 = bci_MR2_2[2] > theta2 & bci_MR2_2[1] < theta2,
fusiom = bci_2_semo[2] < 0 | bci_2_semo[1] > 0, fusios = bci_s1[2] < 0 | bci_s1[1] > 0, ivw = b_ivw + 1.96*se_ivw < 0 | b_ivw - 1.96*se_ivw > 0, egger = b_egger + 1.96*se_egger < 0 | b_egger - 1.96*se_egger > 0, cml = b_cml + 1.96*se_cml < 0 | b_cml - 1.96*se_cml > 0, cml_dp = b_cml_dp + 1.96*se_cml_dp < 0 | b_cml_dp - 1.96*se_cml_dp > 0, MR2 = p_MR2_2 < 0.05)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
colnames(out) = c('seed', 'K',    )

colnames(out) = c('seed', 'K', 'b_fusiom', 'b_fusios', 'b_ivw', 'b_egger', 'b_cml', 'b_cml_dp', 'b_MR2', 
'se_fusiom', 'se_fusios', 'se_ivw', 'se_egger', 'se_cml', 'se_cml_dp', 'se_MR2', 
'cover_fusiom', 'cover_fusios', 'cover_ivw', 'cover_egger', 'cover_cml', 'cover_cml_dp', 'cover_MR2', 
'fusiom', 'fusios', 'ivw', 'egger', 'cml', 'cml_dp', 'MR2'
)

out_file_name = paste0("out/sim_m_", m, "_nx_", nx, "_ny1_", ny1, "_ny2_", ny2, "_bgamma_", b_gamma, "_quhp1_", q_uhp1, "_quhp2_", q_uhp2, "_balpha1_", b_alpha1, "_balpha2_", b_alpha2, "_rhoalpha_", rho_alpha, "_theta1_", theta1, "_theta2_", theta2, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)

cat("File saved!\n", out_file_name)




