library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(MendelianRandomization)
library(Rcpp)

source("dgf/dgm4.R")
source("functions/set_var_prior.R")
source("functions/set_init_seso.R")
source("code/label_flip.R")
Rcpp::sourceCpp("code/gibbs_seso_uhpchp.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")
Rcpp::sourceCpp("code/fastlm.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")

args <- commandArgs(trailingOnly = T)
param = read.csv(args[1], header = T)
offset = as.numeric(args[2])
args = as.numeric(unlist(param[offset,]))
print(args)

#param = read.csv("param_sim3.csv", header = T)
#offset = 1
#args = as.numeric(unlist(param[offset,]))

m = args[1]
nx = args[2]
ny = args[3]
a_gamma = args[4]*(-1)
b_gamma = args[4]
a_f = args[5]
b_f = args[6]
q_uhp = args[7]
a_alpha = args[8]*(-1)
b_alpha = args[8]
theta = args[9]
p_cutoff = args[10]
seed_k = args[11]
q_chp = args[12]
a_phi = args[13]*(-1)
b_phi = args[13]
kappa_gamma = args[14]
kappa_theta = args[15]
c_gamma = args[16]
c_theta = args[17]

out = NULL
nrep = 200
for (ii in 1:nrep) {
set.seed(ii+200*(seed_k-1))
if (ii%%10 == 0) print(ii)
# generate summary statistics
GWAS_summary = dgm4(m, nx, ny, a_gamma, b_gamma, a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp)
# get summary statistics 
b_exp = c(GWAS_summary$b_exp); se_exp = c(GWAS_summary$se_exp)
b_out = c(GWAS_summary$b_out); se_out = c(GWAS_summary$se_out)
# select IV
b_exp_raw = c(GWAS_summary$b_exp); se_exp_raw = c(GWAS_summary$se_exp)
b_out_raw = c(GWAS_summary$b_out); se_out_raw = c(GWAS_summary$se_out)
z_gamma = abs(b_exp_raw) / se_exp_raw
p_gamma = 2*(1 - pnorm(z_gamma))
iv_keep = p_gamma < p_cutoff
K = sum(iv_keep)
if (K>=5) {
# get summary statistics after selection
b_out = b_out_raw[iv_keep]; se_out = se_out_raw[iv_keep]
b_exp = b_exp_raw[iv_keep]; se_exp = se_exp_raw[iv_keep]

# run MR
# FusioMR empirical
niter = 20000
# parameters for variance priors
# use all SNPs before pval thresholding
vp = set_variance_priors(ghat = b_exp_raw, gse = se_exp_raw, Ghat = b_out_raw, Gse = se_out_raw, beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = kappa_gamma, kappa_theta = kappa_theta)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
a_q_prior = b_q_prior = 1
# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K, alpha_init = 1, beta_init = theta, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean))
# rcpp
res2 = gibbs_rcpp_nopr(niter, K, beta_tk = start_val$beta_tk, alpha_tk = start_val$alpha_tk,
eta_tk = start_val$eta_tk, theta_tk = start_val$theta_tk, gamma_tk = start_val$gamma_tk, 
q_tk = start_val$q_tk, b_out, b_exp, (se_out)^2, (se_exp)^2, 
sigma2_gamma_tk = start_val$sigma2_gamma_tk, sigma2_theta_tk = start_val$sigma2_theta_tk,
a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior, a_q_prior, b_q_prior) 

ids = (niter/2 + 1):niter
post_flip = label_flip(niter, res2, rep(0, K))
bhat2 = post_flip$b_mean
se_bhat2 = post_flip$b_sd
# 1-coverage, bCI
# bci = quantile(res2$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)
bci = post_flip$bci
# prop_neg = mean(res2$beta_tk[ids] < 0) # F_hat(0)
# pseudo_p = 2*min(prop_neg, 1-prop_neg)
# prop_small = mean(abs(res2$beta_tk[ids]) < se_bhat2*0.1)

# competing methods
mr.obj = MendelianRandomization::mr_input(bx = b_exp, bxse = se_exp, by = b_out, byse = se_out)
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
#cml_dp = MendelianRandomization::mr_cML(mr.obj, MA = T, DP = T, num_pert = 200, n = nx)
#b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError
b_cml_dp = 0; se_cml_dp = 1

# output
tmp = c(seed = ii, K = K, s2_gamma_theo = (b_gamma-a_gamma)^2/12, s2_theta_theo = (b_alpha-a_alpha)^2/12, s2_gamma_emp = max(1e-3, (var(b_exp_raw) - mean(se_exp_raw)^2)), s2_theta_emp = max(1e-6, (var(b_out_raw)-mean(se_out_raw^2))), 
b_fusio = bhat2, b_ivw = b_ivw, b_egger = b_egger, b_cml = b_cml, b_cml_dp = b_cml_dp, 
se_fusio = se_bhat2, se_ivw = se_ivw, se_egger = se_egger, se_cml = se_cml, se_cml_dp = se_cml_dp,
cover_fusio = bci[2] > theta & bci[1] < theta, cover_ivw = b_ivw + 1.96*se_ivw > theta & b_ivw - 1.96*se_ivw < theta, cover_egger = b_egger + 1.96*se_egger > theta & b_egger - 1.96*se_egger < theta, cover_cml = b_cml + 1.96*se_cml > theta & b_cml - 1.96*se_cml < theta, cover_cml_dp = b_cml_dp + 1.96*se_cml_dp > theta & b_cml_dp - 1.96*se_cml_dp < theta,
fusio = bci[2] < 0 | bci[1] > 0, ivw = b_ivw + 1.96*se_ivw < 0 | b_ivw - 1.96*se_ivw > 0, egger = b_egger + 1.96*se_egger < 0 | b_egger - 1.96*se_egger > 0, cml = b_cml + 1.96*se_cml < 0 | b_cml - 1.96*se_cml > 0, cml_dp = b_cml_dp + 1.96*se_cml_dp < 0 | b_cml_dp - 1.96*se_cml_dp > 0)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
colnames(out) = c('seed', 'K', 's2_gamma_theo', 's2_theta_theo', 's2_gamma_emp', 's2_theta_emp', 'b_fusio', 'b_ivw', 'b_egger', 'b_cml', 'b_cml_dp', 'se_fusio', 'se_ivw', 'se_egger', 'se_cml', 'se_cml_dp', 'cover_fusio', 'cover_ivw', 'cover_egger', 'cover_cml', 'cover_cml_dp', 'fusio', 'ivw', 'egger', 'cml', 'cml_dp')

out_file_name = paste0("out/sim7c_m_", m, "_nx_", nx, "_ny_", ny, "_bgamma_", b_gamma, "_quhp_", q_uhp, "_balpha_", b_alpha, "_theta_", theta, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_qchp_", q_chp, "_bphi_", b_phi, "_kgamma_", kappa_gamma, "_ktheta_", kappa_theta, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)

cat("File saved!\n", out_file_name)



