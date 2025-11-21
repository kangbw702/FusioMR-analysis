library(data.table)
library(tidyverse)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(Rcpp)

source("dgf/dgm4.R")
source("functions/functions.R")
Rcpp::sourceCpp("code/gibbs_seso_uhp_only.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")
Rcpp::sourceCpp("code/fastlm.cpp", cacheDir ="/gpfs/data/linchen-lab/Bowei/")

args <- commandArgs(trailingOnly = T)
param = read.csv(args[1], header = T)
offset = as.numeric(args[2])
args = as.numeric(unlist(param[offset,]))
print(args)

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
c_gamma = args[14]
c_theta = args[15]


out = NULL
nrep = 1000
for (ii in 1:nrep) {
set.seed(ii+1000*(seed_k-1))
if (ii%%1 == 0) print(ii)
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
vp = set_variance_priors(ghat = b_exp_raw, gse = se_exp_raw, Ghat = b_out_raw, Gse = se_out_raw, beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.0, kappa_gamma = 1, kappa_theta = 1.05)
a_gamma_prior = vp$gamma$a
a_theta_prior = vp$theta$a
b_gamma_prior = vp$gamma$b
b_theta_prior = vp$theta$b
# starting values
# Use prior mean as init for sigma_gamma and sigma_theta
start_val = init_setup(niter, K, alpha_init = 1, beta_init = theta, sigma_gamma_init = sqrt(vp$gamma$prior_mean), sigma_theta_init = sqrt(vp$theta$prior_mean))
# rcpp
res2 = gibbs_seso_uhp_only_cpp(
niter, K, start_val$beta_tk, start_val$theta_tk, start_val$gamma_tk,
start_val$sigma2_gamma_tk, start_val$sigma2_theta_tk,
b_out, b_exp, (se_out)^2, (se_exp)^2,
a_gamma_prior, b_gamma_prior, a_theta_prior, b_theta_prior
)
ids = (niter/2 + 1):niter
bhat2 = mean(res2$beta_tk[ids], na.rm = T)
se_bhat2 = sd(res2$beta_tk[ids], na.rm = T)
# sigma2_gamma2 = mean(res2$sigma2_gamma_tk[ids])
# sigma2_theta2 = mean(res2$sigma2_theta_tk[ids])
# 1-coverage, bCI
bci = quantile(res2$beta_tk[ids], probs=c(0.025,0.975), na.rm=T)
# prop_neg = mean(res2$beta_tk[ids] < 0) # F_hat(0)
# pseudo_p = 2*min(prop_neg, 1-prop_neg)
# prop_small = mean(abs(res2$beta_tk[ids]) < se_bhat2*0.1)


# output
tmp = c(seed = ii, K = K, s2_gamma_theo = (b_gamma-a_gamma)^2/12, s2_theta_theo = (b_alpha-a_alpha)^2/12, s2_gamma_emp = max(1e-3, (var(b_exp_raw) - mean(se_exp_raw)^2)), s2_theta_emp = max(1e-6, (var(b_out_raw)-mean(se_out_raw^2))), b_fusio_emp = bhat2, se_fusio_emp = se_bhat2, cover_fusio_emp = bci[2] > theta & bci[1] < theta, fusio_emp = bci[2] < 0 | bci[1] > 0)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
colnames(out) = c('seed', 'K', 's2_gamma_theo', 's2_theta_theo', 's2_gamma_emp', 's2_theta_emp', 'b_fusio_emp', 'se_fusio_emp', 'cover_fusio_emp', 'fusio_emp')

out_file_name = paste0("out/sim2h_m_", m, "_nx_", nx, "_ny_", ny, "_bgamma_", b_gamma, "_quhp_", q_uhp, "_balpha_", b_alpha, "_theta_", theta, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_qchp_", q_chp, "_bphi_", b_phi, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)

cat("File saved!\n", out_file_name)



