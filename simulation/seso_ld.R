library(data.table)
library(tidyverse)
library(MASS)
library(mvtnorm)
library(LaplacesDemon)
library(dirmult)
library(invgamma)
library(MendelianRandomization)
library(Rcpp)
library(cisMRcML)

source("dgf/dgm4_ld.R")
source("functions/set_var_prior.R")
source("functions/set_init.R")
Rcpp::sourceCpp("fusiomr/fusiomr_s.cpp"")
Rcpp::sourceCpp("dgf/fastlm.cpp"")

args <- commandArgs(trailingOnly = T)

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
rho_ld = args[14]
c_gamma = args[15]
c_theta = args[16]

tt1 = Sys.time()
out = NULL
nrep = 20
for (ii in 1:nrep) {
set.seed(ii+20*(seed_k-1))
if (ii%%1 == 0) print(ii)
# generate summary statistics
GWAS_summary = dgm4_ld(m, nx, ny, a_gamma, b_gamma, a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp, beta_XU = 1, beta_YU = 1, rho_ld)
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
vp = set_variance_priors(ghat = b_exp_raw, gse = se_exp_raw, Ghat = b_out_raw, Gse = se_out_raw, beta0 = NULL, K = K, Kmin = 5, Kmax = 20, rho_ov = 0, c_gamma = c_gamma, c_theta = c_theta, global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL)
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

# ciscml
R = GWAS_summary$R_ld
R = R[iv_keep, iv_keep]
LD_inv = solve(R)
b_exp_cond = LD_inv %*% b_exp; b_out_cond = LD_inv %*% b_out
Sig_exp = LD_inv %*% (R * crossprod(t(se_exp))) %*% LD_inv;
Sig_out = LD_inv %*% (R * crossprod(t(se_out))) %*% LD_inv;
Sig_exp_inv = solve(Sig_exp); Sig_out_inv = solve(Sig_out)
ciscML_res = cismr_cML_DP(b_exp_cond, b_out_cond, Sig_exp_inv, Sig_out_inv, maxit=200, n = nx, random_start = 5, random_start_pert = 5)

ciscML0_b = ciscML_res$BIC_theta; ciscML0_se = ciscML_res$BIC_se
ciscMLDP0_b = ciscML_res$BIC_DP_theta; ciscMLDP0_se = ciscML_res$BIC_DP_se

# other competing methods
mrinput = MendelianRandomization::mr_input(bx = as.vector(b_exp), bxse = as.vector(se_exp), by = as.vector(b_out), byse = as.vector(se_out), correlation = R)
# IVW
ivw_res = MendelianRandomization::mr_ivw(mrinput, correl = TRUE)
b_ivw = ivw_res@Estimate; se_ivw = ivw_res@StdError
# MR-Egger
egger_res = try(MendelianRandomization::mr_egger(mrinput, correl=TRUE))
if(class(egger_res) != 'try-error' & !is.null(egger_res)) {
b_egger = egger_res@Estimate; se_egger = egger_res@StdError.Est
} else {b_egger = se_egger = NA}
# cML
cml = MendelianRandomization::mr_cML(mrinput, MA = TRUE, DP = FALSE, num_pert = 200, n = nx)
b_cml = cml$Estimate; se_cml = cml$StdError
# cML-DP
cml_dp = MendelianRandomization::mr_cML(mrinput, MA = TRUE, DP = TRUE, num_pert = 200, n = nx)
b_cml_dp = cml_dp$Estimate; se_cml_dp = cml_dp$StdError

# output
tmp = c(seed = ii, K = K, b_fusio_emp = bhat2, b_ivw = b_ivw, b_egger = b_egger, b_cml = b_cml, b_cml_dp = b_cml_dp, b_ciscml = ciscML0_b, b_ciscml_dp = ciscMLDP0_b, 
se_fusio_emp = se_bhat2, se_ivw = se_ivw, se_egger = se_egger, se_cml = se_cml, se_cml_dp = se_cml_dp, se_ciscml = ciscML0_se, se_ciscml_dp = ciscMLDP0_se,
cover_fusio_emp = bci[2] > theta & bci[1] < theta, cover_ivw = b_ivw + 1.96*se_ivw > theta & b_ivw - 1.96*se_ivw < theta, cover_egger = b_egger + 1.96*se_egger > theta & b_egger - 1.96*se_egger < theta, cover_cml = b_cml + 1.96* se_cml > theta & b_cml - 1.96* se_cml < theta, cover_cml_dp = b_cml_dp + 1.96* se_cml_dp > theta & b_cml_dp - 1.96* se_cml_dp < theta, cover_ciscml = ciscML0_b + 1.96* ciscML0_se > theta & ciscML0_b - 1.96* ciscML0_se < theta, cover_ciscml_dp = ciscMLDP0_b + 1.96* ciscMLDP0_se > theta & ciscMLDP0_b - 1.96* ciscMLDP0_se < theta,
fusio_emp = bci[2] < 0 | bci[1] > 0, ivw = b_ivw + 1.96*se_ivw < 0 | b_ivw - 1.96*se_ivw > 0, egger = b_egger + 1.96*se_egger < 0 | b_egger - 1.96*se_egger > 0, cml = b_cml + 1.96* se_cml < 0 | b_cml - 1.96* se_cml > 0, cml_dp = b_cml_dp + 1.96* se_cml_dp < 0 | b_cml_dp - 1.96* se_cml_dp > 0, ciscml = ciscML0_b + 1.96* ciscML0_se < 0 | ciscML0_b - 1.96* ciscML0_se > 0, ciscml_dp = ciscMLDP0_b + 1.96* ciscMLDP0_se < 0 | ciscMLDP0_b - 1.96* ciscMLDP0_se > 0)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
colnames(out) = c('seed', 'K', 'b_fusio_emp', 'b_ivw', 'b_egger', 'b_cml', 'b_cml_dp', 'b_ciscml', 'b_ciscml_dp', 'se_fusio_emp', 'se_ivw', 'se_egger', 'se_cml', 'se_cml_dp', 'se_ciscml', 'se_cml_cisdp', 'cover_fusio_emp', 'cover_ivw', 'cover_egger', 'cover_cml', 'cover_cml_dp', 'cover_ciscml', 'cover_ciscml_dp', 'fusio_emp', 'ivw', 'egger', 'cml', 'cml_dp', 'ciscml', 'ciscml_dp')

out_file_name = paste0("out/sim_m_", m, "_nx_", nx, "_ny_", ny, "_bgamma_", b_gamma, "_quhp_", q_uhp, "_balpha_", b_alpha, "_theta_", theta, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_qchp_", q_chp, "_bphi_", b_phi, "_rho_", rho_ld, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)
cat("File saved!\n", out_file_name, "\n")

tt2 = Sys.time()
cat(difftime(tt2, tt1, units="mins"), "mins\n")



