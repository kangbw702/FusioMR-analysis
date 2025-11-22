library(data.table)
library(tidyverse)
library(Rcpp)
library(cisMRcML)

source("dgf/dgm4.R")
Rcpp::sourceCpp("dgf/fastlm.cpp")

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
c_gamma = args[14]
c_theta = args[15]


out = NULL
nrep = 200
for (ii in 1:nrep) {
set.seed(ii+200*(seed_k-1))
if (ii%%1 == 0) print(ii)
# generate summary statistics
GWAS_summary = dgm4(m, nx, ny, a_gamma, b_gamma, a_f, b_f, a_alpha, b_alpha, a_phi, b_phi, theta, q_uhp, q_chp)
# get summary statistics 
#b_exp = c(GWAS_summary$b_exp); se_exp = c(GWAS_summary$se_exp)
#b_out = c(GWAS_summary$b_out); se_out = c(GWAS_summary$se_out)
# select IV
b_exp_raw = c(GWAS_summary$b_exp); se_exp_raw = c(GWAS_summary$se_exp)
b_out_raw = c(GWAS_summary$b_out); se_out_raw = c(GWAS_summary$se_out)
z_gamma = abs(b_exp_raw) / se_exp_raw
z_Gamma = abs(b_out_raw) / se_out_raw
p_gamma = 2*(1 - pnorm(z_gamma))
p_Gamma = 2* exp(pnorm(z_Gamma, lower.tail =FALSE, log.p = TRUE))
#sum(p_Gamma<5e-8)
ivx = p_gamma < p_cutoff; 
#ivy = p_Gamma < 1e-20; 
ivy = (p_Gamma < 5e-8) & (rank(p_Gamma, ties.method = "first") <= 2)
iv_keep = (ivx|ivy)
K = sum(iv_keep); Kexp = sum(ivx)
print(K)
if (K>=5) {
# get summary statistics after selection
b_out = b_out_raw[iv_keep]; se_out = se_out_raw[iv_keep]
b_exp = b_exp_raw[iv_keep]; se_exp = se_exp_raw[iv_keep]

# run MR
# ciscml-dp
R = diag(K)
LD_inv = solve(R)
b_exp_cond = LD_inv %*% b_exp; b_out_cond = LD_inv %*% b_out
Sig_exp = LD_inv %*% (R * crossprod(t(se_exp))) %*% LD_inv;
Sig_out = LD_inv %*% (R * crossprod(t(se_out))) %*% LD_inv;
Sig_exp_inv = solve(Sig_exp); Sig_out_inv = solve(Sig_out)
ciscML_res = cismr_cML_DP(b_exp_cond, b_out_cond, Sig_exp_inv, Sig_out_inv, maxit=200, n = nx, random_start = 5, random_start_pert = 5)

ciscML0_b = ciscML_res$BIC_theta; ciscML0_se = ciscML_res$BIC_se
ciscMLDP0_b = ciscML_res$BIC_DP_theta; ciscMLDP0_se = ciscML_res$BIC_DP_se


# output
tmp = c(seed = ii, K = K, Kexp = Kexp, s2_gamma_theo = (b_gamma-a_gamma)^2/12, s2_theta_theo = (b_alpha-a_alpha)^2/12, s2_gamma_emp = max(1e-3, (var(b_exp_raw) - mean(se_exp_raw)^2)), s2_theta_emp = max(1e-6, (var(b_out_raw)-mean(se_out_raw^2))), b_ciscml = ciscML0_b, b_ciscml_dp = ciscMLDP0_b, se_ciscml = ciscML0_se, se_ciscml_dp = ciscMLDP0_se, cover_ciscml = ciscML0_b + 1.96* ciscML0_se > theta & ciscML0_b - 1.96* ciscML0_se < theta, cover_ciscml_dp = ciscMLDP0_b + 1.96* ciscMLDP0_se > theta & ciscMLDP0_b - 1.96* ciscMLDP0_se < theta, ciscml = ciscML0_b + 1.96* ciscML0_se < 0 | ciscML0_b - 1.96* ciscML0_se > 0, ciscml_dp = ciscMLDP0_b + 1.96* ciscMLDP0_se < 0 | ciscMLDP0_b - 1.96* ciscMLDP0_se > 0)
out = rbind(out, tmp)
} # end if
} # end loop

out = as.data.frame(out)
if (length(out)) {
colnames(out) = c('seed', 'K', 'Kexp', 's2_gamma_theo', 's2_theta_theo', 's2_gamma_emp', 's2_theta_emp', 'b_ciscml', 'b_ciscml_dp', 'se_ciscml', 'se_cml_cisdp', 'cover_ciscml', 'cover_ciscml_dp', 'ciscml', 'ciscml_dp')

out_file_name = paste0("out/simcis2g_m_", m, "_nx_", nx, "_ny_", ny, "_bgamma_", b_gamma, "_quhp_", q_uhp, "_balpha_", b_alpha, "_theta_", theta, "_pcut_", p_cutoff, "_cgamma_", c_gamma, "_ctheta_", c_theta, "_qchp_", q_chp, "_bphi_", b_phi, "_seed_", seed_k, ".txt")
fwrite(out, out_file_name)
} else {
cat("Empty!\n")
}

cat("File saved!\n", out_file_name)



