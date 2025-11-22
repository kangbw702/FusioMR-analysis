# Set variance priors for FusioMRs (both sigma^2_gamma and sigma^2_theta)
# LD-free; uses clumped, ~independent SNPs
#
# Inputs:
# ghat, gse   : QTL/eQTL effects and SEs (exposure)
# Ghat, Gse   : GWAS effects and SEs (outcome)
# beta0       : optional initial beta (if NULL, GLS init on outcome SEs)
# K           : number of IVs after filtering
# Kmin, Kmax  : floor and cap of K when calculating a
# rho_ov      : sampling correlation due to overlap in [-1,1] (default 0)
# c_gamma     : prior weight per IV for sigma^2_gamma (a_gamma = 1 + c_gamma*K/2)
# c_theta     : prior weight per IV for sigma^2_theta (a_theta = 1 + c_theta*K/2)
# global_mean_gamma, global_mean_theta : global EB centers (for hybrid mode)
# hybrid      : if TRUE, prior mean = eta*local + (1-eta)*global (eta = K/(K+kappa_hybrid))
# kappa_hybrid : pooling control (>0)
# z_thresh    : optional |Z_gamma| selection threshold used to pick QTLs (winner’s-curse fix).
#             Example: for p=5e-8 two-sided, use qnorm(1 - 5e-8/2).
# trim        : tail prob. for winsor
# kappa_gamma : tunning parameter for prior mean of sigma2_gamma 
# kappa_theta : tunning parameter for prior mean of sigma2_theta 
#
# Returns: list with sublists $gamma and $theta:
# $gamma/$theta: list(mom_local, a, b, prior_mean)
# plus shared K, eta, beta0, notes

set_variance_priors <- function(ghat, gse, Ghat, Gse, beta0 = NULL, K = NULL, Kmin = 5, Kmax = 20, 
rho_ov = 0, c_gamma = 0.5, c_theta = 0.8, global_mean_gamma = NULL, global_mean_theta = NULL, 
hybrid = FALSE, kappa_hybrid = 5, z_thresh = NULL, trim = 0.1, kappa_gamma = 1, kappa_theta = 1
){
stopifnot(length(ghat) == length(gse), length(Ghat) == length(Gse), length(ghat) == length(Ghat))
if (is.null(K)) K = length(ghat) # summary statistics are after filtering
notes = c()
# helper: truncated-normal inflation for Z^2 under two-sided |Z|>c selection
kappa_TN <- function(c) {
if (is.null(c) || c <= 0) return(1)
tail = 1 - pnorm(c)
1 + c * dnorm(c) / tail
}
infl_g = kappa_TN(z_thresh)

winsor <- function(x, p = trim) {
if (length(x) < 5) return(x)
qs = quantile(x, c(p, 1-p), na.rm = TRUE)
pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p = trim) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p = trim) median(winsor(x, p), na.rm = TRUE)

if (!is.null(z_thresh) && z_thresh > 0)
notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)", z_thresh, infl_g))

# sigma^2_gamma: local MoM (de-noised, TN-corrected if requested)
# gamma: winsorized mean(z^2) - 1, rescaled by robust se^2 scale
zg2 = (ghat / gse)^2
mom_gamma = wmed(gse^2, trim) * max(wmean(zg2, trim) - 1, 0)
#Zg = ghat / gse; mom_gamma = mean((gse * Zg)^2 - (gse^2) * infl_g); mom_gamma = max(mom_gamma, 0)

# beta0 (for theta); GLS init if not provided
if (is.null(beta0)) {
wG = 1 / (Gse^2)
xtx = sum((ghat^2) * wG); xty = sum((ghat * Ghat) * wG)
beta0 = if (xtx > 0) xty/xtx else 0
notes = c(notes, "beta0 estimated by GLS on outcome SEs")
} else {
notes = c(notes, "beta0 provided by user")
}
# sigma^2_theta: residual MoM (de-noised, overlap-aware, TN-correct gamma part)
r = Ghat - beta0 * ghat
sv = (Gse^2) + (beta0^2) * (gse^2) - 2 * beta0 * rho_ov * (Gse * gse)
if (!is.null(z_thresh) && z_thresh > 0) {
# conservative: apply TN correction to the gamma-driven sampling component
sv = sv + (beta0^2) * gse^2 * (infl_g - 1)
notes = c(notes, "TN correction applied to gamma component inside residual variance")
}
zr2 = (r^2) / sv
mom_theta = wmed(sv, trim) * max(wmean(zr2, trim) - 1, 0)

#ZG = Ghat / Gse; r = (Gse * ZG) - beta0 * (gse * Zg)
#samp_var_r = (Gse^2) + (beta0^2) * (gse^2) - 2 * beta0 * rho_ov * (Gse * gse)
#if (!is.null(z_thresh) && z_thresh > 0) {
## conservative: apply TN correction to the gamma-driven sampling component
#samp_var_r = samp_var_r + (beta0^2) * gse^2 * (infl_g - 1)
#notes = c(notes, "TN correction applied to gamma component inside residual variance")
#}
#mom_theta = mean( r^2 - samp_var_r )
#mom_theta = max(mom_theta, 0)
  
# shapes (pseudo-df) and hybrid prior means
a_gamma = 1 + c_gamma * min(max(K, Kmin), Kmax) / 2
a_theta = 1 + c_theta * min(max(K, Kmin), Kmax) / 2
if (a_gamma <= 1) a_gamma = 1 + 1e-3
if (a_theta <= 1) a_theta = 1 + 1e-3

if (isTRUE(hybrid)) {
if (is.null(global_mean_gamma) || is.null(global_mean_theta))
stop("hybrid=TRUE requires global_mean_gamma and global_mean_theta.")
eta = K / (K + kappa_hybrid)
m_gamma = eta * mom_gamma + (1 - eta) * global_mean_gamma
m_theta = eta * mom_theta + (1 - eta) * global_mean_theta
notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f", eta))
} else {
eta = NA
m_gamma = mom_gamma
m_theta = mom_theta
}

# conservative inflation 
m_gamma = m_gamma * kappa_gamma
m_theta = m_theta * kappa_theta

# scales b to make E[sigma^2] = prior mean (valid for a>1)
b_gamma = (a_gamma - 1) * m_gamma
b_theta = (a_theta - 1) * m_theta
  
# assemble outputs
list(
gamma = list(mom_local  = mom_gamma, a = a_gamma, b = b_gamma, prior_mean = b_gamma / (a_gamma - 1)),
theta = list(mom_local  = mom_theta, a = a_theta, b = b_theta, prior_mean = b_theta / (a_theta - 1)),
K = K, eta = eta, beta0 = beta0, notes  = unique(notes)
)
}



# Set variance priors for FusioMRm (shared exposure, 2 outcomes), LD-free
#
# Inputs:
# ghat, gse          : vectors (length K) of QTL effects and SEs (exposure)
# Ghat_mat, Gse_mat  : K x 2 matrices of GWAS effects and SEs for two outcomes (cols 1,2)
# beta0              : optional length-2 vector of initial betas (if NULL, GLS init per outcome)
# K, Kmin, Kmax      : number of IVs after filtering (auto from ghat if NULL); floor/cap for a/nu
# rho12              : sampling correlation between outcome 1 and 2 (due to sample overlap)
# rho1g, rho2g       : sampling correlations between outcome j and QTL (default 0)
# c_gamma, c_theta   : prior weight per IV for sigma^2_gamma and Sigma_theta
# global_mean_gamma  : optional scalar global EB center for sigma^2_gamma (for hybrid)
# global_mean_theta  : optional 2x2 global EB center for Sigma_theta (for hybrid)
# hybrid             : if TRUE, prior mean = eta*local + (1-eta)*global, eta = K/(K + kappa_hybrid)
# kappa_hybrid       : pooling control (>0)
# z_thresh           : optional |Z_gamma| threshold used to select QTLs (TN winner’s-curse fix)
# trim               : tail prob. for robust trimmed mean/median (set 0 for no trimming)
# kappa_gamma        : tunning parameter for prior mean of sigma2_gamma 
# kappa_theta        : tunning parameter for prior mean of sigma2_theta 
#
# Returns:
# list(
# gamma = list(mom_local, a, b, prior_mean),
# theta = list(mom_local (2x2), nu, Phi, prior_mean (2x2)),
# K, eta, beta0 (length-2), notes
# )
set_variance_priors_m2 <- function(ghat, gse, Ghat_mat, Gse_mat, beta0 = NULL, K = NULL, 
Kmin = 5, Kmax = 20, rho12 = 0, rho1g = 0, rho2g = 0, c_gamma = 0.5, c_theta = 0.8, 
global_mean_gamma = NULL, global_mean_theta = NULL, hybrid = FALSE, kappa_hybrid = 5, 
z_thresh = NULL, trim = 0.10, kappa_gamma = 1, kappa_theta = 1
) {
  
stopifnot(is.matrix(Ghat_mat), is.matrix(Gse_mat), ncol(Ghat_mat) == 2, ncol(Gse_mat) == 2)
stopifnot(length(ghat) == length(gse), nrow(Ghat_mat) == length(ghat), nrow(Gse_mat) == length(ghat))
if (is.null(K)) K = length(ghat)
notes = c()
  
# helpers 
# TN inflation E[Z^2 | |Z|>c]
kappa_TN <- function(c) {              
if (is.null(c) || c <= 0) return(1)
tail = 1 - pnorm(c)
1 + c * dnorm(c) / tail
}
# winsorize
winsor <- function(x, p) {
if (p <= 0 || length(x) < 5) return(x)
qs = stats::quantile(x, c(p, 1 - p), na.rm = TRUE)
pmin(pmax(x, qs[[1]]), qs[[2]])
}
wmean <- function(x, p) mean(winsor(x, p), na.rm = TRUE)
wmed  <- function(x, p) stats::median(winsor(x, p), na.rm = TRUE)
# PSD projection (eigen floor at 0)
proj_psd <- function(M) {              
M = as.matrix(M); M = (M + t(M)) / 2
ev = eigen(M, symmetric = TRUE)
vals = pmax(ev$values, 0)
ev$vectors %*% diag(vals, nrow = 2) %*% t(ev$vectors)
}
  
infl_g = kappa_TN(z_thresh)
if (!is.null(z_thresh) && z_thresh > 0)
notes = c(notes, sprintf("TN correction for gamma with |Z|>%.2f (infl=%.3f)", z_thresh, infl_g))

# sigma^2_gamma, local MoM (scalar, de-noised, TN-corrected if requested)
Zg = ghat / gse
# de-noised MoM; robust version (trim) or exact subtraction with TN factor
if (trim > 0) {
mom_gamma = wmed(gse^2, trim) * max(wmean(Zg^2, trim) - 1, 0)
} else {
mom_gamma = mean((gse * Zg)^2 - gse^2 * infl_g)
mom_gamma = max(mom_gamma, 0)
}
  
# beta0 for each outcome (GLS init if missing)
if (is.null(beta0)) {
beta0 = numeric(2)
for (j in 1:2) {
wj = 1 / (Gse_mat[, j]^2)
xtx = sum(ghat^2 * wj)
xty = sum(ghat * Ghat_mat[, j] * wj)
beta0[j] = if (xtx > 0) xty / xtx else 0
}
notes = c(notes, "beta0 estimated by GLS per outcome")
} else {
stopifnot(length(beta0) == 2)
notes = c(notes, "beta0 provided by user")
}
  
# Sigma_theta (2x2) local MoM 
# residuals r_k (2x1), sampling covariance S_k (2x2) with overlap and TN fix on gamma parts
sg = gse
s1 = Gse_mat[, 1]; s2 = Gse_mat[, 2]
b1 = beta0[1]; b2 = beta0[2]
Z1 = Ghat_mat[, 1] / s1; Z2 = Ghat_mat[, 2] / s2
r1 = (s1*Z1) - b1 * (sg*Zg)
r2 = (s2*Z2) - b2 * (sg*Zg)
  
S11 = (s1^2) + (b1^2) * (sg^2) - 2 * b1 * rho1g * (s1*sg)
S22 = (s2^2) + (b2^2) * (sg^2) - 2 * b2 * rho2g * (s2*sg)
S12 = rho12 * (s1*s2) + b1*b2 * (sg^2) - b1 * rho2g * (s2*sg) - b2 * rho1g * (s1*sg)
if (!is.null(z_thresh) && z_thresh > 0) {
S11 = S11 + (b1^2) * sg^2 * (infl_g-1)
S22 = S22 + (b2^2) * sg^2 * (infl_g-1)
S12 = S12 + (b1*b2)* sg^2 * (infl_g-1)
notes = c(notes, "TN correction applied to gamma components inside residual sampling covariance")
}
  
# local MoM covariance: average of (r r^T - S_k), PSD-projected
M11 = mean(r1*r1 - S11)
M22 = mean(r2*r2 - S22)
M12 = mean(r1*r2 - S12)
mom_theta_mat = matrix(c(M11, M12, M12, M22), 2, 2)
mom_theta_mat = proj_psd(mom_theta_mat)

# pseudo-df and hybrid centers 
Keff = min(max(K, Kmin), Kmax)
a_gamma = 1 + (c_gamma * Keff) / 2
# For IW(p=2), need nu > p+1 = 3 for mean to exist
nu_theta = 3 + (c_theta * Keff)
if (a_gamma <= 1) a_gamma = 1 + 1e-3
if (nu_theta <= 3) nu_theta = 3 + 1e-3
  
if (isTRUE(hybrid)) {
if (is.null(global_mean_gamma) || is.null(global_mean_theta))
stop("hybrid=TRUE requires global_mean_gamma (scalar) and global_mean_theta (2x2).")
if (!is.matrix(global_mean_theta) || any(dim(global_mean_theta) != c(2, 2)))
stop("global_mean_theta must be a 2x2 matrix.")
eta = K / (K + kappa_hybrid)
m_gamma = eta * mom_gamma + (1 - eta) * global_mean_gamma
m_theta = eta * mom_theta_mat + (1 - eta) * global_mean_theta
notes = c(notes, sprintf("hybrid prior means with eta=K/(K+kappa)=%.3f", eta))
} else {
eta = NA_real_
m_gamma = mom_gamma
m_theta = mom_theta_mat
}

# conservative inflation 
m_gamma = m_gamma * kappa_gamma
m_theta = m_theta * kappa_theta 

# ensure m_theta is PSD
m_theta = proj_psd(m_theta)

# scales so that E[sigma^2_gamma] = m_gamma and E[Sigma_theta] = m_theta
b_gamma = (a_gamma - 1) * m_gamma
Phi = (nu_theta - 3) * m_theta
  
# outputs 
list(
gamma = list(mom_local = mom_gamma, a = a_gamma, b = b_gamma, prior_mean = b_gamma / (a_gamma - 1)),
theta = list(mom_local = mom_theta_mat, nu = nu_theta, Phi = Phi, prior_mean = Phi / (nu_theta - 3)),
K = K, eta = eta, beta0 = beta0, notes = unique(notes)
)
}




