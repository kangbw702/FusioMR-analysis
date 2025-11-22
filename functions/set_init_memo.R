
# function to set up initial values for FusioMR_m method
# memo, uhp + chp

init_setup_memo = function(niter, K, alpha_1_init, alpha_2_init, beta_1_init, 
beta_2_init, eta_1_init, eta_2_init, pst_init) {
# starting values
theta_1_init = theta_2_init = rep(0, K) 
gamma_1_init = gamma_2_init = rep(0, K)
# eta_init = rbern(K, 0.5)
# MCMC trace
theta_1_tk = theta_2_tk = matrix(NA, nrow = niter, ncol = K)
gamma_1_tk = gamma_2_tk = matrix(NA, nrow = niter, ncol = K)
eta_1_tk = eta_2_tk = matrix(NA, nrow = niter, ncol = K)
pst_tk = matrix(NA, nrow = niter, ncol = 4)
alpha_1_tk = beta_1_tk = alpha_2_tk = beta_2_tk = rep(NA, niter)
theta_1_tk[1, ] = theta_1_init
theta_2_tk[1, ] = theta_2_init
gamma_1_tk[1, ] = gamma_1_init
gamma_2_tk[1, ] = gamma_2_init
eta_1_tk[1, ] = eta_1_init
eta_2_tk[1, ] = eta_2_init
alpha_1_tk[1] = alpha_1_init
alpha_2_tk[1] = alpha_2_init
beta_1_tk[1] = beta_1_init
beta_2_tk[1] = beta_2_init
pst_tk[1, ] = pst_init
# output
res = list(theta_1_tk = theta_1_tk, theta_2_tk = theta_2_tk, gamma_1_tk = gamma_1_tk, 
gamma_2_tk = gamma_2_tk, eta_1_tk = eta_1_tk, eta_2_tk = eta_2_tk, alpha_1_tk = alpha_1_tk, 
alpha_2_tk = alpha_2_tk, beta_1_tk = beta_1_tk, beta_2_tk = beta_2_tk, pst_tk = pst_tk)
return(res)
}


