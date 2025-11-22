#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// single population MR


// rinvgamma from R::invgamma

// [[Rcpp::export]]
NumericVector my_rinvgamma(int n, double shape, double rate) {
  Function ff("rinvgamma"); 
  NumericVector res = ff(n, Named("shape") = shape, _["rate"] = rate);
  return res;
}

// [[Rcpp::export]]
NumericVector my_rinvwishart(double nu, mat S) {
  Function ff("rinvwishart"); 
  NumericVector res = ff(Named("nu") = nu, _["S"] = S);
  return res;
}

// [[Rcpp::export]]
NumericVector my_rdirichlet(int n, NumericVector alpha) {
  Function ff("rdirichlet"); 
  NumericVector res = ff(Named("n") = n, _["alpha"] = alpha);
  return res;
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// Gibbs sampling

// [[Rcpp::export]]
  
List gibbs_joint_rcpp_nopr(int niter, int K,
                           vec beta_1_tk, vec beta_2_tk, vec alpha_1_tk, vec alpha_2_tk,
                           mat eta_1_tk, mat eta_2_tk, mat theta_1_tk, mat theta_2_tk, 
                           mat gamma_1_tk, mat gamma_2_tk, mat pst_tk,
                           vec Gamma_hat_1, vec gamma_hat_1, vec s2_hat_Gamma_1, vec s2_hat_gamma_1,
                           vec Gamma_hat_2, vec gamma_hat_2, vec s2_hat_Gamma_2, vec s2_hat_gamma_2,
                           mat Sigma_gamma_init, mat Sigma_theta_init,
                           double m_gamma, mat V_gamma, double m_theta, mat V_theta,
                           vec cc
  ) {
  
  // set current values as initial values
  mat theta_cur(K,2); 
  for (int ii = 0; ii < K; ii ++) { theta_cur(ii,0) = theta_1_tk(0,ii); theta_cur(ii,1) = theta_2_tk(0,ii); }
  mat gamma_cur(K,2); 
  for (int ii = 0; ii < K; ii ++) { gamma_cur(ii,0) = gamma_1_tk(0,ii); gamma_cur(ii,1) = gamma_2_tk(0,ii); }
  vec alpha_cur(2); alpha_cur[0] = alpha_1_tk[0]; alpha_cur[1] = alpha_2_tk[0];
  vec beta_cur(2); beta_cur[0] = beta_1_tk[0]; beta_cur[1] = beta_2_tk[0];
  mat eta_cur(K,2);
  for (int ii = 0; ii < K; ii ++) { eta_cur(ii,0) = eta_1_tk(0,ii); eta_cur(ii,1) = eta_2_tk(0,ii); }
  mat Sigma_gamma_cur = Sigma_gamma_init;
  mat Sigma_theta_cur = Sigma_theta_init;
  vec pst_cur = pst_tk.row(0).t();
  NumericMatrix q_post_1_tk(niter, K);
  NumericVector q_post_1_cur = q_post_1_tk(0, _);
  NumericMatrix q_post_2_tk(niter, K);
  NumericVector q_post_2_cur = q_post_2_tk(0, _);
  std::vector<double> sigma2_gamma1_tk(niter, 0.0);
  std::vector<double> sigma2_gamma2_tk(niter, 0.0);
  sigma2_gamma1_tk[0] = Sigma_gamma_init(0,0);
  sigma2_gamma2_tk[0] = Sigma_gamma_init(0,0);
  std::vector<double> sigma2_theta1_tk(niter, 0.0);
  std::vector<double> sigma2_theta2_tk(niter, 0.0);
  sigma2_theta1_tk[0] = Sigma_theta_init(0,0);
  sigma2_theta2_tk[0] = Sigma_theta_init(0,0);
  
  // iteration
  for (int iter = 0; iter < (niter-1); iter ++) {

    // update gamma_k and theta_k
    for (int k = 0; k < K; k ++) {
      // generate A_k and S_hat_Gamma_k
      mat A_k(2,2); 
      A_k(0,0) = beta_cur[0] + alpha_cur[0] * eta_cur(k,0);
      A_k(1,1) = beta_cur[1] + alpha_cur[1] * eta_cur(k,1);
      mat S_hat_Gamma_k(2,2);
      S_hat_Gamma_k(0,0) = s2_hat_Gamma_1[k]; S_hat_Gamma_k(1,1) = s2_hat_Gamma_2[k];
      mat S_hat_gamma_k(2,2);
      S_hat_gamma_k(0,0) = s2_hat_gamma_1[k]; S_hat_gamma_k(1,1) = s2_hat_gamma_2[k];
      // update gamma_k
      mat precision_mat_post_gamma = trans(A_k) * inv(S_hat_Gamma_k) * A_k + inv(S_hat_gamma_k) + inv(Sigma_gamma_cur);
      vec temp1 = {Gamma_hat_1[k] - theta_cur(k,0), Gamma_hat_2[k] - theta_cur(k,1)};
      vec temp2 = {gamma_hat_1[k], gamma_hat_2[k]};
      vec mu_post_nom_gamma = trans(A_k) * inv(S_hat_Gamma_k) * temp1 + inv(S_hat_gamma_k) * temp2;
      vec mu_post_gamma = inv(precision_mat_post_gamma) * mu_post_nom_gamma;
      mat gamma_sample = mvrnormArma(1, mu_post_gamma, inv(precision_mat_post_gamma));
      gamma_cur(k,0) = gamma_sample(0,0);
      gamma_cur(k,1) = gamma_sample(0,1);
      // update theta_k
      mat precision_mat_post_theta = inv(S_hat_Gamma_k) + inv(Sigma_theta_cur);
      vec temp3 = {Gamma_hat_1[k] - beta_cur(0)*gamma_cur(k,0) - alpha_cur(0)*gamma_cur(k,0)*eta_cur(k,0), 
                   Gamma_hat_2[k] - beta_cur(1)*gamma_cur(k,1) - alpha_cur(1)*gamma_cur(k,1)*eta_cur(k,1)};
      vec mu_post_nom_theta = inv(S_hat_Gamma_k) * temp3;
      vec mu_post_theta = inv(precision_mat_post_theta) * mu_post_nom_theta;
      mat theta_sample = mvrnormArma(1, mu_post_theta, inv(precision_mat_post_theta));
      theta_cur(k,0) = theta_sample(0,0);
      theta_cur(k,1) = theta_sample(0,1);
    }
    
    //update Sigma_gamma
    mat S_gamma(2,2);
    for (int ii = 0; ii < K; ii ++) {
      rowvec rowk = gamma_cur.row(ii);
      S_gamma = S_gamma + trans(rowk) * rowk;
    }
    S_gamma = S_gamma / K; 
    double m_post_gamma = m_gamma + K;
    mat V_post_gamma = S_gamma * K + V_gamma;
    mat Sigma_gamma_sample = my_rinvwishart(m_post_gamma, V_post_gamma);
    double rho_gamma = Sigma_gamma_sample(1,0)/sqrt(Sigma_gamma_sample(0,0))/sqrt(Sigma_gamma_sample(3,0));
    if (Sigma_gamma_sample(0,0) < 1e-4) Sigma_gamma_sample(0,0) = 1e-4;
    if (Sigma_gamma_sample(0,0) > 10) Sigma_gamma_sample(0,0) = 10;
    if (Sigma_gamma_sample(3,0) < 1e-4) Sigma_gamma_sample(3,0) = 1e-4;
    if (Sigma_gamma_sample(3,0) > 10) Sigma_gamma_sample(3,0) = 10;
    Sigma_gamma_cur(0,0) = Sigma_gamma_sample(0,0);
    Sigma_gamma_cur(1,1) = Sigma_gamma_sample(3,0);
    Sigma_gamma_cur(0,1) = Sigma_gamma_sample(0,0)*Sigma_gamma_sample(3,0)*rho_gamma;
    Sigma_gamma_cur(1,0) = Sigma_gamma_cur(0,1);
    sigma2_gamma1_tk[iter+1] = Sigma_gamma_cur(0,0);
    sigma2_gamma2_tk[iter+1] = Sigma_gamma_cur(1,1);
    
    //update Sigma_theta
    mat S_theta(2,2);
    for (int ii = 0; ii < K; ii ++) {
      rowvec rowk = theta_cur.row(ii);
      S_theta = S_theta + trans(rowk) * rowk;
    }
    S_theta = S_theta / K; 
    double m_post_theta = m_theta + K;
    mat V_post_theta = S_theta * K + V_theta;
    mat Sigma_theta_sample = my_rinvwishart(m_post_theta, V_post_theta);
    double rho_theta = Sigma_theta_sample(1,0)/sqrt(Sigma_theta_sample(0,0))/sqrt(Sigma_theta_sample(3,0));
    if (Sigma_theta_sample(0,0) < 1e-8) Sigma_theta_sample(0,0) = 1e-8;
    if (Sigma_theta_sample(0,0) > 10) Sigma_theta_sample(0,0) = 10;
    if (Sigma_theta_sample(3,0) < 1e-8) Sigma_theta_sample(3,0) = 1e-8;
    if (Sigma_theta_sample(3,0) > 10) Sigma_theta_sample(3,0) = 10;
    Sigma_theta_cur(0,0) = Sigma_theta_sample(0,0);
    Sigma_theta_cur(1,1) = Sigma_theta_sample(3,0);
    Sigma_theta_cur(0,1) = Sigma_theta_sample(0,0)*Sigma_theta_sample(3,0)*rho_theta;
    Sigma_theta_cur(1,0) = Sigma_theta_cur(0,1);
    sigma2_theta1_tk[iter+1] = Sigma_theta_cur(0,0);
    sigma2_theta2_tk[iter+1] = Sigma_theta_cur(1,1);
        
    //update eta_k
    for (int k = 0; k < K; k ++) {
      double fk1_part1 = pow(Gamma_hat_1[k] - theta_cur(k,0) - (beta_cur[0] + alpha_cur[0]) * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k] * (-0.5);
      double fk1_part2 = pow(pst_cur[2], 1*(eta_cur(k,1) == 0)) * pow(pst_cur[3], 1*(eta_cur(k,1) == 1)); 
      double fk1 = exp(fk1_part1) * fk1_part2;
      double fk0_part1 = pow(Gamma_hat_1[k] - theta_cur(k,0) - beta_cur[0] * gamma_cur(k,0), 2) / s2_hat_Gamma_1[k] * (-0.5);
      double fk0_part2 = pow(pst_cur[0], 1*(eta_cur(k,1) == 0)) * pow(pst_cur[1], 1*(eta_cur(k,1) == 1)); 
      double fk0 = exp(fk0_part1) * fk0_part2;
      double q1_post = fk1/(fk1 + fk0);
      if (q1_post <= 0.1) { q1_post = 0; }
      if (q1_post >= 0.9) { q1_post = 1; }
      eta_cur(k,0) = rbinom(1, 1, q1_post)[0];
      q_post_1_cur[k] = q1_post;
      eta_1_tk(iter+1, k) = eta_cur(k,0);
      
      double gk1_part1 = pow(Gamma_hat_2[k] - theta_cur(k,1) - (beta_cur[1] + alpha_cur[1]) * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k] * (-0.5);
      double gk1_part2 = pow(pst_cur[1], 1*(eta_cur(k,0) == 0)) * pow(pst_cur[3], 1*(eta_cur(k,0) == 1)); 
      double gk1 = exp(gk1_part1) * gk1_part2;
      double gk0_part1 = pow(Gamma_hat_2[k] - theta_cur(k,1) - beta_cur[1] * gamma_cur(k,1), 2) / s2_hat_Gamma_2[k] * (-0.5);
      double gk0_part2 = pow(pst_cur[0], 1*(eta_cur(k,0) == 0)) * pow(pst_cur[2], 1*(eta_cur(k,0) == 1)); 
      double gk0 = exp(gk0_part1) * gk0_part2;
      double q2_post = gk1/(gk1 + gk0);
      if (q2_post <= 0.1) { q2_post = 0; }
      if (q2_post >= 0.9) { q2_post = 1; }
      eta_cur(k,1) = rbinom(1, 1, q2_post)[0];
      q_post_2_cur[k] = q2_post;
      eta_2_tk(iter+1, k) = eta_cur(k,1);
    }
    
    q_post_1_tk(iter+1, _) = q_post_1_cur;
    q_post_2_tk(iter+1, _) = q_post_2_cur;

    // update pst
    int n00 = 0; int n01 = 0; int n10 = 0; int n11 = 0;
    for (int k = 0; k < K; k ++) {
      if (eta_cur(k,0) == 0 && eta_cur(k,1) == 0) n00 = n00 + 1;
      if (eta_cur(k,0) == 0 && eta_cur(k,1) == 1) n01 = n01 + 1;
      if (eta_cur(k,0) == 1 && eta_cur(k,1) == 0) n10 = n10 + 1;
      if (eta_cur(k,0) == 1 && eta_cur(k,1) == 1) n11 = n11 + 1;
    }
    NumericVector dir_par = {n00 + cc[0], n01 + cc[1], n10 + cc[2], n11 + cc[3]};
    pst_cur = my_rdirichlet(1, dir_par);
    pst_tk(iter+1,0) = pst_cur[0];
    pst_tk(iter+1,1) = pst_cur[1];
    pst_tk(iter+1,2) = pst_cur[2];
    pst_tk(iter+1,3) = pst_cur[3];

    // Matrices used in updating alpha and beta
    mat Omega_hat_Gamma_1(K, K);  
    mat Omega_hat_Gamma_2(K, K); 
    for (int ii = 0; ii < K; ii ++) {
      Omega_hat_Gamma_1(ii, ii) = 1/s2_hat_Gamma_1[ii]; 
      Omega_hat_Gamma_2(ii, ii) = 1/s2_hat_Gamma_2[ii];     
    }
    
    // update alpha_j
    vec U_alpha_1(K);
    vec U_alpha_2(K);
    vec W_alpha_1(K);
    vec W_alpha_2(K);
    for (int ii = 0; ii < K; ii ++) {
      U_alpha_1[ii] = Gamma_hat_1[ii] - theta_cur(ii,0) - beta_cur[0] * gamma_cur(ii,0);
      U_alpha_2[ii] = Gamma_hat_2[ii] - theta_cur(ii,1) - beta_cur[1] * gamma_cur(ii,1);
      W_alpha_1[ii] = eta_cur(ii,0) * gamma_cur(ii,0);
      W_alpha_2[ii] = eta_cur(ii,1) * gamma_cur(ii,1);
    }
    
    if (sum(W_alpha_1) == 0) alpha_cur[0] = 0;
    else {
    mat A_alpha_1  = trans(W_alpha_1) * Omega_hat_Gamma_1 * W_alpha_1;
    mat mu_alpha_1 = trans(W_alpha_1) * Omega_hat_Gamma_1 * U_alpha_1 / A_alpha_1;
    //double alpha_sample_1 = rnorm(1, mu_alpha_1(0,0), sqrt(1/A_alpha_1(0,0)))[0];
    //if (alpha_sample_1 >  0) { alpha_cur[0] = std::min(alpha_sample_1,  10000.0); }
    //if (alpha_sample_1 <= 0) { alpha_cur[0] = std::max(alpha_sample_1, -10000.0); }
    alpha_cur[0] = rnorm(1, mu_alpha_1(0,0), sqrt(1/A_alpha_1(0,0)))[0];
    }
    
    if (sum(W_alpha_2) == 0) alpha_cur[1] = 0;
    else {
    mat A_alpha_2  = trans(W_alpha_2) * Omega_hat_Gamma_2 * W_alpha_2;
    mat mu_alpha_2 = trans(W_alpha_2) * Omega_hat_Gamma_2 * U_alpha_2 / A_alpha_2;
    //double alpha_sample_2 = rnorm(1, mu_alpha_2(0,0), sqrt(1/A_alpha_2(0,0)))[0];
    //if (alpha_sample_2 >  0) { alpha_cur[1] = std::min(alpha_sample_2,  10000.0); }
    //if (alpha_sample_2 <= 0) { alpha_cur[1] = std::max(alpha_sample_2, -10000.0); }
    alpha_cur[1] = rnorm(1, mu_alpha_2(0,0), sqrt(1/A_alpha_2(0,0)))[0];
    }
    
    alpha_1_tk[iter+1] = alpha_cur[0];
    alpha_2_tk[iter+1] = alpha_cur[1];
    
    // update beta_j
    vec U_beta_1(K);
    vec U_beta_2(K);
    vec W_beta_1(K);
    vec W_beta_2(K);
    for (int ii = 0; ii < K; ii ++) {
      U_beta_1[ii] = Gamma_hat_1[ii] - theta_cur(ii,0) - alpha_cur[0] * eta_cur(ii,0) * gamma_cur(ii,0);
      U_beta_2[ii] = Gamma_hat_2[ii] - theta_cur(ii,1) - alpha_cur[1] * eta_cur(ii,1) * gamma_cur(ii,1);
      W_beta_1[ii] = gamma_cur(ii,0);
      W_beta_2[ii] = gamma_cur(ii,1);
    }
    
    mat A_beta_1  = trans(W_beta_1) * Omega_hat_Gamma_1 * W_beta_1;
    mat A_beta_2  = trans(W_beta_2) * Omega_hat_Gamma_2 * W_beta_2;
    mat mu_beta_1 = trans(W_beta_1) * Omega_hat_Gamma_1 * U_beta_1 / A_beta_1;
    mat mu_beta_2 = trans(W_beta_2) * Omega_hat_Gamma_2 * U_beta_2 / A_beta_2;
    //double beta_sample_1 = rnorm(1, mu_beta_1(0,0), sqrt(1/A_beta_1(0,0)))[0];
    //double beta_sample_2 = rnorm(1, mu_beta_2(0,0), sqrt(1/A_beta_2(0,0)))[0];
    //if (beta_sample_1 >  0) { beta_cur[0] = std::min(beta_sample_1,  10.0); }
    //if (beta_sample_1 <= 0) { beta_cur[0] = std::max(beta_sample_1, -10.0); }
    //if (beta_sample_2 >  0) { beta_cur[1] = std::min(beta_sample_2,  10.0); }
    //if (beta_sample_2 <= 0) { beta_cur[1] = std::max(beta_sample_2, -10.0); }
    beta_cur[0] = rnorm(1, mu_beta_1(0,0), sqrt(1/A_beta_1(0,0)))[0];
    beta_cur[1] = rnorm(1, mu_beta_2(0,0), sqrt(1/A_beta_2(0,0)))[0];
    beta_1_tk[iter+1] = beta_cur[0];
    beta_2_tk[iter+1] = beta_cur[1];

  } // end iter
  
  List res = List::create(
    Named("K") = K, Named("alpha_1_tk") = alpha_1_tk, Named("alpha_2_tk") = alpha_2_tk, 
    Named("beta_1_tk") = beta_1_tk, Named("beta_2_tk") = beta_2_tk, 
    Named("eta_1_tk") = eta_1_tk, Named("eta_2_tk") = eta_2_tk, 
    Named("pst_tk") = pst_tk, 
    Named("q_post_1_tk") = q_post_1_tk, Named("q_post_2_tk") = q_post_2_tk,
    Named("sigma2_gamma1_tk") = sigma2_gamma1_tk, 
    Named("sigma2_gamma2_tk") = sigma2_gamma2_tk,
    Named("sigma2_theta1_tk") = sigma2_theta1_tk, 
    Named("sigma2_theta2_tk") = sigma2_theta2_tk);
  
  return res;
  
  //double dd = sigma2_beta_1_cur;
  
  //List res = List::create(Named("niter") = niter, Named("K") = K, Named("dd") = dd);
  //return res;
}

