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


// Gibbs sampling

// [[Rcpp::export]]
List gibbs_rcpp_nopr(int niter, int K, NumericVector beta_tk, NumericVector alpha_tk, 
                         NumericMatrix eta_tk, NumericMatrix theta_tk, NumericMatrix gamma_tk,
                         NumericVector q_tk,
                         NumericVector Gamma_hat, NumericVector gamma_hat, 
                         vec s2_hat_Gamma, NumericVector s2_hat_gamma, 
                         NumericVector sigma2_gamma_tk, NumericVector sigma2_theta_tk,
                         double a_gamma, double b_gamma, double a_theta, double b_theta,
                         double a_q, double b_q) {
  
  // set current values as initial values
  double beta_cur = beta_tk[0];
  double alpha_cur = alpha_tk[0];
  NumericVector theta_cur = theta_tk(0, _);
  NumericVector gamma_cur = gamma_tk(0, _);
  NumericVector eta_cur = eta_tk(0, _);
  double sigma2_gamma_cur = sigma2_gamma_tk[0];
  double sigma2_theta_cur = sigma2_theta_tk[0];
  double q_cur = q_tk[0];
  NumericMatrix q_post_tk(niter, K);
  NumericVector q_post_cur = q_post_tk(0, _);

  
  // iteration
  for (int iter = 0; iter < (niter-1); iter ++) {

    // update gamma_k
    for (int k = 0; k < K; k ++) {
      double Ak_gamma = pow(beta_cur + alpha_cur * eta_cur[k], 2)/s2_hat_Gamma[k] + 1/s2_hat_gamma[k] + 1/sigma2_gamma_cur;
      double Bk_gamma = (beta_cur + alpha_cur * eta_cur[k])*(Gamma_hat[k] - theta_cur[k])/s2_hat_Gamma[k] + gamma_hat[k]/s2_hat_gamma[k];
      gamma_cur[k] = rnorm(1, Bk_gamma/Ak_gamma, sqrt(1/Ak_gamma))[0];
    }
    
    gamma_tk(iter+1, _) = gamma_cur;
    
    // update theta_k
    for (int k = 0; k < K; k ++) {
      double Ak_theta = 1/s2_hat_Gamma[k] + 1/sigma2_theta_cur;
      double Bk_theta = (Gamma_hat[k] - (beta_cur + alpha_cur * eta_cur[k])*gamma_cur[k])/s2_hat_Gamma[k];
      theta_cur[k] = rnorm(1, Bk_theta/Ak_theta, sqrt(1/Ak_theta))[0];
    }
    
    theta_tk(iter+1, _) = theta_cur;
    
    // update sigma2_gamma and sigma2_theta
    sigma2_gamma_cur = my_rinvgamma(1, a_gamma + K/2, b_gamma + 0.5*sum(pow(gamma_cur, 2)))[0];
    sigma2_theta_cur = my_rinvgamma(1, a_theta + K/2, b_theta + 0.5*sum(pow(theta_cur, 2)))[0];
    sigma2_gamma_tk[iter+1] = sigma2_gamma_cur;
    sigma2_theta_tk[iter+1] = sigma2_theta_cur;
    
    // update eta_k
    for (int k = 0; k < K; k ++) {
      double fk1_part1 = pow(Gamma_hat[k] - theta_cur[k] - (beta_cur + alpha_cur) * gamma_cur[k], 2) / s2_hat_Gamma[k] * (-0.5);
      double fk1_part2 = q_cur; 
      double fk1 = exp(fk1_part1) * fk1_part2;
      double fk0_part1 = pow(Gamma_hat[k] - theta_cur[k] - beta_cur * gamma_cur[k], 2) / s2_hat_Gamma[k] * (-0.5);
      double fk0_part2 = 1 - q_cur;
      double fk0 = exp(fk0_part1) * fk0_part2;
      double q_post = fk1/(fk1 + fk0);
      if (q_post <= 0.1) { 
        q_post = 0;
      }
      if (q_post >= 0.9) { 
        q_post = 1; 
      }
      
      eta_cur[k] = rbinom(1, 1, q_post)[0];
      q_post_cur[k] = q_post;
    }
    
    eta_tk(iter+1, _) = eta_cur;
    q_post_tk(iter+1, _) = q_post_cur;
        
    // update q
    int n0 = sum(eta_cur == 0);
    int n1 = sum(eta_cur == 1);
    q_cur = rbeta(1, n1 + a_q, n0 + b_q)[0];
    q_tk[iter+1] = q_cur;
      
    // Matrices used in updating alpha and beta
    mat Omega_hat_Gamma(K, K);  
    for (int ii = 0; ii < K; ii ++) {
      Omega_hat_Gamma(ii, ii) = 1/s2_hat_Gamma[ii];       
    }
    
    // update alpha
    vec U_alpha = Gamma_hat - theta_cur - beta_cur * gamma_cur;
    vec W_alpha = eta_cur * gamma_cur;
    if (sum(W_alpha) == 0) {
    alpha_cur = 0;
    } 
    else {
    //if (sum(W_alpha) == 0) {
    mat A_alpha  = trans(W_alpha) * Omega_hat_Gamma * W_alpha;
    mat mu_alpha = trans(W_alpha) * Omega_hat_Gamma * U_alpha / A_alpha;
    alpha_cur = rnorm(1, mu_alpha(0,0), sqrt(1/A_alpha(0,0)))[0];
    //if (alpha_cur >  0) { alpha_cur = std::min(alpha_cur,  10.0); }
    //if (alpha_cur <= 0) { alpha_cur = std::max(alpha_cur, -10.0); }
    }
    alpha_tk[iter+1] = alpha_cur;
    
    // update beta
    vec U_beta = Gamma_hat - theta_cur - alpha_cur * eta_cur * gamma_cur;
    vec W_beta = gamma_cur;
    mat A_beta  = trans(W_beta) * Omega_hat_Gamma * W_beta;
    mat mu_beta = trans(W_beta) * Omega_hat_Gamma * U_beta / A_beta; 
    beta_cur = rnorm(1, mu_beta(0,0), sqrt(1/A_beta(0,0)))[0];
    //if (beta_cur >  0) { beta_cur = std::min(beta_cur,  10.0); }
    //if (beta_cur <= 0) { beta_cur = std::max(beta_cur, -10.0); }
    beta_tk[iter+1] = beta_cur;
    
  } // end iter
  
  List res = List::create(
    Named("K") = K, Named("alpha_tk") = alpha_tk, 
    Named("beta_tk") = beta_tk, Named("eta_tk") = eta_tk, Named("q_tk") = q_tk,
    Named("gamma_tk") = gamma_tk, Named("theta_tk") = theta_tk,
    Named("sigma2_gamma_tk") = sigma2_gamma_tk, Named("sigma2_theta_tk") = sigma2_theta_tk
    ,
    Named("q_post_tk") = q_post_tk
    );
  
  return res;
  
  
  //List res = List::create(Named("niter") = niter, Named("K") = K, Named("cc") = cc);
  //return res;
}


/*** R
# gibbs_rcpp(niter, K, beta_tk, alpha_tk, eta_tk, theta_tk, gamma_tk, q_tk, Gamma_hat, gamma_hat,
#           s2_hat_Gamma, s2_hat_gamma, sigma2_gamma_tk, sigma2_theta_tk, sigma2_alpha_tk,
#           sigma2_beta_tk, a_gamma, b_gamma, a_theta, b_theta,
#           a_alpha, b_alpha, a_beta, b_beta,
#           a_q, b_q)
*/

