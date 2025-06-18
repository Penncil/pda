
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double rcpp_cc_log_plk(NumericVector beta,
               List covariate_list,
               List failure_position,
               IntegerVector failure_num,
               List risk_sets,
               List risk_set_weights,
               int site_num) {
  
  NumericMatrix X = covariate_list[site_num - 1];
  IntegerVector fail_pos = as<IntegerVector>(failure_position[site_num - 1]);
  int fail_num = failure_num[site_num - 1];
  
  int n_row = X.nrow(), n_col = X.ncol();
  NumericVector eta(n_row);
  for (int i = 0; i < n_row; ++i) {
    for (int j = 0; j < n_col; ++j) {
      eta[i] += X(i, j) * beta[j];
    }
  }
  
  NumericVector exp_eta = exp(eta);
  double res = 0.0;
  for (int i = 0; i < fail_pos.size(); ++i) {
    res += eta[fail_pos[i] - 1];
  }
  
  List site_risk_sets = risk_sets[site_num - 1];
  List site_risk_weights = risk_set_weights[site_num - 1];
  for (int j = 0; j < fail_num; ++j) {
    IntegerVector idx = as<IntegerVector>(site_risk_sets[j]);
    NumericVector weights = as<NumericVector>(site_risk_weights[j]);
    double temp_sum = 0.0;
    for (int i = 0; i < idx.size(); ++i) {
      temp_sum += exp_eta[idx[i] - 1] * weights[i];
    }
    res -= log(temp_sum + 1e-12);
  }
  
  return res;
}

// [[Rcpp::export]]
double rcpp_cc_pool_fun(NumericVector beta,
                List covariate_list,
                List failure_position,
                IntegerVector failure_num,
                List risk_sets,
                List risk_set_weights,
                int K) {
  double total = 0.0;
  for (int i = 1; i <= K; ++i) {
    total += rcpp_cc_log_plk(beta, covariate_list, failure_position, failure_num, risk_sets, risk_set_weights, i);
  }
  return total;
}

// [[Rcpp::export]]
NumericVector rcpp_cc_grad_plk(NumericVector beta,
                       List covariate_list,
                       List failure_position,
                       IntegerVector failure_num,
                       List risk_sets,
                       List risk_set_weights,
                       int site_num) {
  
  NumericMatrix X = covariate_list[site_num - 1];
  IntegerVector fail_pos = as<IntegerVector>(failure_position[site_num - 1]);
  int p = X.ncol();
  int n = X.nrow();
  
  NumericVector eta(n, 0.0);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < p; ++j)
      eta[i] += X(i, j) * beta[j];
  NumericVector exp_eta = exp(eta);
  
  NumericVector grad(p, 0.0);
  for (int i = 0; i < fail_pos.size(); ++i) {
    int row = fail_pos[i] - 1;
    for (int j = 0; j < p; ++j)
      grad[j] += X(row, j);
  }
  
  List site_risk_sets = risk_sets[site_num - 1];
  List site_risk_weights = risk_set_weights[site_num - 1];
  
  for (int j = 0; j < failure_num[site_num - 1]; ++j) {
    IntegerVector idx = as<IntegerVector>(site_risk_sets[j]);
    NumericVector weights = as<NumericVector>(site_risk_weights[j]);
    
    double denom = 0.0;
    NumericVector weighted_sum(p, 0.0);
    
    for (int i = 0; i < idx.size(); ++i) {
      int row = idx[i] - 1;
      double w = weights[i] * exp_eta[row];
      denom += w;
      for (int k = 0; k < p; ++k)
        weighted_sum[k] += w * X(row, k);
    }
    
    for (int k = 0; k < p; ++k)
      grad[k] -= weighted_sum[k] / (denom + 1e-12);
  }
  
  return grad;
}


// [[Rcpp::export]]
NumericMatrix rcpp_cc_hess_plk(NumericVector beta,
                       List covariate_list,
                       List failure_position,
                       IntegerVector failure_num,
                       List risk_sets,
                       List risk_set_weights,
                       int site_num) {
  NumericMatrix X = covariate_list[site_num - 1];
  int p = X.ncol();
  int n_row = X.nrow();
  
  NumericVector eta(n_row);
  for (int i = 0; i < n_row; ++i) {
    for (int j = 0; j < p; ++j) {
      eta[i] += X(i, j) * beta[j];
    }
  }
  
  NumericVector exp_eta = exp(eta);
  NumericMatrix H(p, p);
  
  List site_risk_sets = risk_sets[site_num - 1];
  List site_risk_weights = risk_set_weights[site_num - 1];
  int fail_num = failure_num[site_num - 1];
  
  for (int j = 0; j < fail_num; ++j) {
    IntegerVector idx = as<IntegerVector>(site_risk_sets[j]);
    NumericVector weights = as<NumericVector>(site_risk_weights[j]);
    NumericVector weighted_sum(p);
    double denom = 0.0;
    
    for (int i = 0; i < idx.size(); ++i) {
      double w = exp_eta[idx[i] - 1] * weights[i];
      denom += w;
      for (int k = 0; k < p; ++k) {
        weighted_sum[k] += w * X(idx[i] - 1, k);
      }
    }
    
    NumericVector mu(p);
    for (int k = 0; k < p; ++k) {
      mu[k] = weighted_sum[k] / (denom + 1e-12);
    }
    
    for (int k = 0; k < p; ++k) {
      for (int l = 0; l < p; ++l) {
        double cov_term = 0.0;
        for (int i = 0; i < idx.size(); ++i) {
          double w = exp_eta[idx[i] - 1] * weights[i];
          cov_term += w * X(idx[i] - 1, k) * X(idx[i] - 1, l);
        }
        H(k, l) += mu[k] * mu[l] - cov_term / (denom + 1e-12);
      }
    }
  }
  
  return H;
}
