#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double rcpp_cc_log_plk(NumericVector beta,
               NumericMatrix covariate,
               IntegerVector failure_position,
               int failure_num,
               List risk_sets) {
  // log partial Lik, no negate or average
  int n = covariate.nrow();
  int p = covariate.ncol();
  
  NumericVector eta(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
      acc += covariate(i, j) * beta[j];
    eta[i] = acc;
  }
  
  NumericVector exp_eta = exp(eta);
  
  double res = 0.0;
  int nf = failure_position.size();
  for (int i = 0; i < nf; ++i) {
    int row = failure_position[i] - 1;
    res += eta[row];
  }
  
  for (int j = 0; j < failure_num; ++j) {
    IntegerVector idx = as<IntegerVector>(risk_sets[j]);
    double s = 0.0;
    for (int k = 0; k < idx.size(); ++k)
      s += exp_eta[idx[k] - 1];
    res -= std::log(s + 1e-12);
  }
  
  return res;
}

// // [[Rcpp::export]]
// double rcpp_cc_pool_fun(NumericVector beta,
//                 List covariate_list,
//                 List failure_position,
//                 IntegerVector failure_num,
//                 List risk_sets,
//                 List risk_set_weights,
//                 int K) {
//   double total = 0.0;
//   for (int i = 1; i <= K; ++i) {
//     total += rcpp_cc_log_plk(beta, covariate_list, failure_position, failure_num, risk_sets, risk_set_weights, i);
//   }
//   return total;
// }

// [[Rcpp::export]]
NumericVector rcpp_cc_grad_plk(NumericVector beta,
                               NumericMatrix X,
                               IntegerVector failure_position,
                               int failure_num,
                               List risk_sets) {
  
  int n = X.nrow();
  int p = X.ncol();
  
  NumericVector eta(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
      acc += X(i, j) * beta[j];
    eta[i] = acc;
  }
  NumericVector exp_eta = exp(eta);
  
  NumericVector grad(p, 0.0);
  int nf = failure_position.size();
  for (int i = 0; i < nf; ++i) {
    int row = failure_position[i] - 1;
    for (int j = 0; j < p; ++j)
      grad[j] += X(row, j);
  }
  
  for (int j = 0; j < failure_num; ++j) {
    IntegerVector idx = as<IntegerVector>(risk_sets[j]);
    
    double denom = 0.0;
    NumericVector weighted_sum(p, 0.0);
    
    for (int k = 0; k < idx.size(); ++k) {
      int row = idx[k] - 1;
      double w = exp_eta[row];
      denom += w;
      for (int t = 0; t < p; ++t)
        weighted_sum[t] += w * X(row, t);
    }
    
    for (int t = 0; t < p; ++t)
      grad[t] -= weighted_sum[t] / denom;
  }
  
  return grad;
}


// [[Rcpp::export]]
NumericMatrix rcpp_cc_hess_plk(NumericVector beta,
                               NumericMatrix X,
                               int failure_num,
                               List risk_sets) {
  
  int n = X.nrow();
  int p = X.ncol();
  
  NumericVector eta(n, 0.0);
  for (int i = 0; i < n; ++i) {
    double acc = 0.0;
    for (int j = 0; j < p; ++j)
      acc += X(i, j) * beta[j];
    eta[i] = acc;
  }
  NumericVector exp_eta = exp(eta);
  
  NumericMatrix H(p, p);
  
  for (int j = 0; j < failure_num; ++j) {
    IntegerVector idx = as<IntegerVector>(risk_sets[j]);
    
    double denom = 0.0;
    NumericVector mean_vec(p, 0.0);
    
    for (int k = 0; k < idx.size(); ++k) {
      int row = idx[k] - 1;
      double w = exp_eta[row];
      denom += w;
      for (int t = 0; t < p; ++t)
        mean_vec[t] += w * X(row, t);
    }
    
    NumericMatrix cov_mat(p, p);
    for (int k = 0; k < idx.size(); ++k) {
      int row = idx[k] - 1;
      double w = exp_eta[row];
      for (int a = 0; a < p; ++a) {
        double xa = X(row, a);
        for (int b = 0; b < p; ++b)
          cov_mat(a, b) += w * xa * X(row, b);
      }
    }
    
    double denom2 = denom * denom;
    
    for (int a = 0; a < p; ++a) {
      for (int b = 0; b < p; ++b) {
        double term1 = (mean_vec[a] * mean_vec[b]) / denom2;
        double term2 = cov_mat(a, b) / denom;
        H(a, b) += term1 - term2;
      }
    }
  }
  
  return H;
}
