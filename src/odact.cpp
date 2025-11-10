#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp; 
using namespace std;
 
// cumulative sum in possibly reverse order
// [[Rcpp::export]]
arma::vec cum_sum(arma::vec& x, bool reversely = false)
{
  // if cumsum reversely
  if (reversely) {
    int n_x {static_cast<int>(x.n_rows)};
    arma::vec res {arma::zeros(n_x)};
    double tmp {0.0};
    for (size_t i {1}; i <= n_x; ++i) {
      tmp += x[n_x - i];
      res[n_x - i] = tmp;
    }
    return res;
  }
  // otherwise, using arma::cumsum
  return arma::cumsum(x);
}

 
// column-wise cumulative sum in possibly reverse order
// [[Rcpp::export]]
arma::mat cum_sum_cols(arma::mat& x, bool reversely = false)
{
  // if cumsum reversely
  if (reversely) {
    int n_x = {static_cast<int>(x.n_rows)};
    arma::mat tmp {arma::zeros(1, x.n_cols)};
    arma::mat res {x};
    for (size_t i {1}; i <= n_x; ++i) {
      tmp += x.row(n_x - i);
      res.row(n_x - i) = tmp;
    }
    return res;
  }
  // otherwise, using arma::cumsum
  return arma::cumsum(x, 0);
}

 

// [[Rcpp::export]]
arma::mat Xotimes2(arma::mat& X){
  // t(apply(X, 1, function(x) c(x %o% x )))
  int px = {static_cast<int>(X.n_cols)}, n = {static_cast<int>(X.n_rows)}, i;   
  arma::mat X2(n, px*px, arma::fill::zeros), m; 
  arma::rowvec Xi; 
  
  for(i=0; i<n; i++){
    Xi = X.row(i);
    m = Xi.t() * Xi;
    X2.row(i) = m.as_row(); 
  }
  
  return X2;
}
