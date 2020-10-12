// Copyright (2020) Chongliang Luo, Rui Duan, Jiayi Tong and Yong Chen
//     
//     This file is part of pda
//     
//     Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//     
//     http://www.apache.org/licenses/LICENSE-2.0
//     
//     Unless required by applicable law or agreed to in writing, software
//     distributed under the License is distributed on an "AS IS" BASIS,
//     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
//     limitations under the License.


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include "rcpp_coxph.h"

// Cox partial (neg) log-likelihood
// [[Rcpp::export]]
double rcpp_coxph_logL(const arma::vec& beta,
               const arma::vec& time,
               const arma::vec& event,
               const arma::mat& z)
{
  return Intsurv::RcppCoxph(time, event, z).objective(beta);
}

// gradient
// [[Rcpp::export]]
arma::vec rcpp_coxph_logL_gradient(const arma::vec& beta,
               const arma::vec& time,
               const arma::vec& event,
               const arma::mat& z)
{
  return Intsurv::RcppCoxph(time, event, z).gradient(beta);
}

// [[Rcpp::export]]
arma::vec rcpp_coxph_logL_hessian(const arma::vec& beta,
                                  const arma::vec& time,
                                  const arma::vec& event,
                                  const arma::mat& z)
{
  return Intsurv::RcppCoxph(time, event, z).hessian(beta);
}


// use efron's approximation for tie correction
// [[Rcpp::export]]
double rcpp_coxph_logL_efron(const arma::vec& beta,
                                    const arma::vec& time,
                                    const arma::vec& event,
                                    const arma::mat& z)
{
   return Intsurv::RcppCoxph(time, event, z).objective_efron(beta);
}
 
// gradient with efron's
// [[Rcpp::export]]
arma::vec rcpp_coxph_logL_gradient_efron(const arma::vec& beta,
                                   const arma::vec& time,
                                   const arma::vec& event,
                                   const arma::mat& z)
{
  return Intsurv::RcppCoxph(time, event, z).gradient_efron(beta);
}


// gradient with efron's with distributed machines, 
// assuming limited communication of sum_Djs_exp(Xb), sum_Rjs_exp(Xb), sum_Djs_Xexp(Xb) and sum_Rjs_Xexp(Xb) ...
// [[Rcpp::export]]
arma::vec rcpp_coxph_logL_gradient_efron_dist(const arma::vec& beta,
                                              const arma::vec& ind_machine_, 
                                              bool& useLocal,        // transfer quantities from Local to remote to help approx the denominator of the likelihood?
                                              const int dj_cutoff,   // if djs < cutoff, don't use Efron's for that machine that time
                                              const arma::vec& time,
                                              const arma::vec& event,
                                              const arma::mat& z)
{
  return Intsurv::RcppCoxph(time, event, z).gradient_efron_dist(beta, ind_machine_, useLocal, dj_cutoff);
}


// [[Rcpp::export]]
arma::mat rcpp_aggregate(const arma::mat& x,
                              const arma::vec& indices,
                              const bool simplify = true,
                              const bool cumulative = false,
                              const bool reversely = false)
{
  return Intsurv::aggregate_sum_cols(x, indices, simplify, cumulative, reversely);
}
