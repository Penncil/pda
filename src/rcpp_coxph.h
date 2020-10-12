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
#include "utils.h"

namespace Intsurv {
// define class for inputs and outputs
class RcppCoxph {
private:
  arma::vec time;           // observed times
  arma::vec event;          // event indicators
  arma::mat z;              // design matrix
  bool hasTies {false};     // if there exists ties on event times
  arma::uvec uni_event_ind; // the first unique index of event times
  arma::uvec uni_time_ind;  // the first unique index of times
  arma::uvec event_ind;     // indices of event times
  // arma::vec d_time;         // distinct event times
  arma::vec d_time0;         // distinct event times  (event times with event)
  arma::mat dz;            // design matrix multiple event
  arma::mat d_z;            // design matrix aggregated at d_time
  arma::vec delta_n;        // event counts at d_time
  arma::uvec s_time_ind;     // time order ind
    
public:
  double partial_logL {0}; // partial log-likelihood
  arma::vec coef;          // covariate coefficient estimates
  arma::mat h0;            // baseline hazard estimates
  
  // constructors
  RcppCoxph(const arma::vec time_,
            const arma::vec event_,
            const arma::mat z_)
  {
    // sort event and z based on time
    s_time_ind = arma::sort_index(time_);
    time = time_.elem(s_time_ind);
    event = event_.elem(s_time_ind);
    z = z_.rows(s_time_ind);
    dz = z;      // z multiple by event, 
    dz.rows(arma::find(event == 0)).zeros();
    
    event_ind = arma::find(event > 0);
    // check if there exists ties on time
    // hasTies = any_duplicated(time.elem(event_ind));
    hasTies = any_duplicated(time);
    uni_time_ind = find_unique(time);          
    
    if (hasTies) {
      // d_time0 = time.elem(uni_time_ind);
      // arma::vec d_time0 {time.elem(event_ind)};
      // d_time0 = time.elem(event_ind);
      // uni_time_ind = find_unique(time);
      // uni_event_ind = vec_intersection(uni_time_ind, event_ind);
      // d_time = time.elem(uni_event_ind);   
      // aggregate at distinct event times
      // delta_n = event.elem(event_ind);
      delta_n = aggregate_sum(event, time);      // dj, can be 0
      // d_z = z.rows(event_ind);
      d_z = aggregate_sum_cols(dz, time);         // sum_{Dj}X, can be 0
    } else {
      // arma::vec d_time0 {time.elem(event_ind)};
      // d_time0 = time.elem(event_ind);
      // d_time = time.elem(event_ind);
      delta_n = event;  // arma::ones(time.n_elem);  // wrong!!! 
      d_z = dz;                                // z.rows(event_ind);
    }
  }
  
  // function members with overloads
  
  //// For Breslow's formula:
  // function that computes objective function only
  inline double objective(const arma::vec& beta) const;
  
  // function that computes gradients only
  inline arma::vec gradient(const arma::vec& beta) const;
  
  // function that computes objective and overwrite gradients
  inline double objective(const arma::vec& beta, arma::vec& grad) const;
  
  // function that computes hessian (vectorized) only
  inline arma::vec hessian(const arma::vec& beta) const;
  
  // //// For Efron's tie correction formula:
  // // function that computes objective function only
  inline double objective_efron(const arma::vec& beta) const;
  // 
  // // function that computes gradients of Efron's tie correction only
  inline arma::vec gradient_efron(const arma::vec& beta) const;
  // 
  // // function that computes gradients of Efron's tie correction in distributed version
  inline arma::vec gradient_efron_dist(const arma::vec& beta, const arma::vec& ind_machine, bool& useLocal, const int dj_cutoff) const;
  
};


// the negative loglikelihood function based on the broslow's formula
inline double RcppCoxph::objective(const arma::vec& beta) const
{
  const arma::vec dz_beta {d_z * beta};                   // sum_{Dj}Xb
  const arma::vec exp_z_beta {arma::exp(z * beta)};       // eXb
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};   
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};     // sum_{Rj}eXb
  
  return - (arma::sum(dz_beta - delta_n % arma::log(h0_denom_event)));
}


// the negative gradient of negative loglikelihood function
inline arma::vec RcppCoxph::gradient(const arma::vec& beta) const
{
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};     // sum_{Rj}eXb
  
  arma::mat numer_mat {arma::zeros(z.n_rows, z.n_cols)};
  for (size_t i {0}; i < z.n_rows; ++i) {
    numer_mat.row(i) = exp_z_beta(i) * z.row(i);            // XeXb
  }
  numer_mat = cum_sum_cols(numer_mat,  true);      
  numer_mat = numer_mat.rows(uni_time_ind);                // sum_{Rj}XeXb
  for (size_t j {0}; j < z.n_cols; ++j) {
    numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;       // dj * sum_{Rj}XeXb / sum_{Rj}eXb
  }
  // if (hasTies) {
  //   arma::vec h0_denom_event {h0_denom.elem(uni_event_ind)};
  //   numer_mat = numer_mat.rows(uni_event_ind);
  //   for (size_t j {0}; j < z.n_cols; ++j) {
  //     numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;
  //   }
  //   return - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
  // } else {
  //   arma::vec h0_denom_event {h0_denom.elem(event_ind)};
  //   numer_mat = numer_mat.rows(event_ind);
  //   for (size_t j {0}; j < z.n_cols; ++j) {
  //     numer_mat.col(j) = numer_mat.col(j) / h0_denom_event;
  //   }
  //   return - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
  // }
  return - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
}

// the negative loglikelihood function based on the broslow's formula
inline double RcppCoxph::objective(const arma::vec& beta,
                                   arma::vec& grad) const
{
  const arma::vec dz_beta {d_z * beta};
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta, true)};
  arma::mat numer_mat {arma::zeros(z.n_rows, z.n_cols)};
  for (size_t i {0}; i < z.n_rows; ++i) {
    numer_mat.row(i) = exp_z_beta(i) * z.row(i);
  }
  numer_mat = cum_sum_cols(numer_mat, true);
  if (hasTies) {
    arma::vec h0_denom_event {h0_denom.elem(uni_event_ind)};
    numer_mat = numer_mat.rows(uni_event_ind);
    for (size_t j {0}; j < z.n_cols; ++j) {
      numer_mat.col(j) = numer_mat.col(j) % delta_n / h0_denom_event;
    }
    // overwrite grad
    grad = - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
    return - arma::sum(dz_beta - delta_n % arma::log(h0_denom_event));
  } else {
    arma::vec h0_denom_event {h0_denom.elem(event_ind)};
    numer_mat = numer_mat.rows(event_ind);
    for (size_t j {0}; j < z.n_cols; ++j) {
      numer_mat.col(j) = numer_mat.col(j) / h0_denom_event;
    }
    // overwrite grad
    grad = - (arma::sum(d_z, 0) - arma::sum(numer_mat, 0)).t();
    return - (arma::sum(dz_beta) -
              arma::sum(arma::log(h0_denom_event)));
  }
}

// the hessian matrix of negative loglikelihood function
inline arma::vec RcppCoxph::hessian(const arma::vec& beta) const
{
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};     // sum_{Rj}eXb
  
  arma::mat numer_mat {arma::zeros(z.n_rows, z.n_cols)};
  arma::mat numer_mat2 {arma::zeros(z.n_rows, z.n_cols*z.n_cols)};
  // std::cout << "chkpt-0 " << size(numer_mat2) << vectorise(z.row(1).t() * z.row(1));
  for (size_t i {0}; i < z.n_rows; ++i) {
    numer_mat.row(i) = exp_z_beta(i) * z.row(i);                          // XeXb
    numer_mat2.row(i) = exp_z_beta(i) * vectorise(z.row(i).t() * z.row(i)).t();   // XX'eXb
  }
  // std::cout << "chkpt-1 " << numer_mat2.row(1);
  numer_mat = cum_sum_cols(numer_mat,  true);      
  numer_mat = numer_mat.rows(uni_time_ind);                   // sum_{Rj}XeXb
  numer_mat2 = cum_sum_cols(numer_mat2,  true);      
  numer_mat2 = numer_mat2.rows(uni_time_ind);                 // sum_{Rj}XX'eXb
  // std::cout << "chkpt-2 " << numer_mat2.row(1);
  
  for (size_t j {0}; j < z.n_cols; ++j) {
    numer_mat.col(j) = numer_mat.col(j) % sqrt(delta_n) / h0_denom_event;       // sqrt(dj) * sum_{Rj}XeXb / sum_{Rj}eXb
  }
  for (size_t j {0}; j < numer_mat2.n_cols; ++j) {
    numer_mat2.col(j) = numer_mat2.col(j) % delta_n / h0_denom_event;       // dj * sum_{Rj}XeXb / sum_{Rj}eXb
  }
  
  // std::cout << "chkpt-3 " << numer_mat2.row(1);
  arma::mat mat3 {arma::zeros(numer_mat.n_rows, z.n_cols*z.n_cols)};
  for (size_t i {0}; i < numer_mat.n_rows; ++i) {
    mat3.row(i) = vectorise(numer_mat.row(i).t() * numer_mat.row(i)).t();                           
  }
  // std::cout << "chkpt-4 " << numer_mat2.row(1) << mat3.row(1);
  // vec V(X.memptr(), X.n_elem, false, false);
  
  return arma::sum(numer_mat2 - mat3, 0).t();
}




// // the negative loglikelihood function based on the Efron's formula
inline double RcppCoxph::objective_efron(const arma::vec& beta) const
{
  // const arma::vec dz_beta {dz * beta};
  const arma::vec d_z_beta {d_z * beta};
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};        
  arma::vec exp_dz_beta {exp_z_beta};     
  exp_dz_beta.elem(arma::find(event == 0)).zeros();
  exp_dz_beta = aggregate_sum(exp_dz_beta, time);                       // sum_{Dj}eXb
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};         // sum_{Rj}eXb

  arma::vec efron {arma::zeros(d_z.n_rows)};
  for (size_t j {0}; j < d_z.n_rows; ++j) {
    double tmp {0};
    for (size_t h {0}; h < delta_n(j); ++h){
      if(delta_n(j) > 0) tmp = tmp + log(h0_denom_event(j) - h*exp_dz_beta(j) / delta_n(j));    // += log(sum_{Rj}eXb - h * sum_{Dj}eXb / dj)
    }
    efron(j) = tmp;
  }

  return - (arma::sum(d_z_beta) - arma::sum(efron));
}

// // // the negative gradient of negative loglikelihood function
inline arma::vec RcppCoxph::gradient_efron(const arma::vec& beta) const
{
  // const arma::vec dz_beta {dz * beta};
  // const arma::vec d_z_beta {d_z * beta};
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};        
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};         // sum_{Rj}eXb
  arma::vec exp_dz_beta {exp_z_beta};     
  exp_dz_beta.elem(arma::find(event == 0)).zeros();
  exp_dz_beta = aggregate_sum(exp_dz_beta, time);                       // sum_{Dj}eXb
  
  arma::mat Z_mat {arma::zeros(z.n_rows, z.n_cols)};
  for (size_t i {0}; i < z.n_cols; ++i) Z_mat.col(i) = exp_z_beta % z.col(i);    // XeXb
  
  arma::mat RZeZb {cum_sum_cols(Z_mat,  true)};            
  RZeZb = RZeZb.rows(uni_time_ind);                                     // sum_{Rj}XeXb
  arma::mat DZeZb {Z_mat};    
  DZeZb.rows(arma::find(event == 0)).zeros();
  DZeZb = aggregate_sum_cols(DZeZb, time);                                   // sum_{Dj}XeXb 

  arma::mat efron {arma::zeros(d_z.n_rows, z.n_cols)};
  for (size_t j {0}; j < d_z.n_rows; ++j) {
    arma::vec tmp {arma::zeros(z.n_cols)};
    for (size_t h {0}; h < delta_n(j); ++h){            // += (sum_{Rj}XeXb - h * (sum_{Dj}XeXb) / dj) / (sum_{Rj}eXb - h * (sum_{Dj}eXb) / dj)
      if(delta_n(j) > 0) tmp = tmp + (RZeZb.row(j) - h*DZeZb.row(j) / delta_n(j)).t()  / (h0_denom_event(j) - h*exp_dz_beta(j) / delta_n(j));
    }
    efron.row(j) = tmp.t() ;
  }

  return - arma::sum(d_z-efron, 0).t();
}


// // the negative gradient of negative loglikelihood function, using efron's formula
// distributed in K machines, with denom weighted by M1 and Mk
inline arma::vec RcppCoxph::gradient_efron_dist(const arma::vec& beta,
                                                const arma::vec& ind_machine_,
                                                bool& useLocal,
                                                const int dj_cutoff=0) const
{
  arma::vec uni_machine {arma::unique(ind_machine_)};
  // for now assume the 1st machine is the largest machine / local machine
  arma::uword K {uni_machine.n_elem};
  arma::vec ind_machine = ind_machine_.elem(s_time_ind);    // ind_machine sorted by increasing time 
  
  const arma::vec uni_time {arma::unique(time)};
  const arma::vec exp_z_beta {arma::exp(z * beta)};
  const arma::vec h0_denom {cum_sum(exp_z_beta,  true)};        
  const arma::vec h0_denom_event {h0_denom.elem(uni_time_ind)};         // sum_{Rj}eXb
  arma::vec exp_dz_beta {exp_z_beta};     
  exp_dz_beta.elem(arma::find(event == 0)).zeros();
  exp_dz_beta = aggregate_sum(exp_dz_beta, time);                       // sum_{Dj}eXb
  
  arma::mat Z_mat {arma::zeros(z.n_rows, z.n_cols)};
  for (size_t i {0}; i < z.n_cols; ++i) Z_mat.col(i) = exp_z_beta % z.col(i);    // XeXb
  
  arma::mat RZeZb {cum_sum_cols(Z_mat,  true)};            
  RZeZb = RZeZb.rows(uni_time_ind);                                     // sum_{Rj}XeXb, used in denom
  arma::mat DZeZb {Z_mat};
  DZeZb.rows(arma::find(event == 0)).zeros();
  DZeZb = aggregate_sum_cols(DZeZb, time);                              // sum_{Dj}XeXb

  // get sum_Dj0_eZb from the local (first) machine
  // this will be used by other machines
  arma::vec efron {arma::zeros(z.n_cols)};            // cumulate the Efron terms through: (k, j, h)
  arma::vec delta_n0, exp_dz_beta0, uni_time0;        // store local machine and send to other machines for calculating gradients
 
  for(size_t k {0}; k < K; ++k){
    // std::cout << "k=" << k;
    // get data of each machine and calculate the aggregated quantities
    // notice ind_machine is already sorted by time (s_time_ind)
    arma::uvec idx {arma::find(ind_machine==k)};
    //if(accu(idx)==0) std::cout << "machine not exist: " << k << std::endl;
    arma::mat zk  {z.rows(idx)};
    arma::vec timek  {time.elem(idx)};
    arma::vec uni_timek {arma::unique(timek)};               // uni_timek may be smaller than uni_time
    arma::vec eventk  {event.elem(idx)};
    arma::vec delta_nk = aggregate_sum(eventk, timek);      // djk, can be 0
    arma::uvec uni_time_indk {find_unique(timek)};
    arma::mat d_zk = zk;                                     // zk multiple by eventk, 
    d_zk.rows(arma::find(eventk == 0)).zeros();
    d_zk = aggregate_sum_cols(d_zk, timek);
    // quantities involving beta
    arma::vec exp_z_betak {arma::exp(zk * beta)};
    arma::vec exp_dz_betak {exp_z_betak};     
    exp_dz_betak.elem(arma::find(eventk == 0)).zeros();
    exp_dz_betak = aggregate_sum(exp_dz_betak, timek);                        // sum_{Djk}eXb
    
    arma::mat Z_matk {arma::zeros(zk.n_rows, zk.n_cols)};
    for (size_t i {0}; i < zk.n_cols; ++i) Z_matk.col(i) = exp_z_betak % zk.col(i);      // XexpXb
    arma::mat RZeZbk {cum_sum_cols(Z_matk,  true)};            
    RZeZbk = RZeZbk.rows(uni_time_indk);                                      // sum_{Rjk}XeXb
    arma::mat DZeZbk {Z_matk};    
    DZeZbk.rows(arma::find(eventk == 0)).zeros();
    DZeZbk = aggregate_sum_cols(DZeZbk, timek);                               // sum_{Djk}XeXb 
   
    if(k==0) {                        // store local machine and send to other machines for calculating gradients
      delta_n0 = delta_nk;            // djk
      exp_dz_beta0 = exp_dz_betak;    // sum_{Djk}eXb
      uni_time0 = uni_timek;
    }
    
    int jk, j0, djk, dj0;     //
    arma::vec num1k=arma::zeros(z.n_cols), num2k=arma::zeros(z.n_cols);  // , djkz=arma::zeros(z.n_cols)
    double denom2k = 0;
    for (size_t jj {0}; jj < d_z.n_rows; ++jj) {
      int j = d_z.n_rows - jj - 1;
      if(delta_n(j) > 0) {                      // no contribution to gradient if tj has no event  
        if(any( uni_timek == uni_time(j) )){    // the jth uni_time  exist in machine k, def num1k, num2k, denom2k
          jk = arma::conv_to<int>::from(arma::find( uni_timek == uni_time(j), 1));
          djk = delta_nk(jk);
          num1k = RZeZbk.row(jk).t();
          num2k = DZeZbk.row(jk).t();
          
          // djkz = d_zk.row(jk).t();
          // approx sum_{Dj}eXb/dj by sum_{Djk}eXb/djk and (if useLocal) sum_{Dj0}eXb/dj0
          if(useLocal){                         // denom2k: use machine 0 + machine k to approx sum_{Dj}eXb / dj at machine k
            if(any(uni_time0==uni_time(j))) {
              j0 = arma::conv_to<int>::from(arma::find(uni_time0==uni_time(j), 1));
              dj0 =  delta_n0(j0);               
            } else{
              j0 = -1;
              dj0 = 0;
            }
            // double tmp = delta_nk.elem(jk)+delta_n0.elem(j0);
            if(dj0 > 0)
              denom2k = (exp_dz_betak(jk)+exp_dz_beta0(j0)) / (djk+dj0);
            else  
              denom2k = (exp_dz_betak(jk)) / djk;   // tjk not exist in uni_time0, can't use Local machine 0
          } else {
              denom2k =  (exp_dz_betak(jk)) / djk;  // only use remote machine k to approx sum_{Dj}eXb / dj at machine k
          }
          
          // if(arma::is_finite(denom2k)==false) denom2k = 0;    // this happens if djk=0
          if(djk <= dj_cutoff){                         // if djk too small, avoid using Efron's ...
            if(useLocal==false || dj0 <= dj_cutoff){
              denom2k = 0;
              num2k.fill(0);
            }
          }
          
        } else{                                     // the jth uni_time not exist in machine k
          // num2k = sum_{Djk}XeXb = 0, but num1k = sum_{Rjk}XeXb and denom2k still exist! use values from last step
          num2k = arma::zeros(z.n_cols);            // 
          // out.row(j).subvec(2*k+1+2*K+4, 2*k+2+2*K+4) = num1k.t();
          // djkz = arma::zeros(z.n_cols);
        } // end def num1k, num2k, denom2k
        
        // arma::vec tmp = arma::zeros(z.n_cols);
        for (size_t h {0}; h < delta_n(j); ++h){     // h from 0 to dj-1 !!!
          // tmp = tmp + (num1k - h*num2k/delta_n(j)) / (h0_denom_event(j) - h*denom2k);
          arma::vec tmp = (num1k - h*num2k/delta_n(j)) / (h0_denom_event(j) - h*denom2k);
          if(tmp.is_finite()==false) {
            // std::cout << "NaN efron term ! (k, j)=" << k << j <<  "denom2k=" << denom2k << std::endl;
            tmp.fill(0);
          }
          efron = efron + tmp;
        }
      }   // end if(delta_n(j) > 0)
    }     // end of j
  }       // end of k

  return - (arma::sum(d_z, 0).t() - efron);
}

}       // end namespace Intsurv
