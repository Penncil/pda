//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace Eigen; 
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double soft_c(double x, double thres) {
  double out;
  if (x > thres) {
    out=x-thres;
  } else if(x<-thres){
    out=x+thres;
  } else{
    out=0;
  }
  return out;
}


// [[Rcpp::export]]
List coordi_cho(VectorXd atilde, MatrixXd B, VectorXd betainit, double lambda, int iter_max) {
  int iter=1, p=betainit.size();
  double dif=1, a, b, c;
  VectorXd beta1;
  while(iter<=iter_max && dif>1e-4){
    beta1=betainit;
    for(int i=0; i<p; i++){
      a=B(i,i);
      c=betainit.adjoint()*B.col(i);
      b=atilde(i) + c - a*betainit(i);
      if(i==0){
        betainit(i)=-b/a;
      }else{
      betainit(i)=-soft_c(b/a,lambda/a);
      }
    }
    a=(betainit-beta1).adjoint()*(betainit-beta1);
    dif = sqrt(a);
    iter=iter+1;
  }
  List out;
  out["betainit"]=betainit;
  out["iter"]=iter;
  return out;
}


// [[Rcpp::export]]
List coordi_c(VectorXd atilde, MatrixXd B, VectorXd betainit, double lambda, int iter_max) {
  int iter=1, p=betainit.size();
  double dif=1, a, b, c;
  VectorXd beta1;
  while(iter<=iter_max && dif>1e-4){
    beta1=betainit;
    for(int i=0; i<p; i++){
      a=B(i,i);
      c=betainit.adjoint()*B.col(i);
      b=atilde(i) + c - a*betainit(i);
      betainit(i)=-soft_c(b/a,lambda/a);
    }
    a=(betainit-beta1).adjoint()*(betainit-beta1);
    dif = sqrt(a);
    iter=iter+1;
  }
  List out;
  out["betainit"]=betainit;
  out["iter"]=iter;
  return out;
}

