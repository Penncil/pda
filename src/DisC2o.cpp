#include <Rcpp.h>
using namespace Rcpp;  
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
  } else {
    out=0;
  } 
  return out;
} 


// [[Rcpp::export]]
List coordi_cho(NumericVector atilde, NumericMatrix B, NumericVector betainit, double lambda, int iter_max) {
  int iter=1, p=betainit.size();
  double dif=1, a, b, c;
  NumericVector beta1;
  while(iter<=iter_max && dif>1e-4){
    beta1=betainit;
    for(int i=0; i<p; i++){
      a=B(i,i);
      c= Rcpp::sum(betainit * B(_, i));
      b=atilde(i) + c - a*betainit(i);
      if(i==0){
        betainit(i)=-b/a;
      } else{
        betainit(i)=-soft_c(b/a,lambda/a);
      }
    }
    a= Rcpp::sum((betainit-beta1) * (betainit-beta1));
    dif = sqrt(a);
    iter=iter+1;
  } 
  List out;
  out["betainit"]=betainit;
  out["iter"]=iter;
  return out;
} 


// [[Rcpp::export]]
List coordi_c(NumericVector atilde, NumericMatrix B, NumericVector betainit, double lambda, int iter_max) {
  int iter=1, p=betainit.size();
  double dif=1, a, b, c;
  NumericVector beta1;
  while(iter<=iter_max && dif>1e-4){
    beta1=betainit;
    for(int i=0; i<p; i++){
      a=B(i,i);
      c= Rcpp::sum(betainit * B(_, i));  
      b=atilde(i) + c - a*betainit(i);
      betainit(i)=-soft_c(b/a,lambda/a);
    } 
    a= Rcpp::sum((betainit-beta1) * (betainit-beta1));
    dif = sqrt(a);
    iter=iter+1;
  }
  List out;
  out["betainit"]=betainit;
  out["iter"]=iter;
  return out;
} 

