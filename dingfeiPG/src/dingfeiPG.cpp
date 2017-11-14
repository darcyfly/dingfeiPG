#include <RcppArmadillo.h>
#include<Rcpp.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
inline arma::mat POmegaY(arma::mat& Y,arma::vec& index){
  for(int i = 0;i < index.n_rows;i++){
    Y(index(i)-1) = 0;
  }
  return Y;
}
inline arma::mat PCOmegaY(arma::mat& Y,arma::vec& Omega){
  for(int i = 0;i < Omega.n_rows;i++){
    Y(Omega(i)-1) = 0;
  }
  return Y;
}
inline arma::mat SVThresholding(arma::mat& X, double& lambda){
  mat U;
  vec s;
  mat V;
  svd(U,s,V,X);
  mat Xu = U.head_cols(s.n_elem);
  mat TV = V.t();
  mat Xv = TV.head_rows(s.n_elem);
  for(int i = 0;i < s.n_elem;i++){
    if((s(i)-lambda) > 0){
      s(i) = s(i)-lambda;
    }else{
      s(i) = 0;
    }
  }
  mat SVTlambda = Xu*diagmat(s)*Xv;
  return SVTlambda;
  }
//[[Rcpp::export]]
arma::mat dingfeiPG(arma::mat& X,  double& lambda,arma::vec& Omega, arma::vec& index, const int& maxinteration ){
  mat Y ;
  int n = X.n_cols;
  int m = X.n_rows;
  Y = zeros(size(X));
  mat MDold ;
  MDold = zeros<mat>(m,n);
  mat MDnew ;
  int i = 1;
  int told = 1;
  int tnew;
  mat Input;
  do {
    Input=PCOmegaY(Y,Omega)+X;
    MDnew=SVThresholding(Input,lambda);
    tnew = 1/2*(1+sqrt(1+4*told*told));
    Y = MDnew+((told-1)/tnew)*(MDnew-MDold);
    i = i +1;
  } while (i < maxinteration);
  for(i=0;i<Y.n_elem;i++){
    if(Y(i)> 1){
      Y(i) = 1;
    }else if(Y(i)<0){
      Y(i)=0;
    }
  }
  return Y;
}