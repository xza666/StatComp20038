#include <Rcpp.h>
using namespace Rcpp;

//' @title A Metropolis sampler using Rcpp
//' @description A Metropolis sampler using Rcpp
//' @param sigma Variances  used for the proposal distribution
//' @param x0 Initial value
//' @param N the number of samples
//' @return a random sample of size \code{n} and the rejection rate
//' @examples
//' \dontrun{
//' rw=rwMC(2,25,1000)
//' plot(rw$x,type="l",lwd=2,xlab="sigma=2",ylab="X")
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix rwMC(double sigma,double x0,int N) {
  NumericMatrix mat(N, 2);
  mat(0,0)=x0;
  mat(0,1)=0;
  double y=0,u=0;
  for(int i = 2; i < N+1; i++) {
    y=rnorm(1,mat(i-2,0),sigma)[0];
    u=runif(1,0,1)[0];
    if(u<=exp(-abs(y))/exp(-abs(mat(i-2,0)))){
      mat(i-1,0)=y;
      mat(i-1,1)=mat(i-2,1);
    }
    else{
      mat(i-1,0)=mat(i-2,0);
      mat(i-1,1)=mat(i-2,1)+1;
    }
  }
  return(mat);
}
