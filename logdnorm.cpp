#include <Rcpp.h>
#include <cmath>        
#include <valarray>     

using namespace Rcpp;

// [[Rcpp::export]]
double logdnormCpp(NumericMatrix X, NumericMatrix mu, double sigma2) {
  int nrows = X.nrow();
  int ncolumns = X.ncol();
  double res = 0;
  double temp;
  for (int i = 0; i < nrows; i++){
    for (int j = 0; j < ncolumns; j++){
      temp = X(i, j)-mu(i, j);
      res -= 0.5*temp*temp/sigma2;
    }
  }
  res -= 0.5*nrows*ncolumns*log(2*PI*sigma2);
  return(res);
}

