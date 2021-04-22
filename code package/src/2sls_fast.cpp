//based on https://github.com/AdaemmerP/lpirfs/blob/master/src/newey_west_tsls.cpp
#include <RcppArmadillo.h>
using namespace Rcpp;
//' @name newey_west_tsls
//' @title Compute 2SLS parameters and robust standard errors based on Newey-West
//' @description  Compute 2SLS parameters and robust standard errors based on Newey and West (1987).
//' Part of the function is based on the Matlab code by James P. LeSage.
//' @param y Numeric vector.
//' @param x Numeric matrix.
//' @param z Numeric matrix.
//' @param h Integer.
//' @return A vector, the 2SLS-Huber-White robust t-statistics for the IV coefficients.
//' @keywords internal
//' @references
//' Newey, W.K., and West, K.D. (1987). “A Simple, Positive-Definite, Heteroskedasticity and
//' Autocorrelation Consistent Covariance Matrix.” \emph{Econometrica}, 55, 703–708.
//' Wooldridge, J.M. (2002), Econometric Analysis of Cross Section and Panel Data, The MIT Press.
// [[Rcpp::export]]
arma::vec huber_white_tsls_tstats(arma::vec y, arma::mat x, arma::mat z){


  // 2SLS
  arma::mat xx_one   = arma::ones<arma::mat>(x.n_rows, 1); // Insert ones for constant
  x.insert_cols(0, xx_one);

  z.insert_cols(0, xx_one);

  // Build x_hat matrix with instrument matrix
  arma::mat  xx_hat   = z*inv(z.t()*z)*z.t()*x;

  // Estimate beta_iv
  arma::mat xpxi_iv  = inv(xx_hat.t()*xx_hat);
  arma::vec beta_iv  = xpxi_iv*xx_hat.t()*y;

  // Estimate corrected residuals
  arma::vec resids   = y - x*beta_iv;


  // Estimate normal cov-matrix of iv_estimators
  arma::vec resids_sq_iv  = resids%resids;


  // Use fitted data for the huber white estimator
  arma::mat bread     = inv(xx_hat.t()*xx_hat);
  arma::mat meat = arma::zeros(xx_hat.n_cols,xx_hat.n_cols);
  for(int i = 0; i<xx_hat.n_rows; ++i){
    arma::mat innerstes=xx_hat.row(i).t()*xx_hat.row(i);
    meat = meat + resids_sq_iv(i)*innerstes;
  }

  arma::mat V = bread*meat*bread;

  return (beta_iv/sqrt(V.diag()));
}
