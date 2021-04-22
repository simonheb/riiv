#include <RcppArmadillo.h>
#include <tictoc.h>
#include <ri_iv_test_cpp.h>
using namespace arma;



/*This is a simple function to efficiently compute a Mann–Whitney U test*/
// [[Rcpp::export]]
int ranktestU(arma::vec y, arma::vec t) {
  vec T1=y( find(t > 0));
  vec T2=y( find(t <= 0));
  int ret=0;
  for(uword i=0; i<T1.n_elem; i++) {
    for(uword j=0; j<T2.n_elem; j++) {
      ret += 0.5*(1+ (T1(i) < T2(j)) - (T2(j) < T1(i)));
    }
  }
  return (ret);
}

/*This function can return a number of different things, that may be useful as test statistics)*/
double test_statistic(double b, vec y, vec x, vec z, int mode=0) {
  vec adjusted_a = y-b*x;
  switch (mode) {
  case 0:
    return sum((z-mean(z)) % (adjusted_a - mean(adjusted_a))) / sum(square(z-mean(z)));
    break;
  case 1:
    return ranktestU(adjusted_a,z);
    break;
  case 2:
    mat rr=cov(z,adjusted_a);
    return rr(0);
    break;
  }
}

/*This function returns an RI-p-value, after imposing a null through b. You can specify which test statistic is used via the mode-varaible*/
// [[Rcpp::export]]
double ri_iv_test_p(double b, arma::vec  y, arma::vec x, arma::vec z, int r=1000, int seed=1234, int mode=0) {

  /*compute the value test statistic prior to resampling*/
  double sample_test_statistic = test_statistic(b,y,x,z,mode);

  /*set the  random seed*/
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));

  /*prepare the vector where the realizations throught the repetitions are store*/
  vec ri_test_statstics = zeros(r);

  vec z_permuted;
  for(uword i=0; i<r; i++) {
    /*draw permuted treatment*/
    z_permuted=randi<vec>( z.n_elem, distr_param(0,1));
    /*compute test statics*/
    ri_test_statstics(i) = test_statistic(b,y,x,z_permuted,mode);
  }

  /*compute p-value*/
  vec abshow_often_was_the_statistic_larger = zeros(size(ri_test_statstics));
  abshow_often_was_the_statistic_larger.elem( find(abs(ri_test_statstics) > fabs(sample_test_statistic)) ).ones();

  /*return p-value*/
  return(mean(abshow_often_was_the_statistic_larger));
}

