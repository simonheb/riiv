#include <RcppArmadillo.h>
#include <Rcpp.h>
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
      //std::cout << ret << std::endl;
      ret += 0.5*(1+ (T1(i) < T2(j)) - (T2(j) < T1(i)));
    }
  }
  return (ret);
}


// [[Rcpp::export]]
double ri_iv_test_p(double b, arma::vec  y, arma::vec x, arma::vec z, int r=1000, int seed=12345, bool ranktest=false) {
  vec adjusted_a = y-b*x;
  double sampleT;
  if (ranktest) {
    sampleT = ranktestU(adjusted_a,z);
  } else {
    sampleT = sum((z-mean(z)) % (adjusted_a - mean(adjusted_a))) / sum(square(z-mean(z)));
  }

  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));

  vec diffdist = zeros(r);

  vec z_permuted;
  for(uword i=0; i<r; i++) {
    z_permuted=randi<vec>( z.n_elem, distr_param(0,1)  );

    if (ranktest) {
      diffdist(i) = ranktestU(adjusted_a,z_permuted);
    } else {
      vec z_permuted_demeaned=z_permuted - mean(z_permuted);
      diffdist(i) = sum(z_permuted_demeaned % (adjusted_a - mean(adjusted_a))) / sum(square(z_permuted_demeaned));
    }
  }
  vec larger = zeros(size(diffdist));
  larger.elem( find(diffdist > sampleT) ).ones();
  return(1-2*(fabs(mean(larger)-0.5)));
}
