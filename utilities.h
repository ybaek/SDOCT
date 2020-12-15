#ifndef __utilities__
#define __utilities__

#include <RcppArmadillo.h>

arma::vec runif_v(const int N);
arma::vec rnorm_v(const int N);
arma::mat runif_v(const int M, const int N);
arma::mat rnorm_v(const int M, const int N);
arma::mat ldnorm_v(const arma::mat& x, const arma::mat& mu_m, const double sigma);
arma::mat ldnorm_v(const arma::mat& x, const arma::mat& mu_m, const arma::vec& sigma);
double rgigRcpp(const double Lambda, const double Chi, const double Psi);
arma::mat rWishart(const double df, const arma::mat r_inverse);
arma::mat backsub(const arma::mat& r, const arma::mat& x);
arma::mat forwardsub(const arma::mat& l, const arma::mat& x);
arma::mat forbacksolve(const arma::mat& r, const arma::mat& x);
arma::mat altbacksolve(const arma::mat& r1, const arma::mat& r2, const arma::mat& x);

#endif // __utilities__