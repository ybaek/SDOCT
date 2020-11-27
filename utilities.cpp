#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
vec runif_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands(i) = R::runif(0.0, 1.0);
    return rands;
}

// [[Rcpp::export]]
vec rnorm_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands(i) = R::rnorm(0.0, 1.0);
    return rands;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat backsub(const mat& r, const mat& x) {
    // solves for Y in RY = X with R upper triangular
    return solve(trimatu(r), x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat forwardsub(const mat& l, const mat& x) {
    // solves for Y in LY = X with L lower triangular
    return solve(trimatl(l), x);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat forbacksolve(const mat& r, const mat& x) {
    // solves for Y in (R'R)Y = X with R upper triangular
    mat ry = forwardsub(r.t(), x);
    return backsub(r, ry);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat altbacksolve(const mat& r1, const mat& r2, const mat& x) {
    // solves for Y in R_1YR_2' = X with R_1, R_2 upper triangular
    mat yr2t = backsub(r1, x);
    mat yt = backsub(r2, yr2t.t());
    return yt.t();
}
