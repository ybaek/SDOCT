#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
mat forbacksolve(const mat& r, const mat& x) {
    // solves for Y in (R'R)Y = X with R upper triangular
    return solve(trimatu(r), solve(trimatl(r.t()), x));
}

// [[Rcpp::depends(RcppArmadillo)]]
mat altbacksolve(const mat& r1, const mat& r2, const mat& x) {
    // solves for Y in R_1YR_2 = X with R_1, R_2 upper triangular
    return trans(solve(trimatu(r2), trans(solve(trimatu(r1), x))));
}

vec runif_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands(i) = R::runif(0.0, 1.0);
    return rands;
}

vec rnorm_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands(i) = R::rnorm(0.0, 1.0);
    return rands;
}

// [[Rcpp::depends(RcppArmadillo)]]
vec conjugateMNormal_blocked(const mat& X, const mat& Y, const mat& cov0inv, const double prec) {
    // TODO: check the dimensions
    int P = static_cast<int>(X.n_cols);
    mat xr = chol(cov0inv + X.t() * X);
    vec result = forbacksolve(xr, X.t() * Y) + solve(trimatu(xr), rnorm_v(P)) / prec;
    return result;
}

// [[Rcpp::depends(RcppArmadillo)]]
vec conjugateMNormal_unblocked(const mat& X, const mat& Y, const mat& cov0inv, const double prec) {
    // TODO: check the dimensions
    int P = static_cast<int>(X.n_cols);
    mat xr = chol(cov0inv + prec * X.t() * X);
    vec result = prec * forbacksolve(xr, X.t() * Y) + solve(trimatu(xr), rnorm_v(P));
    return result;
}
