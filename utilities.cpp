#include <RcppArmadillo.h>

using namespace arma;

// Vectorized pseudorandom & Gaussian RNGs, overloaded for vecs / mats

// [[Rcpp::depends(RcppArmadillo)]]
vec runif_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands[i] = R::runif(0.0, 1.0);
    return rands;
}

// [[Rcpp::depends(RcppArmadillo)]]
vec rnorm_v(const int N) {
    vec rands(N);
    for (int i=0; i < N; i++) 
        rands[i] = R::rnorm(0.0, 1.0);
    return rands;
}

// [[Rcpp::depends(RcppArmadillo)]]
mat runif_v(const int M, const int N) {
    mat rands(M, N);
    for (int i=0; i < M; i++) {
        for (int j=0; j < N; j++) {
            rands[i,j] = R::runif(0.0, 1.0);
        }
    }
    return rands;
}

// [[Rcpp::depends(RcppArmadillo)]]
mat rnorm_v(const int M, const int N) {
    mat rands(M, N);
    for (int i=0; i < M; i++) {
        for (int j=0; j < N; j++) {
            rands[i,j] = R::rnorm(0.0, 1.0);
        }
    }
    return rands;
}

// Normal log-density function, taking mat and returning mat
// (Does NOT take in vec/mat standard deviation)

// [[Rcpp::depends(RcppArmadillo)]]
mat ldnorm_v(const mat& x, const mat& mu_m, const double sigma) {
    const int M = x.n_rows;
    const int N = x.n_cols;
    mat densities(M, N);
    for (int i=0; i < M; i++) {
        for (int j=0; j < N; j++) {
            densities[i,j] = R::dnorm(x[i,j], mu_m[i,j], sigma, 1);
        }
    }
    return densities;
}

// Call Generalized inverse Gaussian RNG from R

double rgigRcpp(const double Lambda, const double Chi, const double Psi) {
    Rcpp::Environment GIGrvg = Rcpp::Environment::namespace_env("GIGrvg");
    Rcpp::Function rgig = GIGrvg["rgig"];
    SEXP rgigSEXP = rgig(1, Rcpp::Named("lambda")=Lambda, Rcpp::Named("chi")=Chi, Rcpp::Named("psi")=Psi);
    return Rcpp::as<double>(rgigSEXP);
}

// Linear algebra functions suitable for 
// multivariate Gaussian full conditional sampling

// [[Rcpp::depends(RcppArmadillo)]]
mat backsub(const mat& r, const mat& x) {
    // solves for Y in RY = X with R upper triangular
    return solve(trimatu(r), x);
}

// [[Rcpp::depends(RcppArmadillo)]]
mat forwardsub(const mat& l, const mat& x) {
    // solves for Y in LY = X with L lower triangular
    return solve(trimatl(l), x);
}

// [[Rcpp::depends(RcppArmadillo)]]
mat forbacksolve(const mat& r, const mat& x) {
    // solves for Y in (R'R)Y = X with R upper triangular
    mat ry = forwardsub(r.t(), x);
    return backsub(r, ry);
}

// [[Rcpp::depends(RcppArmadillo)]]
mat altbacksolve(const mat& r1, const mat& r2, const mat& x) {
    // solves for Y in R_1YR_2' = X with R_1, R_2 upper triangular
    mat yr2t = backsub(r1, x);
    mat yt = backsub(r2, yr2t.t());
    return yt.t();
}
