#include <RcppArmadillo.h>
#include "utilities.h"

// Model-specific component samplers (i.e., full conditional moves that sample per each iteration)
// "Mixed variance" model -- not useful for predictions, but could be a useful intermediary step
// Less space complexity, cuz I omit mixed effects totally for now

using namespace arma;

// 1. Fixed (vectorized) coefficient matrix + fixed precision

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_betaSigma(vec& beta, double& sig2inv, vec& mnorms_v, mat& beta_m,
                        mat X, const mat& Y, const vec& tau2invs, const mat& theta,
                        const vec& beta_diags, double rho) {
    const uword N = X.n_rows;
    const uword P = X.n_cols;
    const uword Q = Y.n_cols;
    double sig_tr_total = 0.0; // Accumulates all the needed traces for sampling sig2inv
    mat centered = Y - theta;
    centered.each_col() %= 1.0 / arma::sqrt(tau2invs);
    const vec beta_umean = vectorise(X.t() * centered);
    vec inorms = rnorm_v(static_cast<int>(P*Q));
    // Added step: weight each column by patient-wise tau's
    X.each_col() /= arma::sqrt(tau2invs); // X was passed in by value
    // Iterate the diagonal Q blocks and accumulate Cholesky factors separately
    for (uword q=1; q < Q+1; ++q) {
        uword start = (q-1)*P;
        uword end = q*P-1;
        // Decompose R'R = Prior precision[inds] + (1-rho)*Z'DZ
        mat w_xtdx = (1.0-rho) * X.t() * X; // self-product of X, weighted by 1-rho
        w_xtdx.diag() += 1.0/beta_diags(span(start, end));
        mat chol_m = chol(w_xtdx);
        // Solve for R'x = vec Z'(Y-all the rest mean) and store
        vec quad_v = vectorise(forwardsub(chol_m.t(), beta_umean(span(start, end))));
        // Accumulate (1-rho)^2*||quads_v||^2
        sig_tr_total += std::pow(1.0-rho, 2) * accu(quad_v);
        // Solve for R(posterior mean) = (1-rho)*x 
        beta(span(start, end)) = (1.0-rho)*backsub(chol_m, quad_v);
        // Solve for R(MVN vector, not scaled by sigma^2) = x
        mnorms_v(span(start, end)) = vectorise(backsub(chol_m, inorms(span(start, end))));
    }
    // Sample sigma^-2 ~ Gamma(.5+.5*N*Q, .5+.5*((1-rho)*(Y-rest)'(Y-rest) - accumulated (1-rho)^2*||quads_v||^2))
    // Remember scale parameter takes the inverse of rate (usually used in R)
    sig2inv = R::rgamma(.5+static_cast<double>(.5*N*Q), 
                        1.0/(.5+.5*((1.0-rho)*trace(centered.t()*centered) - sig_tr_total)));
    // Set beta <- posterior mean + MVN vector / sigma^-2
    beta = beta + mnorms_v / sig2inv;
    beta_m = reshape(beta, P, Q);
}

// 2. TODO: Sample patient-specific precisions
// Assumes an N-length vector to fill in those precisions (J unique)

// 3. Sampling hierarchical parameters for beta (c2)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_c2(vec& c2, vec& beta_diags, 
                 const vec& beta, const double sig2inv, const uvec small_inds,
                 const double c0, const double beta_var0) {
    for (uword i = 0; i < small_inds.n_elem; ++i) {
        uword bInd = small_inds(i); 
        double Chi = std::pow(beta(bInd),2) * sig2inv; // Argument passed to GIG RNG
        c2(i) = rgigRcpp(0.0, Chi, std::pow(c0, -2));
        beta_diags(small_inds(i)) = c2(i) * beta_var0;
    }
}

// 4. Sampling theta (data augmentation for less expensive CAR sampling)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_theta(mat& theta, const mat& Y, const mat& X, const mat& beta_m, const vec& tau2invs, 
                    const double sig2inv, const double rho, const mat& r_y) {
    const uword N = Y.n_rows;
    const uword Q = Y.n_cols;
    const mat centered = Y - X * beta_m;
    mat post_mean_t = forbacksolve(r_y, centered.t()) * (1.0-rho);
    mat post_mnorms_t = backsub(r_y, rnorm_v(Q, N)) / sig2inv;
    // Added step: weight each column by patient-wise tau's
    post_mnorms_t.each_col() /= arma::sqrt(tau2invs);
    theta = post_mean_t.t() + post_mnorms_t.t();
    vec rowMeans = arma::mean(theta, 1); // For sum-to-zero constraint of each image effect
    theta.each_col() -= rowMeans;
}

// 5. Impute missing elements in target based on CAR full conditionals

// [[Rcpp::depends(RcppArmadillo)]]
void impute_car(mat& target, const mat& means, const vec& tau2invs, const double rho,
                const umat& mis_inds, const uvec& n_ns, const uvec& neighbors) {
    // No explicit NA handling, they must have been "filled out" from R as it is now
    // O/w erroneous results -- it'll be best to pre-process in R 
    uword counter = 0;
    uvec mis_rows = mis_inds.col(0);
    for (uword j = 0; j < mis_inds.n_rows; ++j) { // iterate across locations
        const uword mis_loc = mis_inds(j,1);
        const uvec curr_n = neighbors(span(counter,counter+n_ns(j)-1))-1;
        const double denom = n_ns(j)*rho+1-rho;
        for (const uword& i : mis_rows) { // iterate across patients
            const rowvec target_i = target.row(i);
            const rowvec curr_n_i = target_i.cols(curr_n); 
            const double mis_mean = means(i,mis_loc);
             target(i,mis_loc) = (rho*accu(curr_n_i) + (1.0-rho)*mis_mean) / denom + 
                // Added step: instead of a single precision, patient-specific precision at each row i
                 R::rnorm(0.0,1.0) / std::sqrt(tau2invs(i) * denom);
         }
         counter += n_ns(j);
     }
 }

// TODO: modify the Gibbs sampler
