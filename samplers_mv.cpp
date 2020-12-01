#include <RcppArmadillo.h>
#include "utilities.h"

// Model-specific component samplers (i.e., full conditional moves that sample per each iteration)
// "Mixed variance" model -- not useful for predictions, but could be a useful intermediary step
// Less space complexity, cuz I omit mixed effects totally for now

using namespace arma;

// 1. Fixed (vectorized) coefficient matrix (independent prior spec.)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_beta(vec& beta, mat& beta_m,
                    const mat& X, const mat& Y, const vec& tau2invs, const mat& theta,
                    const vec& beta_diags, double rho) {
    // Should I just pass X by value?
    const uword N = X.n_rows;
    const uword P = X.n_cols;
    const uword Q = Y.n_cols;
    // Added step: weight each column by patient-wise tau's
    mat Xd = X; Xd.each_col() %= arma::sqrt(tau2invs);
    mat centered = Y - theta; centered.each_col() %= arma::sqrt(tau2invs);
    const vec beta_umean = vectorise(Xd.t() * centered);
    vec inorms = rnorm_v(static_cast<int>(P*Q));
    // Iterate the diagonal Q blocks and accumulate Cholesky factors separately
    for (uword q=1; q < Q+1; ++q) {
        uword start = (q-1)*P;
        uword end = q*P-1;
        // Decompose R'R = Prior precision[inds] + (1-rho)*Z'DZ
        mat w_xtdx = (1.0-rho) * Xd.t() * Xd; // self-product of X, weighted by 1-rho
        w_xtdx.diag() += 1.0/beta_diags(span(start, end));
        mat chol_m = chol(w_xtdx);
        // Solve for R'x = vec Z'(Y-all the rest mean) and store
        vec quad_v = vectorise(forwardsub(chol_m.t(), beta_umean(span(start, end))));
        // Solve for R(posterior mean) = (1-rho)*x 
        beta(span(start, end)) = (1.0-rho)*backsub(chol_m, quad_v) + 
            vectorise(backsub(chol_m, inorms(span(start, end))));
    }
    beta_m = reshape(beta, P, Q);
}

// 2. Sample patient-specific precisions
// Assumes an N-length vector to fill in those precisions (J unique)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_tau2invs(vec& tau2invs,
                       const mat& X, const mat& Y, const mat& beta_m, 
                       const uvec& ids, const mat& r_y) {
    const uvec ids_u = arma::unique(ids);
    for (uword j : ids_u) { // Group label ints are from R and thus start from 1
        uvec jinds = find(ids==j);
        const mat X_slice = X.rows(jinds);
        const mat Y_slice = Y.rows(jinds);
        const mat centered = Y_slice - X_slice * beta_m;
        double sse = trace(r_y.t() * r_y * centered.t() * centered);
        // Note Rgamma takes scale, not inverse scale (as is frequently done in R)
        double tau2inv_j = R::rgamma(.5+static_cast<double>(.5*Y_slice.n_elem), 1.0 / (.5 + .5*sse));
        tau2invs.rows(jinds).fill(tau2inv_j); 
    }
}

// 3. Sampling hierarchical parameters for beta (c2)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_c2(vec& c2, vec& beta_diags, 
                 const vec& beta, const uvec small_inds,
                 const double c0, const double beta_var0) {
    for (uword i = 0; i < small_inds.n_elem; ++i) {
        uword bInd = small_inds(i); 
        double Chi = std::pow(beta(bInd),2); // Argument passed to GIG RNG
        c2(i) = rgigRcpp(0.0, Chi, std::pow(c0, -2));
        beta_diags(small_inds(i)) = c2(i) * beta_var0;
    }
}

// 4. Sampling theta (data augmentation for less expensive CAR sampling)

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_theta(mat& theta, const mat& Y, const mat& X, const mat& beta_m, const vec& tau2invs, 
                    const double rho, const mat& r_y) {
    const uword N = Y.n_rows;
    const uword Q = Y.n_cols;
    const mat centered = Y - X * beta_m;
    mat post_mean_t = forbacksolve(r_y, centered.t()) * (1.0-rho);
    mat post_mnorms_t = backsub(r_y, rnorm_v(Q, N));
    // Added step: weight each column by patient-wise tau's (Here VARIANCE not precisions!)
    post_mnorms_t.each_row() /= arma::sqrt(tau2invs).as_row();
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

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat mainSampler_mv(const Rcpp::List& data, const Rcpp::List& inits, const Rcpp::List& hyper, const Rcpp::List& mcmc) {
    // Constants
    const mat X = Rcpp::as<mat>(data["X"]);
    const uvec ids = Rcpp::as<uvec>(data["ids"]);
    const uvec n_ns = Rcpp::as<uvec>(data["n_ns"]);
    const uvec neighbors = Rcpp::as<uvec>(data["neighbors"]);
    const umat mis_inds = Rcpp::as<umat>(data["mis_inds"])-1; // Indices are subtracted 1 from R
    const uvec small_inds = Rcpp::as<uvec>(data["small_inds"])-1; // Indices are subtracted 1 from R
    
    const mat prec_y = Rcpp::as<mat>(hyper["prec_y"]);
    const mat r_y = chol(prec_y); // Pass in a decomposed matrix
    // const mat prec_x = Rcpp::as<mat>(hyper["prec_x"]);
    const double rho = hyper["rho"];
    const double c0 = hyper["c0"];
    const double beta_var0 = hyper["beta_var0"];
    
    const int I = mcmc["I"];
    const int burnin = mcmc["burnin"];

    // Pass-in-references
    mat Y = Rcpp::as<mat>(data["Y"]);
    vec beta = Rcpp::as<vec>(inits["beta"]);
    vec c2 = Rcpp::as<vec>(inits["c2"]);
    // cube gamma = Rcpp::as<cube>(inits["gamma"]);
    mat theta = Rcpp::as<mat>(inits["theta"]);
    vec tau2invs = Rcpp::as<vec>(inits["tau2invs"]);
    // double sig2inv = inits["sig2inv"];
    // double tau2inv = inits["tau2inv"];

    // Auxiliary containers for sampling
    vec beta_diags(beta.n_rows);
    beta_diags.fill(beta_var0);
    beta_diags.rows(small_inds) %= c2;
    mat beta_m(X.n_cols, Y.n_cols, fill::zeros);
    //vec mnorms_v(beta.n_rows, fill::zeros);
    // mat gamma_m(Y.n_rows, Y.n_cols, fill::zeros);
    vec lpd(Y.n_rows, fill::zeros); // Log point-wise likelihoods
    
    mat out = mat(beta.n_rows+small_inds.n_rows+Y.n_rows+Y.n_rows, I-burnin, fill::zeros);
    for (uword iter = 0; iter < static_cast<uword>(I); ++iter) {
        nextmove_beta(beta, beta_m, X, Y, tau2invs, theta, beta_diags, rho);
        nextmove_tau2invs(tau2invs, X, Y, beta_m, ids, r_y);
        nextmove_c2(c2, beta_diags, beta, small_inds, c0, beta_var0);
        nextmove_theta(theta, Y, X, beta_m, tau2invs, rho, r_y);
        impute_car(Y, X * beta_m + theta, tau2invs, rho, mis_inds, n_ns, neighbors);
        lpd = arma::sum(ldnorm_v(Y, X * beta_m + theta, arma::pow(tau2invs, -.5)), 1); // Update log pointwise densities
        
        if (iter >= static_cast<uword>(burnin)) {
            uword keep_iter = iter - burnin;
            out(span(0, beta.n_rows-1), keep_iter) = beta;
            out(span(beta.n_rows, beta.n_rows+small_inds.n_rows-1), keep_iter) = c2;
            out(span(beta.n_rows+small_inds.n_rows, beta.n_rows+small_inds.n_rows+Y.n_rows-1), keep_iter) = tau2invs;
            out(span(beta.n_rows+small_inds.n_rows+Y.n_rows, out.n_rows-1), keep_iter) = lpd;
        }
    }
    return out;
}
