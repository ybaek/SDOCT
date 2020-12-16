#include <RcppArmadillo.h>
#include "utilities.h"

// Model-specific component samplers (i.e., full conditional moves that sample per each iteration)

using namespace arma;

// 1. Fixed (vectorized) coefficient matrix + observation precision

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_betaSigma(vec& beta, double& sig2inv, vec& mnorms_v, mat& beta_m,
                        const mat& X, const mat& Y, const mat& gamma_m, const mat& theta,
                        const vec& beta_diags, double rho) {
    const uword N = X.n_rows;
    const uword P = X.n_cols;
    const uword Q = Y.n_cols;
    double sig_tr_total = 0.0; // Accumulates all the needed traces for sampling sig2inv
    const mat centered = Y - gamma_m - theta;
    const vec beta_umean = vectorise(X.t() * centered);
    vec inorms = rnorm_v(static_cast<int>(P*Q));
    // Iterate the diagonal Q blocks and accumulate Cholesky factors separately
    for (uword q=1; q < Q+1; ++q) {
        uword start = (q-1)*P;
        uword end = q*P-1;
        // Decompose R'R = Prior precision[inds] + (1-rho)*Z'Z
        mat w_xtx = (1.0-rho) * X.t() * X; // self-product of X, weighted by 1-rho
        w_xtx.diag() += 1.0/beta_diags(span(start, end));
        mat chol_m = chol(w_xtx);
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

// 2. Random coefficient matrices + hierarchical precision matrix
// Assumes a random effects matrix derived by Z(j,) * gamma(,,j) 
// Precision matrix has a Wishart prior

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_gamma(cube& gamma, mat& gamma_m, const mat& psi,
                    const mat& X, const mat& Y, const mat& beta_m, const double sig2inv,
                    const uvec& ids, const mat& r_y) {
    const uword P = X.n_cols;
    const uword Q = Y.n_cols;
    const uword J = ids.max();
    for (uword j = 1; j < J+1; ++j) { // Group label ints are from R and thus start from 1
        uvec jinds = find(ids==j);
        const mat X_slice = X.rows(jinds);
        const mat Y_slice = Y.rows(jinds);
        const mat centered = Y_slice - X_slice * beta_m;
        // Decompose R'R = tau^-2 * precision_row + Z[j,]'Z[j,] (row-wise covariance)
        mat post_prec_chol = chol(psi + X_slice.t() * X_slice);
        // Solve for R'R(posterior mean) = vec Z'(Y-rest)
        // Also solve for (R(MN matrix)R_y')' = (PxQ normal matrix) where R_y is static (col-wise covariance Cholesky)
        // and set gamma[,,j] <- posterior mean + MN matrix / (sigma^-2 * tau^-2)
        gamma.slice(j-1) = altbacksolve(post_prec_chol, r_y, rnorm_v(P, Q)) / sig2inv + 
            forbacksolve(post_prec_chol, X_slice.t() * centered);
        gamma_m.rows(jinds) = X_slice * gamma.slice(j-1);
    }
}

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_psi(mat& psi, const cube& gamma, 
                  const double sig2inv, const mat& cov_x, const mat& r_y) {
    const uword P = gamma.n_rows;
    const uword Q = gamma.n_cols;
    const uword K = gamma.n_slices;
    double post_df = static_cast<double>(1 + P  + .5 * (K * Q));
    mat post_prec(P, P, fill::zeros);
    for (uword k = 0; k < K; ++k) {
        post_prec += gamma.slice(k) * r_y * r_y.t() * gamma.slice(k).t();
    }
    post_prec *= sig2inv;
    mat post_r = chol(post_prec + cov_x); // Inverse scale is the (sparse) covariance
    psi = rWishart(post_df, post_r);
}

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
void nextmove_theta(mat& theta, const mat& Y, const mat& X, const mat& beta_m, const mat& gamma_m, 
                    const double sig2inv, const double rho, const mat& r_y) {
    const uword N = Y.n_rows;
    const uword Q = Y.n_cols;
    const mat centered = Y - X * beta_m - gamma_m;
    mat post_mean_t = forbacksolve(r_y, centered.t()) * (1.0-rho);
    mat post_mnorms_t = backsub(r_y, rnorm_v(Q, N)) / sig2inv;
    theta = post_mean_t.t() + post_mnorms_t.t();
    vec rowMeans = arma::mean(theta, 1); // For sum-to-zero constraint of each image effect
    theta.each_col() -= rowMeans;
}

// 5. Impute missing elements in target based on CAR full conditionals

// [[Rcpp::depends(RcppArmadillo)]]
void impute_car(mat& target, const mat& means, const double prec, const double rho,
                const umat& mis_inds, const uvec& n_ns, const uvec& neighbors) {
    // No explicit NA handling, they must have been "filled out" from R as it is now
    // O/w erroneous results -- it'll be best to pre-process in R 
    uword counter = 0;
    uvec mis_rows = mis_inds.col(0);
    uvec mis_cols = mis_inds.col(1);
    for (uword j = 0; j < mis_inds.n_rows; ++j) { // iterate across locations
        const uword mis_unit = mis_rows(j);
        const uword mis_loc = mis_cols(j);
        const uvec curr_n = neighbors(span(counter,counter+n_ns(j)-1))-1;
        const double denom = n_ns(j)*rho + 1.0 - rho;
        const rowvec target_i = target.row(mis_unit);
        const rowvec curr_n_i = target_i.cols(curr_n); 
        const double mis_mean = means(mis_unit, mis_loc);
        target(mis_unit, mis_loc) = (rho*accu(curr_n_i) + (1.0-rho)*mis_mean) / denom + 
            R::rnorm(0.0, 1.0) / std::sqrt(prec * denom);
        counter += n_ns(j);
     }
 }

// Main routine for matrix LMM 

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
mat mainSampler_lmm(const Rcpp::List& data, const Rcpp::List& inits, const Rcpp::List& hyper, const Rcpp::List& mcmc) {
    // Constants
    const mat X = Rcpp::as<mat>(data["X"]);
    const uvec ids = Rcpp::as<uvec>(data["ids"]);
    const uvec n_ns = Rcpp::as<uvec>(data["n_ns"]);
    const uvec neighbors = Rcpp::as<uvec>(data["neighbors"]);
    const umat mis_inds = Rcpp::as<umat>(data["mis_inds"])-1; // Indices are subtracted 1 from R
    const uvec small_inds = Rcpp::as<uvec>(data["small_inds"])-1; // Indices are subtracted 1 from R
    
    const mat prec_y = Rcpp::as<mat>(hyper["prec_y"]);
    const mat r_y = chol(prec_y); // Pass in a decomposed matrix
    const mat cov_x = Rcpp::as<mat>(hyper["cov_x"]);
    const double rho = hyper["rho"];
    const double c0 = hyper["c0"];
    const double beta_var0 = hyper["beta_var0"];
    
    const int I = mcmc["I"];
    const int burnin = mcmc["burnin"];

    // Pass-in-references
    mat Y = Rcpp::as<mat>(data["Y"]);
    vec beta = Rcpp::as<vec>(inits["beta"]);
    vec c2 = Rcpp::as<vec>(inits["c2"]);
    cube gamma = Rcpp::as<cube>(inits["gamma"]);
    mat psi = Rcpp::as<mat>(inits["psi"]);
    mat theta = Rcpp::as<mat>(inits["theta"]);
    double sig2inv = inits["sig2inv"];

    // Auxiliary containers for sampling
    vec beta_diags(beta.n_rows);
    beta_diags.fill(beta_var0);
    beta_diags.rows(small_inds) %= c2;
    vec mnorms_v(beta.n_rows, fill::zeros);
    mat beta_m(X.n_cols, Y.n_cols, fill::zeros);
    mat gamma_m(Y.n_rows, Y.n_cols, fill::zeros);
    vec lpd(Y.n_rows, fill::zeros); // Log point-wise likelihoods
    
    // This line can take different values 
    // depending on what parameters are stored
    // TODO: give user the control!
    int par1_length = beta.n_rows;
    int par2_length = 1;
    int par3_length = small_inds.n_rows; 
    // Store the lower diagonal of psi, including main diagonal
    uvec lower_inds = trimatl_ind(size(psi));
    int par4_length = lower_inds.n_rows;
    int par5_length = Y.n_rows;
    int out_length = par1_length + par2_length + par3_length + par4_length + par5_length;
    mat out = mat(out_length, I-burnin, fill::zeros);

    for (uword iter = 0; iter < static_cast<uword>(I); ++iter) {
        nextmove_betaSigma(beta, sig2inv, mnorms_v, beta_m, X, Y, gamma_m, theta, beta_diags, rho);
        nextmove_gamma(gamma, gamma_m, psi, X, Y, beta_m, sig2inv, ids, r_y);
        nextmove_psi(psi, gamma, sig2inv, cov_x, r_y);
        nextmove_c2(c2, beta_diags, beta, sig2inv, small_inds, c0, beta_var0);
        nextmove_theta(theta, Y, X, beta_m, gamma_m, sig2inv, rho, r_y);
        impute_car(Y, X * beta_m + gamma_m + theta, sig2inv, rho, mis_inds, n_ns, neighbors);
        lpd = arma::sum(ldnorm_v(Y, X * beta_m + gamma_m + theta, std::pow(sig2inv, -.5)), 1); // Update log pointwise densities
        
        if (iter >= static_cast<uword>(burnin)) {
            // FIXME: unsigned int comparison. burnin=0 leads to error
            uword keep_iter = iter - burnin;
            out(span(0, par1_length-1), keep_iter) = beta;
            out(par1_length, keep_iter) = sig2inv;
            out(span(par1_length+par2_length, par1_length+par2_length+par3_length-1), keep_iter) = c2;
            out(span(par1_length+par2_length+par3_length, out_length-par5_length-1), keep_iter) = psi(lower_inds);
            out(span(out_length-par5_length, out_length-1), keep_iter) = lpd;
        }
    }
    return out;
}
