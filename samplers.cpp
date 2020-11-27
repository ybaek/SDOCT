#include <RcppArmadillo.h>

using namespace arma;

// Model-specific component samplers (i.e., full conditional moves that sample per each iteration)

// 1. Fixed (vectorized) coefficient matrix + observation precision

// [[Rcpp::depends(RcppArmadillo)]]
void nextmove_betaSigma(vec& beta, double& sig2inv, mat& chols_m, vec& mnorms_v,
                        const mat& X, const mat& Y, const mat& means, const vec& beta_diags, double rho) {
    const uword N = X.n_rows;
    const uword P = X.n_cols;
    const uword Q = Y.n_cols;
    double sig_tr_total = 0; // Accumulates all the needed traces for sampling sig2inv
    mat centered = Y - means;
    vec beta_umean = vectorise(X.t() * centered);
    // TODO include header
    vec inorms = rnorm_v(P*Q);
    // Iterate the diagonal Q blocks and accumulate Cholesky factors separately
    for (int q=0; q < Q; q++) {
        mat chol_m(P, Q);
        vec quad_v(P);
        uword start = (q-1)*P;
        uword end = q*P-1;
        // Decompose R'R = Prior precision[inds] + (1-rho)*Z'Z
        mat w_xtx = (1.0-rho) * X.t() * X; // self-product of X, weighted by 1-rho
        w_xtx.diag() += 1.0/beta_diags.subvec(start, end);
        chols_m = chol(w_xtx);
        // Solve for R'x = vec Z'(Y-all the rest mean) and store
        quad_v = forwardsub(chols_m.t(), beta_umean.subvec(start, end));
        // Accumulate (1-rho)^2*||quads_v||^2
        sig_tr_total += std::pow(1.0-rho, 2) * arma::accu(quad_v);
        // Solve for R(posterior mean) = (1-rho)*x 
        beta.subvec(start, end) = (1.0-rho)*backsub(chol_m, quad_v);
        // Solve for R(MVN vector, not scaled by sigma^2) = x
        mnorms_v.subvec(start, end) = backsub(chol_m, inorms.subvec(start, end));
    }
    // Sample sigma^-2 ~ Gamma(.5+.5*N*Q, .5+.5*((1-rho)*(Y-rest)'(Y-rest) - accumulated (1-rho)^2*||quads_v||^2))
    sig2inv = R::rgamma(.5+static_cast<double>(.5*N*Q), (1.0-rho)*trace(centered.t() * centered) - sig_tr_total);
    // Set beta <- posterior mean + MVN vector / sigma^-2
    beta = beta + mnorms_v / sig2inv;
}

// 2. Random coefficient matrices + hierarchical (nugget/spatial) precision

// 3. Sampling hierarchical parameters for beta (c2)

// 4. Sampling theta (data augmentation for less expensive CAR sampling)

// 5. Impute missing elements in target based on CAR full conditionals

// [[Rcpp::depends(RcppArmadillo)]]
void impute_car(mat& target, const mat& means, const double prec, const double rho,
                const umat& mis_inds, const uvec& n_ns, const uvec& neighbors) {
    double counter = 0;
    for (uword j = 0; j < mis_inds.n_rows; j++) {
        const uword mis_loc = mis_inds(j,1);
        const uvec curr_n = neighbors.subvec(counter,counter+n_ns(j)-1);
        const double denom = n_ns(j)*rho+1-rho;
        for (const uword& i : mis_inds.col(0)) {
            const rowvec target_i = target.row(i);
            const rowvec curr_n_i = target_i.elem(neighbors);
            const double mis_mean = means(i,mis_loc);
            target(i,mis_loc) = (rho*mis_mean + (1-rho)*arma::accu(curr_n_i)) / denom + 
                R::rnorm(0.0,1.0) / std::sqrt(prec * denom);
        }
        counter += n_ns(j);
    }
}
