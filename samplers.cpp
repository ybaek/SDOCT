#include <RcppArmadillo.h>

using namespace arma;

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