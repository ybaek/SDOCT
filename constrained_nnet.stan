data {
    int<lower=1> N;             // no. of data points
    int<lower=1> J;             // no. of groups
    int<lower=1> P;             // dim. of mv response 1
    int<lower=1> Q;             // dim. of mv response 2
    int<lower=1> Nknots1;      // no. of knots 1
    int<lower=1> Nknots2;      // no. of knots 2
    int<lower=1> L;             // latent dimension
    int<lower=1,upper=J> jj[N]; // group vector
    real x[Nknots1];             // knot locations (scalar)
    vector[2] s[Nknots2];        // knot locations (R2 plane)
    matrix[N, P] z;             // data (1D image)
    matrix[N, Q] y;             // data (2D image)
    matrix[Nknots1, P] Kt1;      // kernel matrix 1(transposed format)
    matrix[Nknots2, Q] Kt2;      // kernel matrix 2(transposed format)
}
parameters {
    vector[L] xi[J];               // latent variables (shared across Z and Y)
    matrix[Nknots1, L] weights1[J]; // affine weights for Z
    matrix[Nknots2, L] weights2[J]; // affine weights for Y
    real<lower=0> rho1;            // GP length scale for Z
    real<lower=0> rho2;            // GP length scale for Y
    real<lower=0> tau1;            // GP scale for Z
    real<lower=0> tau2;            // GP scale for Y
    real<lower=0> alpha1;          // pre-tanh scale for Z
    real<lower=0> alpha2;          // pre-tanh scale for Y
    real<lower=0> beta1;           // post-tanh scale for Z
    real<lower=0> beta2;           // post-tanh scale for Y
    real<lower=0> sigma1;          // residual scale for Z
    real<lower=0> sigma2;          // residual scale for Y
}
model {
    /* Since Stan cannot efficiently handle matrix normals yet,
       everything needs to be implemented in vectorized forms
       latent variable eta is introduced in particular so that
       we don't end up having a giant kronecker covariance */

    // Latent Gaussian process mean function
    matrix[J, P] fz;
    matrix[J, Q] fy;
    // Auxiliary variables for GP covariance
    matrix[Nknots1, J] eta_z;
    matrix[Nknots2, J] eta_y;
    {
        matrix[J, Nknots1] latent_mean1;
        matrix[J, Nknots2] latent_mean2;
        matrix[Nknots1, Nknots1] latent_cov1 = cov_exp_quad(x, tau1, rho1);
        matrix[Nknots1, Nknots1] L_lc1 = cholesky_decompose(latent_cov1);
        matrix[Nknots2, Nknots2] latent_cov2 = cov_exp_quad(s, tau2, rho2);
        matrix[Nknots2, Nknots2] L_lc2 = cholesky_decompose(latent_cov2);
        for (j in 1:J) {
            latent_mean1[j, ] = beta1 * tanh(weights1[j] * xi[j] / alpha1)';
            latent_mean2[j, ] = beta2 * tanh(weights2[j] * xi[j] / alpha2)';
        }
        fz = latent_mean1 + (L_lc1 * eta_z)' * Kt1;
        fy = latent_mean2 + (L_lc2 * eta_y)' * Kt2;
    }
    // Data conditional likelihood is independent
    for (n in 1:N) {
        to_vector(z[n,]) ~ normal(to_vector(fz[jj[n],]), sigma1);
        to_vector(y[n,]) ~ normal(to_vector(fy[jj[n],]), sigma2);
    }

    // Parameter priors
    rho1 ~ inv_gamma(5, 5); 
    rho2 ~ inv_gamma(5, 5); 
    tau1 ~ normal(0, 10);
    tau2 ~ normal(0, 10);
    sigma1 ~ normal(0, 10);
    sigma2 ~ normal(0, 10);
    alpha1 ~ normal(0, .3);
    alpha2 ~ normal(0, .3);
    beta1 ~ normal(0, 10);
    beta2 ~ normal(0, 10);
    to_vector(eta_z) ~ std_normal();
    to_vector(eta_y) ~ std_normal();
    for (j in 1:J) {
        xi[j] ~ std_normal();
        to_vector(weights1[j]) ~ std_normal();
        to_vector(weights2[j]) ~ std_normal();
    } 
}
generated quantities {
}
