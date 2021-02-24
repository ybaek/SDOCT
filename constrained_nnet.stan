data {
    int<lower=1> N;               // no. of data points
    int<lower=1> J;               // no. of groups
    int<lower=1> P;               // dim. of mv response 1
    int<lower=1> Q;               // dim. of mv response 2
    int<lower=1> Nknots1;         // no. of knots 1
    int<lower=1> Nknots2;         // no. of knots 2
    int<lower=1> L;               // latent dimension
    int<lower=1,upper=J> jj[N];   // group vector
    vector[Nknots1-1] x[Nknots1]; // knot locations
    vector[Nknots2-1] s[Nknots2]; // knot locations
    vector[P] z[N];               // data (1D image)
    vector[Q] y[N];               // data (2D image)
    matrix[P, Nknots1] K1;        // kernel matrix 1
    matrix[Q, Nknots2] K2;        // kernel matrix 2
}
parameters {
    vector[L] xi[J];                // latent variables (shared across Z and Y)
    matrix[Nknots1, L] weights1[J]; // affine weights for Z
    matrix[Nknots2, L] weights2[J]; // affine weights for Y
    vector[Nknots1] biases1[J];     // "biases" (intercepts) for Z
    vector[Nknots2] biases2[J];     // "biases" (intercepts) for Y
    real<lower=0> rho1;             // GP length scale for Z
    real<lower=0> rho2;             // GP length scale for Y
    real<lower=0> tau1;             // GP scale for Z
    real<lower=0> tau2;             // GP scale for Y
    real<lower=0> alpha1;           // pre-tanh scale for Z
    real<lower=0> alpha2;           // pre-tanh scale for Y
    real<lower=0> beta1;            // post-tanh scale for Z
    real<lower=0> beta2;            // post-tanh scale for Y
    matrix[P, J] intercept1;        // obs-level intercept for Z
    matrix[Q, J] intercept2;        // obs-level intercept for Y
    real<lower=0> sigma1;           // residual scale for Z
    real<lower=0> sigma2;           // residual scale for Y
}
model {
    // Latent Gaussian process mean function
    matrix[P, J] fz;
    matrix[Q, J] fy;
    matrix[P, P] Cz;
    matrix[Q, Q] Cy;
    {
        matrix[Nknots1, J] latent_mean1;
        matrix[Nknots2, J] latent_mean2;
        // TODO: rstan does not support other kernels yet
        matrix[Nknots1, Nknots1] latent_cov1 = cov_exp_quad(x, tau1, rho1);
        matrix[Nknots2, Nknots2] latent_cov2 = cov_exp_quad(s, tau2, rho2);
        Cz = quad_form_sym(latent_cov1, K1');
        Cy = quad_form_sym(latent_cov2, K2');
        for (j in 1:J) {
            // Nonlinearity thru scaled tanh 
            latent_mean1[,j] = beta1 * tanh((biases1[j] + weights1[j] * xi[j]) / alpha1);
            latent_mean2[,j] = beta2 * tanh((biases2[j] + weights2[j] * xi[j]) / alpha2);
        }
        fz = intercept1 + K1 * latent_mean1;
        fy = intercept2 + K2 * latent_mean2;
    }
    // Data conditional likelihood is independent with residual noise
    for (n in 1:N) {
        Cz[n,n] += square(sigma);
        Cy[n,n] += square(sigma);
        z[n] ~ multi_normal(fz[,jj[n]], Cz);
        y[n] ~ multi_normal(fy[,jj[n]], Cy);
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
    for (j in 1:J) {
        xi[j] ~ std_normal();
        to_vector(weights1[j]) ~ std_normal();
        to_vector(weights2[j]) ~ std_normal();
    } 
}
