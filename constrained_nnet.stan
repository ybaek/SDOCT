// FIXME: Code below as it is won't compile!
data {
    int<lower=1> N;               // no. of data points
    int<lower=1> J;               // no. of groups
    int<lower=1> K;               // no. of mixtures
    int<lower=1> Q;               // dim. of mv response 2
    int<lower=1> Nknots2;         // no. of knots 2
    int<lower=1> L;               // latent dimension
    int<lower=1,upper=J> jj[N];   // group vector
    vector[Nknots2-1] s[Nknots2]; // knot locations
    vector[Q] y[N];               // data (2D image)
    matrix[Q,Nknots2] K2;         // kernel matrix 2
}
parameters {
    simplex[K] pis;                 // mixture weights
    vector[L] xi[K];                // latent variables (shared across Z and Y)
    matrix[Nknots2,L] weights2[K];  // affine weights for Y
    vector[Nknots2] biases2[K];     // "biases" (intercepts) for Y
    real<lower=0> rho2;             // GP length scale for Y
    real<lower=0> tau2;             // GP scale for Y
    real<lower=0> alpha2;           // pre-tanh scale for Y
    real<lower=0> beta2;            // post-tanh scale for Y
    vector[Q] eta[J];               // auxiliary variables (for nested structure)
    vector[Q] intercept2[J];        // obs-level intercept for Y
    real<lower=0> sigma2;           // residual scale for Y
}
model {
    vector[K] log_pis = log(pis); // cache log
    vector[Q] fy[K]; // Latent Gaussian process mean function
    matrix[Q,Q] Cy;  // Latent covariance matrix
    {
        matrix[Nknots2,Nknots2] latent_cov2 = cov_exp_quad(s, tau2, rho2);
        vector[Nknots2] latent_mean2[K];
        Cy = quad_form_sym(latent_cov2, K2');
        // residual noise (otherwise this model is NOT proper!)
        for (q in 1:Q) {
            Cy[q,q] += square(sigma2);
        }
        for (k in 1:K) {
            // Nonlinearity thru scaled tanh 
            latent_mean2[k] = beta2 * tanh((biases2[k] + weights2[k] * xi[k]) / alpha2);
            fy[k] = intercept2 + K2 * latent_mean2[k];
        }
    }
    // Data conditional likelihood
    // Marginalize and increment mixture components
    for (n in 1:N) {
        vector[K] lps = log_pis;
        for (k in 1:K)
            lps[k] += multi_normal_lpdf(y[n] | intercept2[jj[n]] + fy[jj[n]], Cy);
        target += log_sum_exp(lps);
    }
    // Parameter priors
    rho2 ~ normal(0,1); 
    tau2 ~ normal(0,1);
    alpha2 ~ normal(0,.3);
    beta2 ~ normal(0,1);
    sigma2 ~ normal(0,1);
    for (k in 1:K) {
        xi[k] ~ std_normal();
        to_vector(weights2[k]) ~ std_normal();
    }
}
