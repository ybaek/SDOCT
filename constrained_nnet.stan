data {
    int<lower=1> N;
    // int<lower=1> J;
    int<lower=1> P;
    int<lower=1> Q;
    // int<lower=1> K;
    int<lower=1> Nknots;
    // int<lower=1,upper=J> jj[N];
    int<lower=1,upper=Q> kk[Nknots];
    vector[2] coords[Q];
    vector[Nknots-1] s_tf[Nknots];
    vector[P] z[N];
    vector[Q] y[N];
}
parameters {
    real<lower=0> rho;
    real<lower=0> tau;
    real<lower=0> lambda;
    real<lower=0> sigma;
    real<lower=0> alpha;
    real<lower=0> beta;
    vector[Q] intercept;
    vector[Nknots] b;
    matrix[Nknots,P] W;
}
transformed parameters {
    real rho_tf = sqrt(rho / 2);
    matrix[Q,Q] kernel = cov_exp_quad(coords, 1.0, lambda);
    matrix[Nknots,Nknots] latent_cov = cov_exp_quad(s_tf, tau, rho_tf);
}
model {
    {
        vector[Nknots] latent_mean[N];
        vector[Q] obs_mean[N];
        matrix[Q,Nknots] kernel_part = kernel[,kk];
        matrix[Q,Q] obs_cov = quad_form_sym(latent_cov, kernel_part'); 
        for (q in 1:Q)
            obs_cov[q,q] += square(sigma);
        for (n in 1:N) {
            // latent_mean[n] = b0 + rbiases[jj[n]] + W0 * z[n] + rweights[jj[n]] * z[n];
            latent_mean[n] = b + W * z[n];
            obs_mean[n] = beta * kernel_part * tanh(latent_mean[n] / alpha);
            y[n] ~ multi_normal(obs_mean[n], obs_cov);
        }
    }
    // Priors
    b ~ normal(0, 5.);
    to_vector(W) ~ double_exponential(0, 5.);
    alpha ~ normal(0, 5.);
    beta ~ normal(0, 10.);
    tau ~ normal(0, 10.);
    rho ~ inv_gamma(5., 20.);
    lambda ~ inv_gamma(5., 5.);
    sigma ~ normal(0, 10.);
}
