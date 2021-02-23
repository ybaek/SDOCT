data {
    int<lower=1> K; // No. of clusters
    int<lower=1> N; // No. of generated samples
}
generated quantities {  
    real v[K];
    real sigma[K];
    real mu[K];
    vector[K] weights;
    real y[N];
    int z[N];
    
    for (k in 1:(K-1)) {
        v[k]     = beta_rng(1., 10.);
        sigma[k] = lognormal_rng(0., 1.);
        mu[k]    = normal_rng(10., sigma[k]);
    }
    v[K] = 1.;
    sigma[K] = lognormal_rng(0., 1.);
    mu[K]    = normal_rng(10., sigma[K]);
    {
        real cprod = 1.;
        for (k in 1:K) {
            weights[k] = v[k] * cprod;
            cprod *= (1. - v[k]);
        }
    }
    for (n in 1:N) {
        z[n] = categorical_rng(weights); 
        y[n] = normal_rng(mu[z[n]], sigma[z[n]]);
    }
}
