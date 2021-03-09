# Main script for fitting the linear mixed model
# (Script for ARVO submission and beyond)
# Last Updated: Mar. 8th, 21
# Bayesian approximate kernel regression model with modifications
library(Matrix)
source("main_processing.r")
utils <- new.env()
source("utilities.r", local = utils)

# 1. We define a norm informed by the spatial structure
# Involves a length scale tuning parameter (for us it should hardly matter)
sparse_scale <- 10 * (2 * pi) / (QUADRANT_NO * 4) # distance > this, sparse
Cz <- as(utils$wendland_c2(D1, sparse_scale), "symmetricMatrix")
pairs <- matrix(0, P, N * (N - 1) / 2) # N * (N-1)/2 pairwise differences
count <- 1
for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
        pairs[, count] <- z_train[i, ] - z_train[j, ]
        count <- count + 1
    }
}
wns <- colSums(pairs * solve(Cz, pairs))
Dz <- matrix(0, N, N)
Dz[lower.tri(Dz)] <- wns
Dz[upper.tri(Dz)] <- t(Dz)[upper.tri(Dz)]

# 2. We can now form the EXACT kernel matrix
# Note this involves another choice of length scale parameter
lambda <- 10 # can vary 
Kz <- exp(-Dz * lambda)

# 3. The approximate kernel is defined by random Fourier bases
# For details, see the paper itself
rl_Cz <- t(chol(Cz))
Phi <- matrix(0, P, N) # approximate eigenfunctions
for (p in seq(P)) {
    freq <- as.numeric(sqrt(2 * lambda) * forwardsolve(rl_Cz, rnorm(P)))
    phase <- runif(1) * (2 * pi)
    Phi[p, ] <- sqrt(2 / P) * cos(t(z_train %*% freq) + phase)
}
Kz_hat <- crossprod(Phi) # approximate kernel
Kz_ed <- eigen(Kz_hat) # eigendecomposition (which we can afford)
incl_ind <- min(which((cumsum(Kz_ed$values) / sum(Kz_ed$values)) > .9))
U_hat <- Kz_ed$vectors[, 1:incl_ind]
diag_hat <- Kz_ed$values[1:incl_ind]

# 4. Map back to original linear bases are easily obtained
# Since we have full-rank design, there is NO ambiguity contrary to the original paper
beta_map <- solve(crossprod(z_train), crossprod(z_train, U_hat))

# 5. Gibbs sampler for a standard, conjugate multivariate linear model
# Initialization of constants and variable containers
Qstar <- dim(K2)[2]
R <- incl_ind
mu <- colMeans(y_train)
Sig0 <- exp(-K2[knots_inds_y, ]^2 / 4) # sqrt(2) = diagonal 1 block
Ru0 <- chol(Sig0)
Sig0inv <- tcrossprod(solve(Ru0))
Siginv <- rWishart(1, Qstar + 1, tcrossprod(solve(Ru0)))[, , 1]
Ruinv <- chol(Siginv)
Theta <- diag_hat * t(backsolve(Ruinv, matrix(rnorm(Qstar * R), Qstar)))
sig2inv <- 1.
gamma <- sqrt(2.)
Nu <- t(backsolve(Ruinv, matrix(rnorm(Qstar * N), Qstar)))
## (Initialize missing values of response)
for (i in 1:nrow(mis_inds)) {
    r <- mis_inds[i, 1]
    j <- mis_inds[i, 2]
    neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
    neighbors <- neighbors[neighbors > 0 & neighbors < Q]
    nmeans <- mean(y_train[r, neighbors], na.rm = T)
    y_train[r, j] <- ifelse(is.nan(nmeans), mean(y_train[r, ], na.rm = T), nmeans)
}
## (MCMC variables)
iter <- 10000
mh_sd <- 1.7
accepted <- numeric(iter)
## (array storing all relevant parameters)
pars <- matrix(0, length(mu) + length(Theta) + Qstar * (Qstar + 1) / 2 + 2 + nrow(mis_inds), iter)
colnames(pars) <- c(
    paste0("mu[", seq(Q), "]"),
    paste0("Theta[", seq(R * Qstar), "]"),
    paste0("Siginv[", seq(Qstar * (Qstar + 1) / 2), "]"),
    "sig2inv", "gamma"
)
## (Actual sampling)
for (i in seq(iter)) {
    mask <- exp(-.5 * gamma^-2 * K2^2)
    Theta_prec_right <- sig2inv * crossprod(mask) + Siginv
    Theta_ru_right <- chol(Theta_prec_right)
    Theta_mean <- sig2inv *
        (crossprod(U_hat, sweep(y_train, 2, mu) - tcrossprod(Nu, mask))) %*% mask # not yet!
    Theta_mean <- (1. + 1 / diag_hat)^-1 * t(
        backsolve(Theta_ru_right, forwardsolve(t(Theta_ru_right), t(Theta_mean)))
    )
    Theta <- Theta_mean +
        (1. + 1 / diag_hat)^-1 * 
        t(backsolve(Theta_ru_right, matrix(rnorm(Qstar * R), Qstar)))
    linpred <- U_hat %*% Theta_mean
    #
    Nu_mean <- sig2inv * (sweep(y_train, 2, mu) - tcrossprod(linpred, mask)) %*% mask
    Nu_mean <- t(backsolve(Theta_ru_right, forwardsolve(t(Theta_ru_right), t(Nu_mean))))
    Nu <- Nu_mean + t(backsolve(Theta_ru_right, matrix(rnorm(Qstar * N), Qstar)))
    #
    Sig_post <- Sig0inv + crossprod(diag_hat * Theta, Theta) + crossprod(Nu)
    Siginv <- rWishart(1, Qstar + 1 + .5 * (N + R), Sig_post)[, , 1]
    #
    remainder <- y_train - tcrossprod(linpred + Nu, mask)
    mu <- colMeans(remainder) + rnorm(Q) / (sig2inv * N)
    sig2inv <- rgamma(1, .5 * N * Q, .5 * sum(diag(crossprod(remainder)))) # improper prior
    #
    # Imputation
    y_train[mis_inds] <- sweep(tcrossprod(linpred + Nu, mask), 2, mu, "+")[mis_inds]
    #
    # M-H proposal is needed for gamma
    gamma_prop <- exp(log(gamma) + rnorm(1, 0, mh_sd)) # monitor acceptances
    mask_prop <- exp(-.5 * gamma_prop^-2 * K2^2)
    ratio <- sum(dnorm(y_train, tcrossprod(linpred + Nu, mask_prop), sig2inv, log = T)) - 
        sum(dnorm(y_train, tcrossprod(linpred + Nu, mask), sig2inv, log = T)) + 
        dgamma(gamma_prop, 3, 2, log = T) - dgamma(gamma, 3, 2, log = T)
    if (log(runif(1)) <= ratio) {
        gamma <- gamma_prop
        accepted[i] <- 1
    }
    else accepted[i] <- 0
    #
    pars[, i] <- c(mu, c(Theta), c(Siginv[lower.tri(Siginv, diag = T)]), sig2inv, gamma, y_train[mis_inds])
}
