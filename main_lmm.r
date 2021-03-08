#
# Main script for fitting the linear mixed model
# (Script for ARVO submission and beyond)
# Last Updated: Mar. 8th, 21
# Bayesian approximate kernel regression model with modifications
library(MASS)
library(Matrix)
source("main_processing.r")
utils <- new.env()
source("utilities.r", local = utils)

# 1. We define a norm informed by the spatial structure
# Involves a length scale tuning parameter (for us it should hardly matter)
sparse_scale <- 10 * (2 * pi) / (QUADRANT_NO * 4) # norm > this, sparse
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
## A full eigendecomposition is needed
## (Not generalizable when N is big?)
Kz_ed <- eigen(Kz_hat)
min(which((cumsum(Kz_ed$values) / sum(Kz_ed$values)) > .9)) # 57
U_hat <- Kz_ed$vectors[, 1:50]
diag_hat <- Kz_ed$values[1:50]

# 4. Map back to original linear bases are easily obtained
# by using the M-P inverse (but NOT the only way to project)
beta_map <- ginv(z_train, tol = .Machine$double.eps^.5) %*% U_hat

# 5. The great thing about all this is that from hereon,
# everything boils back down to a multivariate linear model!
# Even pure R doesn't scale too terribly here for Gibbs sampling...
Qstar <- dim(K2)[2]
R <- 50
mu <- colMeans(y_train)
Sig0 <- exp(-K2[knots_inds_y, ]^2 / 4)
Ru0 <- chol(Sig0)
Sig0inv <- tcrossprod(solve(Ru0))
Siginv <- rWishart(1, Qstar + 1, tcrossprod(solve(Ru0)))[, , 1]
Ruinv <- chol(Siginv)
Theta <- diag_hat * t(backsolve(Ruinv, matrix(rnorm(Qstar * R), Qstar)))
sig2inv <- 1.
gamma <- sqrt(2.)
Nu <- t(backsolve(Ruinv, matrix(rnorm(Qstar * N), Qstar)))
iter <- 1000
accepted <- numeric(iter)
for (i in seq(iter)) {
    mask <- exp(-.5 * gamma^-2 * K2^2)
    Theta_prec_right <- sig2inv * crossprod(mask) + Siginv
    Theta_ru_right <- chol(Theta_prec_right)
    Theta_mean <- sig2inv * (crossprod(U_hat, sweep(y_train, 2, mu) - tcrossprod(Nu, mask))) %*% mask # not yet!
    Theta_mean <- forwardsolve(t(Theta_ru_right), t(Theta_mean))
    Theta_mean <- (1. + 1 / diag_hat)^-1 * t(backsolve(Theta_ru_right, Theta_mean))
    Theta <- (1. + 1 / diag_hat)^-1 * 
        t(backsolve(Theta_ru_right, matrix(rnorm(Qstar * R), Qstar))) + # random term
        Theta_mean # mean term
    #
    Nu_mean <- sig2inv * (sweep(y_train, 2, mu) - tcrossprod(U_hat %*% Theta_mean, mask)) %*% mask
    Nu <- t(backsolve(Theta_ru_right, matrix(rnorm(Qstar * N), Qstar)))
    #
    remainder <- y_train - tcrossprod(U_hat %*% Theta_mean + Nu, mask)
    mu <- colMeans(remainder) + rnorm(Q) / (sig2inv * N)
    sig2inv <- rgamma(1, .5 * N * Q, .5 * sum(diag(crossprod(remainder)))) # improper prior
    #
    Sig_post <- Sig0inv + crossprod(diag_hat * Theta, Theta) + crossprod(Nu)
    Siginv <- rWishart(1, Qstar + 1 + .5 * (N + R), Sig_post)[, , 1]
    #
    # M-H proposal is needed for gamma
    gamma_prop <- exp(log(gamma) + rnorm(1, 0, 1.5)) # need pilot 
    mask_prop <- exp(-.5 * gamma_prop^-2 * K2^2)
    ratio <- sum(dnorm(y_train, tcrossprod(U_hat %*% Theta_mean + Nu, mask_prop), sig2inv, log = T)) - 
        sum(dnorm(y_train, tcrossprod(U_hat %*% Theta_mean + Nu, mask), sig2inv, log = T)) + 
        dgamma(gamma_prop, 3, 2, log = T) - dgamma(gamma, 3, 2, log = T)
    if (log(runif(1)) <= ratio) {
        gamma <- gamma_prop
        accepted[i] <- 1
    }
    else accepted[i] <- 0
}
