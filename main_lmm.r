# Main script for fitting the linear mixed model
# (Script for ARVO submission and beyond)
# Last Updated: Mar. 8th, 21
# Bayesian approximate kernel regression model with modifications
source("main_processing.r")

## (Actual sampling)
for (i in seq(iter)) {
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
    Nu <- Nu_mean + t(backsolve(Theta_ru_right, matrix(rnorm(Qstar * N_train), Qstar)))
    #
    Sig_post <- Sig0inv + crossprod(diag_hat * Theta, Theta) + crossprod(Nu)
    Siginv <- rWishart(1, Qstar + 1 + .5 * (N_train + R), Sig_post)[, , 1]
    #
    remainder <- y_train - tcrossprod(linpred + Nu, mask)
    mu <- colMeans(remainder) + rnorm(Q) / (sig2inv * N_train)
    sig2inv <- rgamma(1, .5 * N_train * Q, .5 * sum(diag(crossprod(remainder)))) # improper prior
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
    mask <- exp(-.5 * gamma^-2 * K2^2)
    #
    pars[, i] <- c(mu, c(Theta), c(Siginv[lower.tri(Siginv, diag = T)]), sig2inv, gamma, y_train[mis_inds])
    #
    lps[, i] <- rowSums(
        dnorm(
            y_train,
            sweep(tcrossprod(linpred + Nu, mask), 2, mu, "+"),
            sig2inv^-.5,
            log = TRUE
        )
    )
}
burnin <- 2000
pars <- pars[, (burnin + 1):ncol(pars)]
lps  <- lps[, (burnin + 1):ncol(lps)]

saveRDS(pars, "data/samples5.rds")
saveRDS(pars_indices, "data/indices5.rds")
saveRDS(lps, "data/lps5.rds")
