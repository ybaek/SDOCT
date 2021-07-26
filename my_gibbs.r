my_gibbs <- function(Y, U_svd, d_svd, prior_Prec, distMat, knots, mcmc_opts) {
    # Gibbs sampler specific to our purpose and thus non-portable
    # Mostly standard linear regression, though some tweaks exist
    N_samples <- mcmc_opts$N
    delta <- mcmc_opts$delta
    N <- dim(Y)[1]
    Q <- dim(Y)[2]
    D <- dim(U_svd)[2]
    Qstar <- length(knots)
    mis_inds <- which(is.na(Y), arr.ind = T)
    Nmis <- dim(mis_inds)[1]
    #
    npars <- as.integer(
        # Total no. of things sampled using MCMC
        (D * Qstar) + (Qstar * (Qstar + 1) / 2) + Nmis + 1 + 1
    )
    pars <- matrix(N_samples, npars) # Parameter space
    lps  <- matrix(N_samples, N) # Observable units space
    acceptances <- numeric(N_samples) # For Metropolis alg.
    #
    theta0 <- colMeans(Y, na.rm = T)
    theta <- matrix(rnorm(D * Qstar), D)
    theta_covinv <- rWishart(1, Qstar + 1, solve(prior_Prec))
    sig2inv <- sd(Y, na.rm = T)^-2
    gamma <- rgamma(1, 3, 3)
    # Initialize missing values on the fly,
    # using a spatial neighborhood average
    for (i in seq(Nmis)) {
        r <- mis_inds[i, 1]
        j <- mis_inds[i, 2]
        neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
        neighbors <- neighbors[neighbors > 0 & neighbors < Q]
        nmeans <- mean(y_train[r, neighbors], na.rm = T)
        # Conditional flow for when ALL neighbors are also missing
        Y[r, j] <- ifelse(
            is.nan(nmeans), mean(y_train[r, ], na.rm = T), nmeans
        )
    }
    #
    for (n in seq(N_samples)) {
        convMat <- exp(-distMat[, knots]^2 / gamma)
        theta0_prec <- (1e-4 + sig2inv * N)
        Ymean <- colMeans(Y - U_svd %*% theta %*% t(convMat))
        theta0 <- rnorm(Q, sig2inv * N / theta0_prec * Ymean, theta0_prec^-.5)
        #
        theta_prec <- theta_covinv %x% diag(d_svd^-2) +
            sig2inv * (crossprod(convMat) %x% crossprod(U_svd))
        Ycentered <- sweep(Y, 2, theta0)
        Rtheta <- chol(theta_prec)
        theta_det <- backsolve(Rtheta, forwardsolve(t(Rtheta), c(Ycentered)))
        theta_random <- backsolve(Rtheta, rnorm(D * Qstar))
        theta <- matrix(theta_det + theta_random, D, Qstar)
        #
        post_Prec <- prior_Prec + t(theta) %*% diag(d_svd^-2) %*% theta
        theta_covinv <- rWishart((Qstar + 1) + D / 2, solve(post_Prec))
        #
        meanMat <- sweep(U_svd %*% theta %*% t(convMat), 2, theta0, "+")
        sse <- sum((Y - meanMat)^2)
        sig2inv <- rgamma(1, .5 + N * Q / 2, .5 + sse / 2)
        #
        for (i in seq(Nmis)) {
            r <- mis_inds[i, 1]
            j <- mis_inds[i, 2]
            Y[r, j] <- rnorm(1, meanMat[r, j], sig2inv^-.5)
        }
        #
        gamma_log_prop <- rnorm(1, log(gamma), delta)
        gamma_prop <- exp(gamma_log_prop)
        convMat_prop <- exp(-distMat[, knots]^2 / gamma_prop)
        meanMat_prop <- sweep(U_svd %*% theta %*% t(convMat_prop), 2, theta0, "+")
        ll_sum <- sum(dnorm(Y, meanMat, sig2inv^-.5, log = T)) +
            dgamma(gamma, 3, 3, log = T)
        ll_prop_sum <- sum(dnorm(Y, meanMat_prop, sig2inv^-.5, log = T)) +
            dgamma(gamma_prop, 3, 3, log = T)
        log_r <- min(0, ll_prop_sum - ll_sum)
        log_u <- log(runif(1))
        accept <- (log_u <= log_r)
        gamma <- accept * gamma_prop + (1 - accept) * gamma
        convMat <- exp(-distMat[, knots]^2 / gamma)
        #
        covinds <- lower.tri(theta_covinv, diag = TRUE)
        meanMat <- sweep(U_svd %*% theta %*% t(convMat), 2, theta0, "+")
        pars[n, ] <- c(
            theta0, c(theta), c(theta_covinv[covinds]),
            c(Y[mis_inds]), sig2inv, gamma 
        )
        lps[n, ] <- rowSums(dnorm(Y, meanMat, sig2inv^-.5, log = T))
        acceptances[n] <- accept
    }
    list(
        pars = pars,
        lps  = lps,
        acceptances = acceptances
    )
}