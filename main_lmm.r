# Main script for fitting the linear mixed model
utils <- new.env()
source("utilities.r", local = utils)
source("main_processing.r")
source("my_gibbs.r")
mcmc_opts <- list(N = 10000, delta = 1)
prior_Prec <- solve(exp(-distMat^2 * .25)[knots, knots])
kernel_choices <- c(1, 2, 4, 6)
results <- vector(mode = "list", length = 4)
for (iter in 1:4) {
    Kmat <- utils$form_kernel_matrix(z_train, kernel_choices[iter])
    K_svd <- svd(Kmat)
    cumprop <- cumsum(K_svd$d^2) / sum(K_svd$d^2)
    cutoff  <- min(which(cumprop > .99))
    U_svd <- K_svd$u[, 1:cutoff]
    d_svd <- K_svd$d[1:cutoff]
    mcmc_opts$delta <- .1
    results[[iter]] <-
        my_gibbs(y_train, U_svd, d_svd, prior_Prec, distMat, knots, mcmc_opts)
}
saveRDS(results, "fit_results.rds")