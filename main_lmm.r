# Main script for fitting the linear mixed model
utils <- new.env()
source("utilities.r", local = utils)
source("main_processing.r")
source("my_gibbs.r")
mcmc_opts <- list(N = 10000, delta = 1e-2)
prior_Prec <- solve(exp(-distMat^2 * .25)[knots, knots])
kernel_choices <- c(2, 4, 6)
M <- length(kernel_choices)
results <- vector(mode = "list", length = M)
for (iter in seq(M)) {
    Kmat <- utils$form_kernel_matrix(z_train, kernel_choices[iter])
    K_svd <- svd(Kmat)
    cumprop <- cumsum(K_svd$d^2) / sum(K_svd$d^2)
    cutoff  <- min(which(cumprop > .99))
    U_svd <- K_svd$u[, 1:cutoff]
    d_svd <- K_svd$d[1:cutoff]
    # only a heuristic
    # mcmc_opts$delta <- mcmc_opts$delta / kernel_choices[iter] * 2
    results[[iter]] <-
        my_gibbs(y_train, U_svd, d_svd, prior_Prec, distMat, knots, mcmc_opts)
}
# WAIC comparison
B <- 2000 # burn-in
library(loo)
waic_reports <- lapply(results, function(x) waic(x$pars[-seq(B), ]))
choose_m <- which.min(sapply(waic_reports, function(x) x$waic))
# choose the best one
saveRDS(results[[choose_m]], "data/fit_results.rds")