# This package ns must be attached
library(GIGrvg)
source("./adj_grid.R")
source("./utilities.R")
input <- readRDS("./data/fit_data.Rds")
A <- adj_grid(8)

# Pre-processing datasets
hyper <- list(rho = .99, beta_var0 = 100, c0 = 1 / 3)
hyper$prec_y <- hyper$rho * input$Lap +
    diag(1.0 - hyper$rho, 64)
l <- 12 / 768 * (2 * pi)
hyper$cov_x <- wendland_c2(input$distMat, l)
mcmc <- list(I = 2000, burnin = 500)
data <- input[c("small_inds", "zmeans", "ymeans", "labels")]

# Compile the Gibbs sampler routine
library(Rcpp)
sourceCpp("samplers_lmm.cpp")
# Update: a giant loop for K-fold
for (k in 1:length(input$train)) {
    data$X <- input$Z_s[input$train[[k]], ]
    data$Y <- input$Y_s[input$train[[k]], ]
    # ID's (which are just row numbers)
    data$ids <- seq(1, nrow(data$Y))
    data$mis_inds <- which(is.na(data$Y), arr.ind = TRUE)
    data$mis_adjs <- apply(A[unique(data$mis_inds)[,2],], 1, function(x) which(!!x))
    data$Y[data$mis_inds] <- rnorm(dim(data$mis_inds)[1]) * 3 # Initialize the missing values
    data$n_ns <- sapply(data$mis_adjs, length) # no. of neighbors for each missing location
    data$neighbors <- unlist(data$mis_adjs) # all neighbors of missing locations ordered by patients
    inits <- list(beta = c(coef(lm(data$Y ~ data$X - 1))),
                sig2inv = 1,
                c2 = rep(.25^2, length(data$small_inds)),
                gamma = array(0, dim = c(ncol(data$X), ncol(data$Y), length(data$ids))),
                psi = rWishart(1, ncol(data$X) + 1, solve(hyper$cov_x))[, , 1],
                theta = matrix(0, nrow(data$Y), ncol(data$Y)))
    inits$theta <- (data$Y - data$X %*% matrix(inits$beta, ncol(data$X))) *
        matrix(rnorm(nrow(data$Y)*ncol(data$Y)), nrow(data$Y))
    chain <- mainSampler_lmm(data, inits, hyper, mcmc)
    saveRDS(chain, file = paste0("objects/lmmfit", k, ".rds"))
}
