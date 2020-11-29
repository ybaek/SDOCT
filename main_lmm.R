# This package ns must be attached
library(GIGrvg)
source("./utilities.R")

input <- readRDS("./data/fit_data.Rds")
data <- input[c("ids", "mis_inds", "small_inds", "zmeans", "ymeans", "labels")]
data$X <- input$Z_s
data$Y <- input$Y_s
data$Y[data$mis_inds] <- rnorm(dim(data$mis_inds)[1]) * 3 # Initialize the missing values
data$n_ns <- sapply(input$mis_adjs, length) # no. of neighbors for each missing location
data$neighbors <- unlist(input$mis_adjs) # all neighbors of missing locations ordered by patients
hyper <- list(rho = .99, beta_var0 = 100, c0 = 1/3)
l <- 15 / 768 * (2*pi) 
hyper$prec_y <- hyper$rho * input$Lap + diag(1-hyper$rho, ncol(data$Y))
hyper$prec_x <- solve(wendland_c2(input$distMat, l))

# Passing in initial values of all parameters
inits <- list( beta = c(coef(lm(data$Y ~ data$X - 1))),
               sig2inv = 1, tau2inv = 1, c2 = rep(.25^2, length(data$small_inds)),
               gamma = array(0, dim = c(ncol(data$X), ncol(data$Y), max(data$ids))),
               theta = matrix(0, nrow(data$Y), ncol(data$Y)) )
inits$theta <- (data$Y - data$X %*% matrix(inits$beta, ncol(data$X))) * matrix(rnorm(nrow(data$Y)*ncol(data$Y)), nrow(data$Y))

mcmc <- list(I = 1800, burnin = 300)

# Compile the Gibbs sampler routine
library(Rcpp)
sourceCpp("samplers.cpp")
chain1 <- mainSampler_lmm(data, inits, hyper, mcmc)
saveRDS(chain1, "./objects/samples_lmm1_cpp.Rds")

# See if standardization leads to drastically different results
data$X <- scale(input$Z_s, F, T)
data$Y <- scale(input$Y_s, F, T)
chain2 <- mainSampler_lmm(data, inits, hyper, mcmc)
saveRDS(chain2, "./objects/samples_lmm2_cpp.Rds"
