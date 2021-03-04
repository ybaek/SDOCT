#
# Main script for fitting the linear mixed model
# (Script for ARVO submission paper)
# Last Updated: Mar. 3rd, 21
#
dataset <- readRDS("data/dataset.rds")
centering <- readRDS("data/stats.rds")
N <- round(nrow(dataset$Z) * .7)
# K <- 6
J <- max(dataset$id)
QUADRANT_NO <- ncol(dataset$Z) / 4
P <- 3 * QUADRANT_NO
Q <- ncol(dataset$Y)
zkeep_inds <- c(
    seq(QUADRANT_NO * 3 + 1, QUADRANT_NO * 4),
    seq(1, QUADRANT_NO * 2)
)
# Knot selection: default is uniformly spaced, dense
Nknots_z <- P / 8
Nknots_y <- 16
full_z  <- 1:P
full_y  <- as.matrix(expand.grid(1:8, 1:8))
knots_inds_z <- seq(1, P, by = 8)
knots_inds_y <- c(10, 12, 13, 15, 26, 28, 29, 31,
                 34, 36, 37, 39, 50, 52, 53, 55)
knots_z <- full_z[knots_inds_z] # x
knots_y <- full_y[knots_inds_y, ] # s
# Cumbersome thing to do: since rstan does not provide
# access to kernel functions other than quadratic exponential,
# I need to compute sqrt of pairwise distances and embed
# coordinates in R^D for a "de factor" exponential kernel
D1 <- as.matrix(dist(full_z) * 2 * pi / 768)
D2 <- as.matrix(dist(full_y))
K1 <- D1[, knots_inds_z]
K2 <- D2[, knots_inds_y]
M1 <- K1[knots_inds_z, ]
M2 <- K2[knots_inds_y, ]
M1 <- sweep(M1[, 1] - M1, 2, M1[1,], "+") / 2
M2 <- sweep(M2[, 1] - M2, 2, M2[1,], "+") / 2
X1 <- with(eigen(M1), sweep(vectors, 2, values^.5, "*"))[, -Nknots_z]
X2 <- with(eigen(M2), sweep(vectors, 2, values^.5, "*"))[, -Nknots_y]
#
# Imputation (I don't want to bother with missing values--just now)
y <- dataset$Y
mis_inds <- which(is.na(y), arr.ind = T)
for (i in 1:nrow(mis_inds)) {
    r <- mis_inds[i, 1]
    j <- mis_inds[i, 2]
    neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
    neighbors <- neighbors[neighbors > 0 & neighbors < Q]
    nmeans <- mean(y[r, neighbors], na.rm = T)
    y[r, j] <- ifelse(is.nan(nmeans), mean(y[r, ], na.rm = T), nmeans)
}
#
z_train <- dataset$Z[1:N, zkeep_inds]
y_train <- y[1:N, ]
z_train <- sweep(z_train, 2, centering$cp["q5", zkeep_inds])
y_train <- sweep(y_train, 2, centering$m["q5", ])
#
data <- list(
    N = N, 
    P = P,
    Q = Q, 
    Nknots_z = Nknots_z,
    Nknots_y = Nknots_y,
    kk_z = knots_inds_z,
    kk_y = knots_inds_y,
    coords_z = full_z * 2 * pi / 768,
    coords_y = as.matrix(expand.grid(seq(Q^.5), seq(Q^.5))),
    x_tf = X1,
    s_tf = X2,
    z = z_train / 10, 
    y = y_train / 10
)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
m <- stan_model(file = 'constrained_nnet.stan')
o <- optimizing(
    m, data = data,
    algorithm = "LBFGS",
    verbose = TRUE
)   
o_points <- o$par
saveRDS(o_points, "estimates.rds")

o_names  <- names(o_points)
W <- matrix(o_points[grep("W", o_names)], P)
b <- o_points[grep("b\\[", o_names)]
intercept <- o_points[grep("intercept", o_names)]
beta <- o_points["beta"]
alpha <- o_points["alpha"]
sigma <- o_points["sigma"]
tau <- o_points["tau"]
lambda <- o_points["lambda"]
rho <- o_points["rho"]
kernel <- exp(-D^2 / (2 * lambda^2))[, knots_inds]
latent_cov <- tau^2 * exp(-D[knots_inds, knots_inds] / rho)