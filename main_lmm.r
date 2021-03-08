#
# Main script for fitting the linear mixed model
# (Script for ARVO submission paper)
# Last Updated: Feb. 10th, 21
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
# rescale to milimeter
z_train <- z_train / 1000
y_train <- y_train / 1000
#
# New method: start by forming a kernel
# First, we define a norm informed by the spatial structure
library(Matrix)
utils <- new.env()
source("utilities.r", local = utils)
sparse_scale <- 10 * (2 * pi) / (QUADRANT_NO * 4) # norm > , sparse
Cz <- as(utils$wendland_c2(D1, sparse_scale), "symmetricMatrix")
## Compare the distance matrices
Dz1 <- as.matrix(dist(z_train)^2)
## Need N * (N-1)/2 pairwise differences
pairs <- matrix(0, P, N * (N - 1) / 2)
count <- 1
for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
        pairs[, count] <- z_train[i, ] - z_train[j, ]
        count <- count + 1
    }
}
wns <- colSums(pairs * solve(Cz, pairs))
Dz2 <- matrix(0, N, N)
Dz2[lower.tri(Dz2)] <- wns
Dz2[upper.tri(Dz2)] <- t(Dz2)[upper.tri(Dz2)]
# After normalizing, Dz1 and Dz2 are practically indistinguishable...
# This is because the more globally independent they are,
# covariance matrix becomes closer to identity
Kz <- exp(-Dz2 / .1) # can vary
# Approximate kernel: thru random Fourier bases
# see how close the approx is with P bases
Cz_inv <- solve(Cz)
Phi <- matrix(0, P, N)
for (p in seq(P)) {
    freq <- as.numeric(sqrt(2 / .1) * t(chol(Cz_inv)) %*% rnorm(P))
    phase <- runif(1) * (2 * pi)
    Phi[p, ] <- sqrt(2 / P) * cos(t(z_train %*% freq) + phase)
}
Kz_hat <- crossprod(Phi)
Kz_ed <- eigen(Kz_hat)
min(which((cumsum(Kz_ed$values) / sum(Kz_ed$values)) > .9)) # 51
U_hat <- Kz_ed$vectors[, 1:50]
diag_hat <- Kz_ed$values[1:50]
# Inverse map is simply obtained by an M-P inverse
# of the design matrix
library(MASS)
beta_map <- ginv(z_train, tol = .Machine$double.eps^.5) %*% U_hat
