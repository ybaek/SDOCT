# Simulations playing with dataset centering scheme
# and various dimension reduction techniques
dataset <- readRDS("data/dataset.rds")
center_stats <- readRDS("data/stats.rds")
# Let's play around with the 5th percentile
dataset$Z <- sweep(dataset$Z, 2, center_stats[["cp"]]["q5", ])
dataset$Y <- sweep(dataset$Y, 2, center_stats[["macula"]]["q5", ])
# Exclude the nasal quadrant and reorder to be continuous
# direction: from inferior to superior
QUAD_NO <- dim(dataset$Z)[2] / 4
dataset$Z <- dataset$Z[, -seq(QUAD_NO * 2 + 1, QUAD_NO * 3)]
dataset$Z <- cbind(
    dataset$Z[, seq(QUAD_NO * 2 + 1, 3 * QUAD_NO)],
    dataset$Z[, seq(1, QUAD_NO * 2)]
)
# Calculating the Euclidean distances
## For cp region, there is a unit angle dividing 2*pi rad into 768 regions
UNIT_ANGLE <- 2 * pi / (QUAD_NO * 4)
dist_Z <- array(0.0, dim = rep(dim(dataset$Z)[2], 2))
dist_Z[lower.tri(dist_Z)] <- dist(UNIT_ANGLE * seq(dim(dataset$Z)[2]))
dist_Z[upper.tri(dist_Z)] <- t(dist_Z)[upper.tri(dist_Z)]
## For macula, 8 x 8 image with equally spaced pixels
dist_Y <- array(0.0, dim = rep(dim(dataset$Y)[2], 2))
dist_Y[lower.tri(dist_Y)] <- dist(cbind(rep(1:8, each = 8), rep(1:8, 8)))
dist_Y[upper.tri(dist_Y)] <- t(dist_Y)[upper.tri(dist_Y)]

# Knot selection strategy
# Goal: reducing the dimension of observation domains
# Default: Dense set of knots "roughly covering" the area
N <- 3
D  <- 5
M1 <- 3
knot_locs_Z <- seq(1, dim(dataset$Z)[2], by = M1)
bw_rbf_Z <- 0.1
bw_ou_Z  <- 0.1 / M1
## RBF kernel with half the bandwidth is the "limit" of same bw
## Ornstein-Uhlenbeck process (special case of Matern)
kz <- exp(-dist_Z[, knot_locs_Z]^2 / (2.0 * bw_rbf_Z))
## OU kernel for Z yields a VERY ill-conditioned matrix
root_covz <- chol(exp(-dist_Z[knot_locs_Z, knot_locs_Z] / bw_ou_Z))
# Prior predictives when convolving white noise processes (Higdon, 1999)
xi <- matrix(1 + rnorm(N * D), N)
weights1 <- matrix(rnorm(dim(dataset$Z)[2] / M1 * D) + 1, ncol = D)
latent_z <- tanh(tcrossprod(xi, weights1)) +
    matrix(rnorm(N * dim(dataset$Z)[2] / M1), N) %*% root_covz
pp_z <- tcrossprod(latent_z, kz) + matrix(rnorm(N * dim(dataset$Z)[2]), N)

png("images/simulation_Z.png", width = 5, height = 5, units = "in", res = 150)
.env$nice_par()
plot(pp_z[1, ], type = "l", ylim = range(pp_z),
     xlab = "ITS", ylab = "",
     main = "RBF bw 0.2, latent OU process bw 0.03, tanh")
lines(rep(latent_z[1, ], each = M1), col = 2)
for (i in 2:N) {
    lines(pp_z[i, ])
    lines(rep(latent_z[i, ], each = M1), col = 2, pch = 19)
}
dev.off()