#
# Main script for fitting the linear mixed model
# (Script for ARVO submission paper)
# Last Updated: Feb. 24th, 21
#
dataset <- readRDS("data/dataset.rds")
centering <- readRDS("data/stats.rds")
N <- round(nrow(dataset$Z) * .7)
J <- max(dataset$id)
QUADRANT_NO <- ncol(dataset$Z) / 4
P <- 3 * QUADRANT_NO
Q <- ncol(dataset$Y)
L <- 2 # Pretty arbitrary at this point
# Knot selection: default is uniformly spaced, dense
Nknots1 <- 48
Nknots2 <- 16
full1  <- seq(P)
full2  <- as.matrix(expand.grid(1:8, 1:8))
knots1 <- seq(1, P, by = 12) # x
knots2_inds <- c(10, 12, 13, 15, 26, 28, 29, 31,
                 34, 36, 37, 39, 50, 52, 53, 55)
knots2 <- full2[knots2_inds, ] # s
# Cumbersome thing to do: since rstan does not provide
# access to kernel functions other than quadratic exponential,
# I need to compute sqrt of pairwise distances and embed
# coordinates in R^D for a "de facto" exponential kernel
D1 <- matrix(0, P, P)
D1[lower.tri(D1)] <- dist(full1) * (2 * pi / QUADRANT_NO / 4)
D1[upper.tri(D1)] <- t(D1)[upper.tri(D1)]
D2 <- matrix(0, Q, Q)
D2[lower.tri(D2)] <- dist(full2)
D2[upper.tri(D2)] <- t(D2)[upper.tri(D2)]
K1 <- D1[, knots1]
K2 <- D2[, knots2_inds]
M1 <- K1[knots1, ]
M2 <- K2[knots2_inds, ]
M1 <- sweep(M1[, 1] - M1, 2, M1[1,], "+") / 2
M2 <- sweep(M2[, 1] - M2, 2, M2[1,], "+") / 2
X1 <- with(eigen(M1), sweep(vectors, 2, values^.5, "*"))[, -Nknots1]
X2 <- with(eigen(M2), sweep(vectors, 2, values^.5, "*"))[, -Nknots2]
# Forming a kernel matrix
rho_k1 <- 4 * (2*pi / QUADRANT_NO / 4)
rho_k2 <- 1 * (1 + 1)^.5
K1 <- exp(-K1 / rho_k1)
K2 <- exp(-K2 / rho_k2)
# Many coefficients of K1 are numerically indistinguishable from zero
K1[K1 < sqrt(.Machine$double.eps)] <- 0.0
# Actual data (pre-processed: centering + normalization?)
zkeep_inds <- c(seq(QUADRANT_NO * 3 + 1, ncol(dataset$Z)),
                seq(1, QUADRANT_NO * 2))
z <- dataset$Z[, zkeep_inds]
y <- dataset$Y
# Imputation (I don't want to bother with missing values--just now)
mis_inds <- which(is.na(y), arr.ind = T)
for (i in 1:nrow(mis_inds)) {
    r <- mis_inds[i, 1]
    j <- mis_inds[i, 2]
    neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
    neighbors <- neighbors[neighbors > 0 & neighbors < Q]
    nmeans <- mean(y[r, neighbors], na.rm = T)
    y[r, j] <- ifelse(is.nan(nmeans), mean(y[r, ], na.rm = T), nmeans)
}
z_train <- z[1:N, ]
y_train <- y[1:N, ]
z_train <- sweep(z_train, 2, centering$cp["q5", zkeep_inds])
y_train <- sweep(y_train, 2, centering$m["q5", ])
#
data <- list(
    N = N, 
    P = P,
    Q = Q, 
    Nknots1 = Nknots1,
    Nknots2 = Nknots2,
    L = L, 
    x = X1, 
    s = X2, 
    z = z_train / 30, 
    y = y_train / 12, 
    K1 = K1,
    K2 = K2
)
library(rstan)
rstan_options(auto_write = TRUE)
m <- stan_model(file = 'constrained_nnet.stan')
options(mc.cores = parallel::detectCores())
o <- optimizing(
    m, data = data,
    algorithm = "LBFGS"
)   
o_points <- o$par
o_names  <- names(o_points)
xi <- o_points[grep("xi", o_names)]
weights <- matrix(o_points[grep("weights", o_names)], Nknots2)
biases <- o_points[grep("biases", o_names)]
intercepts <- o_points[grep("intercept", o_names)]
alpha <- o_points["alpha2"]
beta <- o_points["beta2"]
tau <- o_points["tau2"]
rho <- o_points["rho2"]^2 * 2 
sigma <- o_points["sigma2"]
# Can it simulate realistic images?
Nsamples <- 100
ypred <- matrix(0, Nsamples, Q)
for (n in 1:Nsamples) {
    eta <- beta * tanh(biases + weights %*% xi / alpha) + 
        t(chol(tau^2 * exp(-D2[knots2_inds, knots2_inds] / rho))) %*% 
        matrix(rnorm(Nknots2))
    ypred[n, ] <- intercepts + K2 %*% eta + sigma * matrix(rnorm(Q))
}
ypred <- ypred * 12
y_test <- y[-(1:N), ]
y_test <- sweep(y_test, 2, centering$m["q5", ])

png("images/sample_realizations.png", 1080, 1080, res = 150)
par(mfrow = c(1, 2))
boxplot(y_test, ylim = c(-25, 105), main = "Held-out data")
boxplot(ypred, ylim = c(-25, 105), main = "100 Samples from Model")
dev.off()

plots <- new.env()
source("plots.r", local = plots)
cov_eta <- tau^2 * exp(-D2[knots2_inds, knots2_inds] / rho)
cov_y   <- K2 %*% cov_eta %*% t(K2) + diag(sigma^2, Q)
library(reshape2)
g1 <- ggplot(144 * melt(cov_eta)) + 
    geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    labs(x = "", y = "", fill = "")
g2 <- ggplot(144 * melt(cov_y)) + 
    geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    labs(x = "", y = "", fill = "")
png("images/covariance_estimates.png", 1080, 720, res = 150)
plots$multiplot(g1, g2, cols = 2)
dev.off()

png("images/latent_visual.png", 1080, 1080, res = 150)
latents <- c(biases + weights %*% xi)
plot(latents, beta * tanh(latents / alpha),
     xlab = "Bias + Weights * xi", ylab = "After scaled tanh",
     pch = 19) 
dev.off()
