# Data pre-processing before feeding into the model
# 1-1. Permuting rows of the dataset before train-test split
# (Common ML practice)
dataset <- readRDS("data/dataset.rds")
centering <- readRDS("data/stats.rds")
utils <- new.env()
source("utilities.r", local = utils)

# Permute rows for possible time effect
set.seed(2021-3-8)
permuted_rows <- sample(nrow(dataset$Z))
z <- dataset$Z[permuted_rows, ]
y <- dataset$Y[permuted_rows, ]
labels <- dataset$group[permuted_rows]
# Partition the data on group (patient) lvl
J <- max(dataset$id) # no. of patients
jj <- dataset$id[permuted_rows]
J_train <- round(J * .7)
jj_train <- jj <= J_train
jj_test  <- jj > J_train
N_train  <- sum(jj_train)
N_test   <- sum(jj_test)
N        <- N_train + N_test
# 1-2. Exclusion of the nasal quadrant (not relevant)
QUADRANT_NO <- ncol(dataset$Z) / 4
P <- 3 * QUADRANT_NO
Q <- ncol(dataset$Y)
zkeep_inds <- c(
    seq(QUADRANT_NO * 3 + 1, QUADRANT_NO * 4),
    seq(1, QUADRANT_NO * 2)
)
z <- z[, zkeep_inds]
# 1-3. "Centering" (though deviations are in fact arbitrary scale change)
z <- sweep(z, 2, centering$cp["q5", zkeep_inds])
y <- sweep(y, 2, centering$m["q5", ])
# 1-4. Rescale to milimeters (For more stability)
z <- z / 1000
y <- y / 1000
# 1-5. Resolution change of cpRNFL image
# (Since macula image is relatively so coarse)
# average by pairs => each location ~.9 angle apart
z <- .5 * z %*% (diag(1, P / 2) %x% c(1, 1))
P <- P / 2
# 1-6. Train-test set split
z_train <- z[jj_train, ]
y_train <- y[jj_train, ]
z_test  <- z[jj_test, ]
y_test  <- y[jj_test, ]

# 2. Knot selection over surface of macula image
# (Design set of knots is itself a tuning parameter)
Nknots_y <- 16
full_y  <- as.matrix(expand.grid(1:8, 1:8))
knots_inds_y <- c(11, 14, 18, 20, 21, 23, 27, 30,
                  35, 38, 42, 44, 45, 47, 51, 54)
knots_y <- full_y[knots_inds_y, ]
D2 <- as.matrix(dist(full_y))
K2 <- D2[, knots_inds_y]

# 3. distance matrix of the cpRNFL
# to be used later as valid weighting
D1 <- as.matrix(dist(1:P * (2*pi) / P))

# 4. Finding out missing values
# We don't need to impute them (if so, only for convenience sake)
mis_inds <- which(is.na(y_train), arr.ind = T)

######
library(Matrix)
# Involves a length scale tuning parameter (for us it should hardly matter)
sparse_scale <- 15 * (2 * pi) / (QUADRANT_NO * 4) # distance > this, sparse
Cz <- as(utils$wendland_c2(D1, sparse_scale), "symmetricMatrix")
ru_Cz <- chol(Cz)
Dz <- as.matrix(dist(t(forwardsolve(t(ru_Cz), t(z_train)))))^2

# Tuning parameter whose importance we should investigate!
lambda <- 10 # can vary 
# The approximate kernel is defined by random Fourier bases
Phi <- matrix(0, P, N_train) # approximate eigenfunctions
for (p in seq(P)) {
    freq <- as.numeric(sqrt(2 * lambda) * forwardsolve(t(ru_Cz), rnorm(P)))
    phase <- runif(1) * (2 * pi)
    Phi[p, ] <- sqrt(2 / P) * cos(t(z_train %*% freq) + phase)
}
Kz_hat <- crossprod(Phi) # approximate kernel
Kz_ed <- eigen(Kz_hat) # eigendecomposition (which we can afford)
incl_ind <- min(which((cumsum(Kz_ed$values) / sum(Kz_ed$values)) > .9))
U_hat <- Kz_ed$vectors[, 1:incl_ind]
diag_hat <- Kz_ed$values[1:incl_ind]

# Gibbs sampler for a standard, conjugate multivariate linear model
# Initialization of constants and variable containers
Qstar <- dim(K2)[2]
R <- incl_ind
Sig0 <- exp(-K2[knots_inds_y, ]^2 / 4) # sqrt(2) = diagonal 1 block
Ru0 <- chol(Sig0)
Sig0inv <- tcrossprod(solve(Ru0))
Siginv <- rWishart(1, Qstar + 1, tcrossprod(solve(Ru0)))[, , 1]
Ruinv <- chol(Siginv)
Theta <- diag_hat * t(backsolve(Ruinv, matrix(rnorm(Qstar * R), Qstar)))
sig2inv <- 1.
gamma <- sqrt(2.)
mask <- exp(-.5 * gamma^-2 * K2^2)
Nu <- t(backsolve(Ruinv, matrix(rnorm(Qstar * N_train), Qstar)))
## (Initialize missing values of response)
for (i in 1:nrow(mis_inds)) {
    r <- mis_inds[i, 1]
    j <- mis_inds[i, 2]
    neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
    neighbors <- neighbors[neighbors > 0 & neighbors < Q]
    nmeans <- mean(y_train[r, neighbors], na.rm = T)
    y_train[r, j] <- ifelse(is.nan(nmeans), mean(y_train[r, ], na.rm = T), nmeans)
}
mu <- colMeans(y_train)
## (MCMC variables)
iter <- 10000
mh_sd <- 1.7
accepted <- numeric(iter)
lps <- matrix(0., N_train, iter)
## (array storing all relevant parameters)
pars <- matrix(0, length(mu) + length(Theta) + Qstar * (Qstar + 1) / 2 + 2 + nrow(mis_inds), iter)
## (For indexing)
pars_size <- c(0, length(mu), length(Theta), Qstar * (Qstar + 1) / 2, 1, 1, nrow(mis_inds), nrow(pars))
pars_size <- rev(rev(cumsum(pars_size))[-1])
pars_indices <- vector(mode = "list", length = length(pars_size) - 1)
for (i in 1:(length(pars_size) - 1)) {
    pars_indices[[i]] <- seq(pars_size[i] + 1, pars_size[i + 1])
}
