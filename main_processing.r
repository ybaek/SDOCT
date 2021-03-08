# Data pre-processing before feeding into the model
# 1-1. Permuting rows of the dataset before train-test split
# (Common ML practice)
dataset <- readRDS("data/dataset.rds")
centering <- readRDS("data/stats.rds")

set.seed(2021-3-8)
permuted_rows <- sample(nrow(dataset$Z))
z <- dataset$Z[permuted_rows, ]
y <- dataset$Y[permuted_rows, ]
N <- round(nrow(dataset$Z) * .7) # training set size
J <- max(dataset$id) # no. of patients
jj <- dataset$id[permuted_rows]
labels <- dataset$group[permuted_rows]
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
z_train <- z[1:N, ]
y_train <- y[1:N, ]

# 2. Knot selection over surface of macula image
# (Design set of knots is itself a tuning parameter)
Nknots_y <- 16
full_y  <- as.matrix(expand.grid(1:8, 1:8))
knots_inds_y <- c(10, 12, 13, 15, 26, 28, 29, 31,
                 34, 36, 37, 39, 50, 52, 53, 55)
knots_y <- full_y[knots_inds_y, ]
D2 <- as.matrix(dist(full_y))
K2 <- D2[, knots_inds_y]

# 3. distance matrix of the cpRNFL
# to be used later as valid weighting
D1 <- as.matrix(dist(1:P * (2*pi) / P))

# 4. Finding out missing values
# We don't need to impute them (if so, only for convenience sake)
mis_inds <- which(is.na(y), arr.ind = T)
# for (i in 1:nrow(mis_inds)) {
#     r <- mis_inds[i, 1]
#     j <- mis_inds[i, 2]
#     neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
#     neighbors <- neighbors[neighbors > 0 & neighbors < Q]
#     nmeans <- mean(y[r, neighbors], na.rm = T)
#     y[r, j] <- ifelse(is.nan(nmeans), mean(y[r, ], na.rm = T), nmeans)
# }
