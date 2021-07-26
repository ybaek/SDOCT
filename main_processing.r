# Data pre-processing before feeding into the model
dataset <- readRDS("data/dataset.rds")
centering <- readRDS("data/stats.rds")
utils <- new.env()
source("utilities.r", local = utils)

# 1-1. Permuting rows of the dataset before train-test split
set.seed(2021-3-8)
permuted_rows <- sample(nrow(dataset$Z))
z <- dataset$Z[permuted_rows, ]
y <- dataset$Y[permuted_rows, ]
labels <- dataset$group[permuted_rows]
J <- max(dataset$id) # no. of patients
jj <- dataset$id[permuted_rows]
J_train <- round(J * .7)
jj_train <- jj <= J_train
jj_test  <- jj > J_train
N_train  <- sum(jj_train)
N_test   <- sum(jj_test)
N        <- N_train + N_test

# 1-2. Exclusion of the nasal quadrant
QUADRANT_NO <- ncol(dataset$Z) / 4
P <- 3 * QUADRANT_NO
Q <- ncol(dataset$Y)
zkeep_inds <- c(
    seq(QUADRANT_NO * 3 + 1, QUADRANT_NO * 4),
    seq(1, QUADRANT_NO * 2)
)
z <- z[, zkeep_inds]

# 1-3. "Centering" and rescaling to mm (to a more interpretable scale)
z <- sweep(z, 2, centering$cp["q5", zkeep_inds])
y <- sweep(y, 2, centering$m["q5", ])
z <- z / 1000
y <- y / 1000

# 1-4. Resolution downscaling of cpRNFL image
# Average by pairs => each location ~.9 angle apart
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
knots <- c(11, 14, 18, 20, 21, 23, 27, 30,
           35, 38, 42, 44, 45, 47, 51, 54)
distMat <- as.matrix(dist(full_y))

# [OUTDATED]
# 3. distance matrix of the cpRNFL
# to be used later as valid weighting

# 4. Finding out missing values (but not impute them)
# dim(which(is.na(y)), arr.ind = TRUE)
mis_inds <- which(is.na(y_train), arr.ind = T)