# Pre-processing before fitting the model
df <- readRDS("./data/smallest.Rds")
pseudomean_Y <- readRDS("./data/control_m.Rds")
pseudomean_Z <- readRDS("./data/control_cp.Rds")
# No need to stratify the IDs because no subunits
# ids <- as.integer(factor(df$patientid, levels = sort(unique(df$patientid))))
labels <- df$group
Z <- sweep(as.matrix(df[,(2+1):(2+288)]), 2, pseudomean_Z, "-")
Y <- sweep(as.matrix(df[,(2+288+1):(2+288+64)]), 2, pseudomean_Y, "-")
# Partition the data on PATIENT strata
set.seed(2020 - 12 - 15)
# J <- max(ids)
N <- dim(df)[1]
train <- sort(sample(N, N*.8))
test <- setdiff(1:N, train)
train_rows_bool <- ids %in% train
test_rows_bool <- ids %in% test
# Prior information
source("./adj_grid.R")
source("./prior_distances.R")
## aggregating by 4 for predictors
prior_mean <- prior_mean[seq(1, 288, by = 8), ]
distMat <- distMat[seq(1, 288, by = 8), seq(1, 288, by = 8)]
A <- adj_grid(8)
Lap <- diag(rowSums(A)) - A
small_inds <- which(!prior_mean)
## Scaling and centering data
Z_a <- Z %*% (diag(1, 36) %x% rep(.125,8))
colnames(Z_a) <- colnames(Z)[seq(1, 288, by = 8)]
Z_s <- scale(Z_a, T, F)
Y_s <- scale(Y, T, F)
zmeans <- attr(Z_s, "scaled:center")
ymeans <- attr(Y_s, "scaled:center")
# I've come to realize they are measured in same units, so standardization is not a good idea
# zsds <- attr(Z_s, "scaled:scale")
# ysds <- attr(Y_s, "scaled:scale")

# mis_adjs <- apply(A[unique(mis_inds)[,2],], 1, function(x) which(!!x))

data <- list(Z_s = Z_s, Y_s = Y_s, zmeans = zmeans, ymeans = ymeans,
             small_inds = small_inds, distMat = distMat, Lap = Lap,
             # ids = ids,
             train = train_rows_bool, test = test_rows_bool, labels = labels)
saveRDS(data, file = "./data/fit_data.Rds")
