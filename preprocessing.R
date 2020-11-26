# Pre-processing before fitting the model
df <- readRDS("../data/macula_cross14-17.Rds")
pseudomean_Y <- readRDS("../data/control_m.Rds")
pseudomean_Z <- readRDS("../data/control_cp.Rds")
ids <- as.integer(factor(df$patientid, levels = sort(unique(df$patientid))))
id_list <- lapply(1:max(ids), function(x) which(ids==x))
labels <- df$group
Z <- sweep(as.matrix(df[,(2+1):(2+288)]), 2, pseudomean_Z, "-")
Y <- sweep(as.matrix(df[,(2+288+1):(2+288+64)]), 2, pseudomean_Y, "-")
mis_inds <- which(is.na(Y), arr.ind = TRUE)
# Prior information
source("./adj_grid.R")
source("./prior_distances.R")
## aggregating by 4 for predictors
prior_mean <- prior_mean[seq(1, 288, by = 8), ]
distMat <- distMat[seq(1, 288, by = 8), seq(1, 288, by = 8)]
small_inds <- which(!prior_mean)
A <- adj_grid(8)
mis_adjs <- apply(A[unique(mis_inds)[,2],], 1, function(x) which(!!x))
Lap <- diag(rowSums(A)) - A
Z_a <- Z %*% (diag(1, 36) %x% rep(.125,8))
colnames(Z_a) <- colnames(Z)[seq(1, 288, by = 8)]
# I've come to realize they are measured in same units, so standardization is not a good idea
Z_s <- scale(Z_a, T, F)
Y_s <- scale(Y, T, F)
zmeans <- attr(Z_s, "scaled:center")
ymeans <- attr(Y_s, "scaled:center")
# zsds <- attr(Z_s, "scaled:scale")
# ysds <- attr(Y_s, "scaled:scale")

data <- list(Z_s = Z_s, Y_s = Y_s, zmeans = zmeans, ymeans = ymeans, 
             small_inds = small_inds, mis_adjs = mis_adjs, mis_inds = mis_inds, 
             distMat = distMat, Lap = Lap, 
             ids = ids, id_list = id_list, labels = labels)
saveRDS(data, file = "../data/fit_data.Rds")
