# Individual smoothing of the images
data <- readRDS("data/fit_data.Rds")
source("adj_grid.R")
# Graph Laplacian
Q <- diag(rowSums(adj_grid(8))) - adj_grid(8)
lambdas <- eigen(Q)$values
smoothed_Ys <- t(
    apply(
        data$Y_s[data$test, ],
        1,
        function(row) (1.0 - rho) *
            solve(rho * Q + diag(1.0 - rho, 64)) %*% row)
)
