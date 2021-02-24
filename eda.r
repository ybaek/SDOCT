plots <- new.env()
source('plots.r', local = plots)

# Smoothing the quantiles (problem: measurements are integers)
# Simple moving average
center <- readRDS("data/stats.rds")
ds <- readRDS("data/dataset.rds")
series <- center[[1]]["q5", ]
ma <- numeric(0)
for (i in seq(length(series))) {
    subseries <- window(series, start = max(i - 7, 1), end = min(i + 7, 768))
    ma <- c(ma, mean(subseries))
}

library(ggplot2)
library(reshape2)
plot(ds$Z[1,], ylim = c(0, max(ds$Z)), type = "l")
lines(ds$Z[1,] - ma, col = 2)
example_m <- matrix(ds$Y[1,], 8)
center_m  <- matrix(center[[2]]["q5", ], 8)
ggplot(melt(example_m)) + geom_tile(aes(x = Var1, y = Var2, fill = value)) + 
    scale_fill_distiller(palette = "Reds", lim = c(0, max(ds$Y, na.rm = T)), direction=1)
ggplot(melt(example_m - center_m)) + geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    scale_fill_distiller(palette = "Reds", lim = c(0, max(ds$Y, na.rm = T)), direction=1)

# I don't know if this analysis is worthwhile or not, but
# let's try a k-median per each patient's image
library(cluster)
Y_centered <- sweep(ds$Y, 2, c(center_m))
# Out-of-shelf clustering algorithms aren't good at handling missing values
Y_centered[is.na(Y_centered)] <- mean(Y_centered, na.rm = T)
# Try different values of k's
ks <- c(4, 8, 16, 32)
clusters <- list(4)
for (i in 1:4) {
    clusters[[i]] <- apply(Y_centered, 1, function(row) cluster::pam(row, k = ks[i])$id.med)
}
clusters <- lapply(clusters, function(m) matrix(table(m), 8))
gList <- vector(mode = "list", length = 4)
for (i in 1:4) {
    gList[[i]] <- ggplot(melt(clusters[[i]])) +
    geom_tile(aes(x = Var1, y = Var2, fill = value)) + 
    labs(x = "", y = "", fill = "Frequency",
         title = paste0("k=", ks[i])) + 
    theme_minimal()
}
png("images/medoids.png", 1080, 1080, res = 150)
plots$multiplot(plotlist = gList, cols = 2)
dev.off()


#####
# How does nonlinearity (tanh) affect the emission model?
rm(list=ls())
gc(reset=T)
N <- 1000
L <- chol(matrix(c(1, .9, .9, 1), 2))
x <- matrix(rnorm(2000), 1000) 
modes <- rbind(c(-2, 2), c(2, -2), c(1, 0), c(0, 1))
groups <- sample(4, 1000, replace = T)
y <- x + modes[groups, ]
png("images/original.png", 1080, 1080, res = 200)
plot(y, xlab = "", ylab = "", main = "Mixture of 4 BVNS")
points(modes, col = 2, cex = 2, pch = 19)
dev.off()

png("images/no_jitter.png", 1080, 1080, res = 200)
par(mfrow = c(2, 2),
    mar = c(3, 3, 2, 1), mgp = c(2, .4, 0))
plot(tanh(y), xlab = "", ylab = "", main = "alpha = 1")
plot(tanh(y/.1), xlab = "", ylab = "", main = "alpha = .1")
plot(tanh(y/2), xlab = "", ylab = "", main = "alpha = 2")
plot(tanh(y/10), xlab = "", ylab = "", main = "alpha = 10")
dev.off()

png("images/jitter.png", 1080, 1080, res = 200)
par(mfrow = c(2, 2),
    mar = c(3, 3, 2, 1), mgp = c(2, .4, 0))
plot(tanh(y) + .1 * matrix(rnorm(2000), ncol = 2) %*% L,
    xlab = "", ylab = "", main = "alpha = 1")
plot(tanh(y/.1) + .1 * matrix(rnorm(2000), ncol = 2) %*% L,
    xlab = "", ylab = "", main = "alpha = .1")
plot(tanh(y/2) + .1 * matrix(rnorm(2000), ncol = 2) %*% L,
    xlab = "", ylab = "", main = "alpha = 2")
plot(tanh(y/10) + .1 * matrix(rnorm(2000), ncol = 2) %*% L,
    xlab = "", ylab = "", main = "alpha = 10")
dev.off()