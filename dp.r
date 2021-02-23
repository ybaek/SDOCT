library(rstan)
m <- stan_model(file = 'dp.stan')
prior_samples <- sampling(
    m,
    data = list(K = 50, N = 100),
    pars = c("weights", "mu", "sigma", "y", "z"),
    algorithm = "Fixed_param",
    chains = 1,
    iter = 100,
    warmup = 0,
)
sample_draws <- as.matrix(prior_samples)
names_y <- grep("y", colnames(sample_draws))
names_z <- grep("z", colnames(sample_draws))
densities <- apply(sample_draws[, names_y], 1, function(row) density(row, bw = .3))
xlim <- c(0.0, 20.0)
ylim <- sapply(densities, function(e) c(min(e$y), max(e$y)))
ylim <- c(min(ylim[1,]), max(ylim[2,]))
plot(densities[[1]],
     main = "Density bw = .3",
     xlim = xlim, ylim = ylim,
     col = scales::alpha("black", .2))
for (i in 2:100) lines(densities[[i]], col = scales::alpha("black", .2))
plot(table(sample_draws[, names_z]))
