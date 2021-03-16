# Inference script
library(coda)
library(pROC)
source("main_processing.r")
pars <- readRDS("data/samples1.rds")
pars_indices <- readRDS("data/indices1.rds")
lps <- readRDS("data/lps1.rds")

# Diagnosis
# Trace plots reveal error precision posterior seems to be bimodal
# (Is that a problem? we will see..)
# (The intercept and the missing values are the only ones w/ ESS < 1000)
ess <- coda::effectiveSize(t(pars))
ess_lp <- coda::effectiveSize(t(lps))
acf(t(pars[c(1:3, 64 + 1:3), ])) # cross-correlation seems to be the problem

# Posterior summaries
mu_hat <- rowMeans(pars[pars_indices[[1]], ])
Theta_hat <- matrix(rowMeans(pars[pars_indices[[2]], ]), R)
Siginv_hat <- matrix(0, Qstar, Qstar)
Siginv_hat[lower.tri(Siginv_hat, diag = TRUE)] <- rowMeans(pars[pars_indices[[3]], ])
Siginv_hat[upper.tri(Siginv_hat)] <- t(Siginv_hat)[upper.tri(Siginv_hat)]

gamma_hat <- mean(pars[pars_indices[[5]], ])
mask_hat <- exp(-.5 * gamma_hat^-2 * t(K2)^2)
f_hat <- U_hat %*% Theta_hat
beta_hat <- solve(crossprod(z_train), crossprod(z_train, f_hat)) # inverse map to linear scale
f_proj <- z_train %*% beta_hat # best linear approximation to f_hat
f_proj_test <- z_test %*% beta_hat

## (Labels were NOT provided in the actual regression model)
labels_train <- labels[jj_train]
labels_test  <- labels[jj_test]

# Compare: standard deviation map (threshold explicitly chosen to yield the BEST for SDM)
mis_inds <- which(is.na(y_test), arr.ind = TRUE)
y_test2  <- y_test
y_test2[mis_inds] <- rnorm(nrow(mis_inds), mean(y_test, na.rm = T), sd(y_test, na.rm = T))
# Threshold parameter is found to be optimizing WITHIN the training set
raw_sdm <- apply(y_test, 1, function(x) mean(x < .023, na.rm = T))
raw_coef <- coef(glm(labels_train ~ y_train, family = "binomial"))
raw_preds <- (1 + exp(-raw_coef[1] - y_test2 %*% raw_coef[-1]))^-1
# PROBLEM: latent factors are not directly interpretable on thickness scale
# the scale in particular is not directly comparable to SDM.
model_sdm <- apply(f_proj_test, 1, function(x) mean(x < 21e-5))
model_coef <- coef(glm(labels_train ~ f_proj, family = "binomial"))
model_preds <- (1 + exp(-model_coef[1] - f_proj_test %*% model_coef[-1]))^-1

#
roc_preds1 <- pROC::roc(labels_test, c(model_preds))
roc_preds2 <- pROC::roc(labels_test, c(raw_preds))
roc_sdm1 <- pROC::roc(labels_test, c(model_sdm))
roc_sdm2 <- pROC::roc(labels_test, c(raw_sdm))
plot(roc_preds1)
plot(roc_preds2, add = T, col = 2)
plot(roc_sdm1)
plot(roc_sdm2, add = T, col = 2)

auc(roc_preds1, partial.auc = c(.85, 1), partial.auc.correct = TRUE, partial.auc.focus = "sp")
auc(roc_preds2, partial.auc = c(.85, 1), partial.auc.correct = TRUE, partial.auc.focus = "sp")
auc(roc_sdm1, partial.auc = c(.85, 1), partial.auc.correct = TRUE, partial.auc.focus = "sp")
auc(roc_sdm2, partial.auc = c(.85, 1), partial.auc.correct = TRUE, partial.auc.focus = "sp")
