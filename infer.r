# Inference script
library(coda)
library(loo)
library(pROC)
source("main_processing.r")
pars <- readRDS("data/samples2.rds")
pars_indices <- readRDS("data/indices2.rds")
lps <- readRDS("data/lps2.rds")
demo_df <- readRDS("data/demo.rds")
# Preprocess race response
demo_df$race_nih[demo_df$race_nih %in%
    c("Chinese", "Eastern", "Indian", "Japanese")] <- "Asian"
demo_df$race_nih[demo_df$race_nih %in%
    c("American Indian or Alaska Native",
      "Declined", "Not Reported", "Unknown")] <- "Others"

# Diagnosis
# Trace plots reveal error precision posterior seems to be bimodal
# (Is that a problem? we will see..)
# (The intercept and the missing values are the only ones w/ ESS < 1000)
# ess <- coda::effectiveSize(t(pars))
# ess_lp <- coda::effectiveSize(t(lps))
# acf(t(pars[c(1:3, 64 + 1:3), ])) # cross-correlation seems to be the problem

# Which model should we choose for a bandwidth parameter?
# WAIC / LOOCV statistic
# loo::loo(lps)$elpd_loo / 1e+6 # log-density scale
# loo::waic(lps)$waic / 1e+7 # deviance scale

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
#
f_proj_full <- f_proj %*% mask_hat
f_proj_test_full <- f_proj_test %*% mask_hat

## (Labels were NOT provided in the actual regression model)
labels_train <- labels[jj_train]
labels_test  <- labels[jj_test]

# Compare: standard deviation map
# How do we choose where to threshold??
mis_inds <- which(is.na(y_test), arr.ind = TRUE)
y_test2  <- y_test
y_test2[mis_inds] <- rnorm(nrow(mis_inds), mean(y_test, na.rm = T), sd(y_test, na.rm = T))
# Threshold parameter is found to be optimizing WITHIN the training set
raw_sdm <- apply(y_test, 1, function(x) mean(x < 0, na.rm = T))
rawdes <- cbind(labels_train, demo_df[jj_train, ], y_train)
rawdespred <- cbind(labels_test, demo_df[jj_test, ], y_test)
rawglm <- glm(labels_train ~ . - race_primary, data = rawdes, family = "binomial")
raw_preds <- predict(rawglm, rawdespred, type = "response")
# PROBLEM: latent factors are not directly interpretable on thickness scale
# the scale in particular is not directly comparable to SDM.
model_sdm1 <- apply(f_proj_test, 1, function(x) mean(x < 0.))
model_sdm2 <- apply(f_proj_test_full, 1, function(x) mean(x < 0.))
## Formatting for logistic regression
des1 <- cbind(labels_train, demo_df[jj_train, ], f_proj)
des2 <- cbind(labels_train, demo_df[jj_train, ], f_proj_full)
despred1 <- cbind(demo_df[jj_test, ], f_proj_test)
despred2 <- cbind(demo_df[jj_test, ], f_proj_test_full)
glmfit1 <- glm(labels_train ~ . - race_primary, data = des1, family = "binomial")
glmfit2 <- glm(labels_train ~ . - race_primary, data = des2, family = "binomial")
model_preds1 <- predict(glmfit1, despred1, type = "response")
model_preds2 <- predict(glmfit2, despred2, type = "response")
#
roc_preds1 <- pROC::roc(labels_test, c(model_preds1))
preds_ci1 <- ci.se(roc_preds1, specifities = seq(0, 1, .01))
roc_preds11 <- pROC::roc(labels_test, c(model_preds2))
preds_ci11 <- ci.se(roc_preds11, specifities = seq(0, 1, .01))
roc_preds2 <- pROC::roc(labels_test, c(raw_preds))
preds_ci2 <- ci.se(roc_preds2, specifities = seq(0, 1, .01))
roc_sdm1 <- pROC::roc(labels_test, c(model_sdm1))
sdm_ci1 <- ci.se(roc_sdm1, specifities = seq(0, .1, .01))
roc_sdm11 <- pROC::roc(labels_test, c(model_sdm2))
sdm_ci11 <- ci.se(roc_sdm11, specifities = seq(0, .1, .01))
roc_sdm2 <- pROC::roc(labels_test, c(raw_sdm))
sdm_ci2 <- ci.se(roc_sdm2, specifities = seq(0, .1, .01))
#
png("images/roc_logistic.png")
plot(roc_preds1, col = 4)
# plot(roc_preds11, col = 3, add = T)
plot(roc_preds2, add = T, col = 2)
plot(preds_ci1, type = "shape", col = scales::alpha(4, .2), no.roc = TRUE)
plot(preds_ci2, type = "shape", col = scales::alpha(2, .2), no.roc = TRUE)
abline(v = .85, lty = 2)
legend(
    "bottomright",
    legend = c("Raw", "Model"),
    col    = c(2, 4),
    lty    = 1
)
dev.off()
#
png("images/roc_sdm.png")
plot(roc_sdm1, col = 4)
plot(roc_sdm2, add = T, col = 2)
plot(sdm_ci1, type = "shape", col = scales::alpha(4, .2), no.roc = TRUE)
plot(sdm_ci2, type = "shape", col = scales::alpha(2, .2), no.roc = TRUE)
abline(v = .85, lty = 2)
legend(
    "bottomright",
    legend = c("Raw", "Model"),
    col    = c(2, 4),
    lty    = 1
)
dev.off()

ci.auc(auc(roc_preds1), method = "boot")
ci.auc(auc(roc_preds2), method = "boot")
ci.auc(auc(roc_sdm1), method = "boot")
ci.auc(auc(roc_sdm2), method = "boot")
ci.auc(auc(roc_preds1, partial.auc = c(.85, 1), partial.auc.correct = FALSE, partial.auc.focus = "sp")) / .15
ci.auc(auc(roc_preds2, partial.auc = c(.85, 1), partial.auc.correct = FALSE, partial.auc.focus = "sp")) / .15
ci.auc(auc(roc_sdm1, partial.auc = c(.85, 1), partial.auc.correct = FALSE, partial.auc.focus = "sp")) / .15
ci.auc(auc(roc_sdm2, partial.auc = c(.85, 1), partial.auc.correct = FALSE, partial.auc.focus = "sp")) / .15

#
# Generalized LASSO post-processing
# Mostly to reveal tricky inferential issues
library(genlasso)
full_fit <- fusedlasso1d(y = f_proj_full[, 35], X = z_train, gamma = 1.0)
throw_fit1 <- fusedlasso1d(y = f_proj_full[, 35], X = z_train[, 1:96], gamma = 1.0)
throw_fit2 <- fusedlasso1d(y = f_proj_full[, 35], X = z_train[, 1:192], gamma = 1.0)
png("images/postproc.png")
plot(full_fit$beta[, 1600], type = "l",
    xlab = "p", ylab = expression(beta[35]), col = scales::alpha(1, .5))
lines(full_fit$beta[, 1500], col = scales::alpha(1, .5))
lines(full_fit$beta[, 1300], col = scales::alpha(1, .5))
lines(full_fit$beta[, 1200], col = scales::alpha(1, .5))
abline(v = c(96, 192), col = scales::alpha(2, .9), lwd = 2)
text(x = c(48, 144, 240), y = .65, labels = c("I", "T", "S"))
dev.off()