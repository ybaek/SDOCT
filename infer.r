library(coda)
library(pROC)
source("main_processing.r")
fit <- readRDS("data/fit_results.rds")
demo_df <- readRDS("data/demo.rds")
# Preprocess race response
demo_df$race_nih[demo_df$race_nih %in%
    c("Chinese", "Eastern", "Indian", "Japanese")] <- "Asian"
demo_df$race_nih[demo_df$race_nih %in%
    c("American Indian or Alaska Native",
      "Declined", "Not Reported", "Unknown")] <- "Others"
# Extract coefficient posterior means
burn_in <- 2000
intercept_pm <- colMeans(fit$pars[-seq(burn_in), 1:64])
theta_pm <- matrix(colMeans(fit$pars[-seq(burn_in), 65:128]), 4)
sigma_pm <- mean(fit$pars[-seq(burn_in), 319])^-.5
gamma_pm <- mean(fit$pars[-seq(burn_in), 320])
convMat  <- exp(-distMat[, knots]^2 / gamma_pm)
# Projection step: to obtain estimated deviations
K_train <- utils$form_kernel_matrix(z_train, 2)
u_train <- svd(K_train)$u[, 1:4]
f_hat   <- u_train %*% theta_pm
f_hat_full <- f_hat %*% t(convMat)
K_test  <- utils$form_kernel_matrix(z_test, 2)
u_test  <- svd(K_test)$u[, 1:4]
f_hat_test  <- u_test %*% theta_pm
f_test_full <- f_hat_test %*% t(convMat)
# Comparison: Using raw data necessarily requires imputation
mis_inds <- which(is.na(y_test), arr.ind = TRUE)
y_test2  <- y_test
y_test2[mis_inds] <- rnorm(nrow(mis_inds), mean(y_test, na.rm = T), sd(y_test, na.rm = T))
# Glaucoma labels: NOT provided at first stage
labels_train <- labels[jj_train]
labels_test  <- labels[jj_test]
#
## Second stage "model" wants to actually predict glaucoma labels
## To prevent data leakage, training data should be "used twice"
# Case 1: using raw data summary statistics / training logistic model
raw_sdm <- apply(y_test2, 1, function(x) mean(x < 0, na.rm = T))
rawdes <- cbind(labels_train, demo_df[jj_train, ], y_train)
rawdespred <- cbind(labels_test, demo_df[jj_test, ], y_test2)
rawglm <- glm(labels_train ~ . - race_primary, data = rawdes, family = "binomial")
raw_preds <- predict(rawglm, rawdespred, type = "response")
# Case 2: using denoised deviations' summary statistics / training logistic model
model_sdm <- apply(f_hat_test, 1, function(x) mean(x < 0.))
des <- cbind(labels_train, demo_df[jj_train, ], f_hat)
despred <- cbind(demo_df[jj_test, ], f_hat_test)
glmfit <- glm(labels_train ~ . - race_primary, data = des, family = "binomial")
model_preds <- predict(glmfit, despred, type = "response")
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