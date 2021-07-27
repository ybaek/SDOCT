library(coda)
library(pROC)
library(scales)
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
f_hat_full <- sweep(f_hat %*% t(convMat), 2, intercept_pm, "+")
K_test  <- utils$form_kernel_matrix(z_test, 2)
u_test  <- svd(K_test)$u[, 1:4]
f_hat_test  <- u_test %*% theta_pm
f_test_full <- sweep(f_hat_test %*% t(convMat), 2, intercept_pm, "+")
# (Using raw data necessarily requires imputation -- make a copy)
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
# raw_sdm <- apply(y_test, 1, function(x) mean(x < .02, na.rm = T))
rawdes <- cbind(labels_train, demo_df[jj_train, ], y_train)
rawdespred <- cbind(labels_test, demo_df[jj_test, ], y_test2)
rawglm <- glm(labels_train ~ . - race_primary, data = rawdes, family = "binomial")
raw_preds <- predict(rawglm, rawdespred, type = "response")
# Case 2: using denoised deviations' summary statistics / training logistic model
# model_sdm <- apply(f_test_full, 1, function(x) mean(x < 0.))
des <- cbind(labels_train, demo_df[jj_train, ], f_hat_full)
despred <- cbind(demo_df[jj_test, ], f_test_full)
glmfit <- glm(labels_train ~ . - race_primary, data = des, family = "binomial")
model_preds <- predict(glmfit, despred, type = "response")
#
roc_model <- pROC::roc(labels_test, c(model_preds))
ci_model  <- ci.se(roc_preds1, specifities = seq(0, 1, .01))
roc_raw <- pROC::roc(labels_test, c(raw_preds))
ci_raw  <- ci.se(roc_preds2, specifities = seq(0, 1, .01))
# roc_sdm1 <- pROC::roc(labels_test, c(model_sdm))
# sdm_ci1 <- ci.se(roc_sdm1, specifities = seq(0, .1, .01))
# roc_sdm2 <- pROC::roc(labels_test, c(raw_sdm))
# sdm_ci2 <- ci.se(roc_sdm2, specifities = seq(0, .1, .01))
#
unadj_pauc <- function(roc) {
    # Unnormalized, un-corrected partial AUC in 15% FPR range
    pROC::auc(roc, partial.auc = c(.85, 1),
        partial.auto.correct = F, partial.auc.focus = "sensitivity")
}
auc_tab <- rbind(
    c(ci.auc(auc(roc_model), method = "boot"),
      ci.auc(auc(roc_raw), method = "boot"))
)
pauc_tab <- rbind(
    c(ci.auc(unadj_pauc(roc_model)) / (1 - .85),
      ci.auc(unadj_pauc(roc_raw)) / (1 - .85))
)
roc.test(auc(roc_model), auc(roc_raw)) # p == .028
roc.test(unadj_pauc(roc_model), unadj_pauc(roc_raw))
#
png("images/roc_logistic.png")
plot(roc_model, col = 4)
plot(roc_raw, add = T, col = 2)
plot(ci_model, type = "shape", col = scales::alpha(4, .2), no.roc = TRUE)
plot(ci_raw, type = "shape", col = scales::alpha(2, .2), no.roc = TRUE)
abline(v = .85, lty = 2)
legend(
    "bottomright",
    legend = c("Raw", "Model"),
    col    = c(2, 4),
    lty    = 1
)
dev.off()