# Inference script
library(caTools)
library(coda)
pars <- readRDS("data/samples.rds")
pars_indices <- readRDS("data/indices.rds")

# Diagnosis
# Trace plots reveal error precision posterior seems to be bimodal
# (Is that a problem? we will see..)
# (The intercept and the missing values are the only ones w/ ESS < 1000)
ess <- coda::effectiveSize(t(pars))
acf(t(pars[c(1:3, 64 + 1:3), ])) # cross-correlation seems to be the problem

# FIXME: I assume implicitly some objects from main_lmm.r
# prob separate code needed
mu_hat <- rowMeans(pars[pars_indices[[1]], ])
Theta_hat <- matrix(rowMeans(pars[pars_indices[[2]], ]), R)
Siginv_hat <- matrix(0, Qstar, Qstar)
Siginv_hat[lower.tri(Siginv_hat, diag = TRUE)] <- rowMeans(pars[pars_indices[[3]], ])
Siginv_hat[upper.tri(Siginv_hat)] <- t(Siginv_hat)[upper.tri(Siginv_hat)]

gamma_hat <- mean(pars[pars_indices[[5]], ])
mask_hat <- exp(-.5 * gamma_hat^-2 * t(K2)^2)
f_hat <- U_hat %*% Theta_hat %*% mask_hat
beta_hat <- solve(crossprod(z_train), crossprod(z_train, f_hat)) # inverse map to linear scale
f_proj <- z_train %*% beta_hat # best linear approximation to f_hat
#
# Defining the held-out test data
z_test <- z[(N + 1):nrow(z), ]
y_test <- z[(N + 1):nrow(y), ]
f_proj_test <- z[(N + 1):nrow(z), ] %*% beta_hat
## (Labels were NOT provided in the actual regression model)
labels_train <- labels[1:N]
labels_test <- labels[(N + 1):nrow(z)]

# Compare: standard deviation map (threshold explicitly chosen to yield the BEST for SDM)
## (Initialize missing values of response)
mis_inds <- which(is.na(y_test))
for (i in 1:nrow(mis_inds)) {
    r <- mis_inds[i, 1]
    j <- mis_inds[i, 2]
    neighbors <- c(j-9, j-8, j-7, j-1, j+1, j+7, j+8, j+9)
    neighbors <- neighbors[neighbors > 0 & neighbors < Q]
    nmeans <- mean(y_train[r, neighbors], na.rm = T)
    y_test[r, j] <- ifelse(is.nan(nmeans), mean(y_test[r, ], na.rm = T), nmeans)
}
model_pred <- predict(
    glm(labels_test ~ prcomp(f_proj_test, scale. = T)$x, family = "binomial"),
    type = "response"
)
#
png("images/roc_kr.png")
auc <- caTools::colAUC(cbind("Projected" = model_pred), labels_test, plotROC = TRUE)
## side effect: plotting
legend("bottomleft", legend = paste0("AUC=", round(auc, 2)), bty = "n")
dev.off()
#
# Revising our goals:
# Using the actual data we can obtain almost perfect prediction
# Technically, though, this is cheating.
# We can compare for now including different principal components
# (Instability of raw data causes numerical instabilities in GLM fit)
# I concluded technically this is cheating and uninteresting
# Must report to Sam
raw_preds <- matrix(0, nrow(z_test), 5)
for (i in 1:5) raw_preds[, i] <- predict(
    glm(labels_test ~ prcomp(z_test, scale. = T)$x[, 1:(i * 10)], family = "binomial"),
    type = "response"
)
colnames(raw_preds) <- paste("PC group", 1:5)
#
png("images/roc_vs_raw.png")
auc <- caTools::colAUC(
    cbind(
        "Projected" = model_pred, 
        raw_preds
    ),
    labels_test, 
    plotROC = TRUE)
dev.off()