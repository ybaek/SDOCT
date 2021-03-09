# Inference script
library(coda)
library(ROCR)
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
gamma_hat <- mean(pars[pars_indices[[5]], ])
f_hat <- U_hat %*% Theta_hat %*% exp(-.5 * gamma_hat^-2 * t(K2)^2)
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
sdm_pred <- apply(y_test, 1, function(x) mean(x < .02, na.rm = T))
model_pred <- predict(glm(labels_test ~ f_proj_test, family = "binomial"), type = "response")
#
model_pred_obj1 <- prediction(sdm_pred, labels_test)
model_pred_obj2 <- prediction(model_pred, labels_test)
model_perf_obj1 <- performance(model_pred_obj1, "tpr", "fpr")
model_perf_obj2 <- performance(model_pred_obj2, "tpr", "fpr")
#
plot(model_perf_obj1)
plot(model_perf_obj2, add = T, col = 2)
