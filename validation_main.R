source("validation.R")

# Testing the function with an unregularized linear fit
data <- readRDS("./data/fit_data.Rds")
lss <- with(data, lm(Y_s[train,] ~ Z_s[train,]-1))
# fields::image.plot(t(coef(lss)), col = hcl.colors(12, "RdBu", rev=T))
preds_lss <- as.data.frame(data$Z_s %*% coef(lss))
lss_rocs <- validate_roc(preds_lss, labels, train, test, data$Z_s, data$Y_s)
plot(lss_rocs$roc, main = "ROC curve for least squares fit")

# Same model but standardized makes no difference
# Z_ss <- scale(data$Z_s, F, T)
# Y_ss <- scale(data$Y_s, F, T)
# lss_s <- lm(Y_ss[train,] ~ Z_ss[train,]-1)
# preds_lss_s <- as.data.frame(data$Z_s %*% coef(lss_s))
# # fields::image.plot(t(coef(lss_s)), col = hcl.colors(12, "RdBu", rev=T))
# lss_s_rocs <- validate_roc(preds_lss_s, labels, train, test, Z_ss, Y_ss)
# lines(lss_s_rocs$roc, col = 2)

bayesfit <- readRDS("./objects/samples_lmm1_cpp.Rds")
coef_b <- matrix(rowMeans(bayesfit[1:2304,]), 36)
preds_b <- as.data.frame(data$Z_s %*% coef_b)
bayes_rocs <- validate_roc(preds_b, labels, train, test, data$Z_s, data$Y_s)
lines(bayes_rocs$roc, col = 2) # You can't tell the difference
