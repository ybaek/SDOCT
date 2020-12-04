source("validation.R")

# Testing the function with an unregularized linear fit
data <- readRDS("./data/fit_data.Rds")
lss <- with(data, lm(Y_s[train,] ~ Z_s[train,]-1))

png("./images/coefficient.png", width = 6, height = 6, units = "in", res = 120)
fields::image.plot(t(coef(lss)), col = hcl.colors(12, "RdBu", rev=T), main = "Least squares estimator")
dev.off()

preds_lss <- as.data.frame(data$Z_s %*% coef(lss))
lss_rocs <- validate_roc(preds_lss, data$labels, data$train, data$test, data$Z_s, data$Y_s)

png("./images/roc.png", width = 6, height = 6, units = "in", res = 120)
plot(lss_rocs$roc, main = "ROC curve for least squares fit")
legend("bottomright",
       legend = c(paste0("AUC=", round(lss_rocs$roc$auc, 2)*100, "%"), 
                    paste0("pAUC=", round(lss_rocs$pauc, 2)*100, "%")),
       bty = "n")
dev.off()

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
fields::image.plot(t(coef_b), col = hcl.colors(12, "RdBu", rev=T))

png("./images/c_densities.png", width = 6, height = 6, units = "in", res = 120)
plot(density(bayesfit[(2304+1+1+1),]^.5, from = 0), 
     col = scales::alpha(1, .2), xlim = c(0, 2), ylim = c(0, 5),
     main = "Posterior KDEs of c")
for (i in 2:1712) {
    lines(density(bayesfit[(2304+1+1+i),]^.5, from = 0), col = scales::alpha(1, .1))
}
dev.off()

# Standardization can lead to uninterpretable model fit
bayesfit_s <- readRDS("./objects/samples_lmm2_cpp.Rds")
png("./images/c_densities_s.png", width = 6, height = 6, units = "in", res = 120)
plot(density(bayesfit_s[(2304+1+1+1),]^.5, from = 0), 
     col = scales::alpha(1, .2), xlim = c(0, 3), ylim = c(0, 3),
     main = "Prior decision incoherent with standardization")
for (i in 2:1712) {
    lines(density(bayesfit_s[(2304+1+1+i),]^.5, from = 0), col = scales::alpha(1, .1))
}
dev.off()

# preds_b <- as.data.frame(data$Z_s %*% coef_b)
# bayes_rocs <- validate_roc(preds_b, data$labels, data$train, data$test, data$Z_s, data$Y_s)
# lines(bayes_rocs$roc, col = 2) # You can't tell the difference

