# A function to validate test set prediction 
# by using the mean thicknesses as predictor for another logistic regression
validate_roc <- function(coef_train, labels, train, test, X, Y) {
  if (!isNamespaceLoaded("pROC")) require(pROC)
  predicted_train <- as.data.frame(X %*% coef_train)
  colnames(predicted_train) <- paste0("V", 1:dim(predicted_train)[2]) # To avoid formula errors
  df_train <- cbind.data.frame(labels = labels[train], predicted_train[train,])
  classifier <- glm(labels ~ . - 1, data = df_train, family = binomial)
  test_probs <- predict(classifier, predicted_train[test,], type = "response")
  entire_roc <- pROC::roc(labels[test] ~ test_probs)
  partial_roc <- pROC::roc(labels[test] ~ test_probs, 
                           partial.auc = c(.85, 1), partial.auc.correct = FALSE, 
                           partial.auc.focus = "sp")
  partial_auc <- pROC::auc(partial_roc,
                           partial.auc = c(.85, 1), partial.auc.correct = FALSE, 
                           partial.auc.focus = "sp")
  return( list(roc = entire_roc, pauc = partial_auc / .15) ) # Return normalized pAUC
}

# Testing the function with an unregularized linear fit
data <- readRDS("./data/fit_data.Rds")
attach(data)
N <- dim(Z_s)[1]; P <- dim(Z_s)[2]; Q <- dim(Y_s)[2]; J <- max(ids)
set.seed(2020-11-25)
train <- sort(sample(N, round(N*.8)))
test <- setdiff(1:N, train)
lss <- lm(Y_s[train,] ~ Z_s[train,]-1)
# fields::image.plot(t(coef(lss)), col = hcl.colors(12, "RdBu", rev=T))
lss_rocs <- validate_roc(coef(lss), labels, train, test, Z_s, Y_s)
plot(lss_rocs$roc, main = "ROC curve for least squares fit")

# Same model but standardized
Z_ss <- scale(Z_s, F, T)
Y_ss <- scale(Y_s, F, T)
lss_s <- lm(Y_ss[train,] ~ Z_ss[train,]-1)
# fields::image.plot(t(coef(lss_s)), col = hcl.colors(12, "RdBu", rev=T))
lss_s_rocs <- validate_roc(coef(lss_s), labels, train, test, Z_ss, Y_ss)
lines(lss_s_rocs$roc, col = 2)

legend("bottomright",
        legend = c("Raw scale (AUC .84 / pAUC .65)",
                   "Standardized (AUC .86 / pAUC .63)"),
        col = 1:2, lty = 1, lwd = 2, cex = .8)