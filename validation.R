# A function to validate test set prediction 
# by using the mean thicknesses as predictor for another logistic regression
validate_roc <- function(preds_train, labels, train, test, X, Y) {
  if (!isNamespaceLoaded("pROC")) require(pROC)
  colnames(preds_train) <- paste0("V", 1:dim(preds_train)[2]) # To avoid formula errors
  df_train <- cbind.data.frame(labels = labels[train], preds_train[train,])
  classifier <- glm(labels ~ . - 1, data = df_train, family = binomial)
  test_probs <- predict(classifier, preds_train[test,], type = "response")
  entire_roc <- pROC::roc(labels[test] ~ test_probs)
  partial_roc <- pROC::roc(labels[test] ~ test_probs, 
                           partial.auc = c(.85, 1), partial.auc.correct = FALSE, 
                           partial.auc.focus = "sp")
  partial_auc <- pROC::auc(partial_roc,
                           partial.auc = c(.85, 1), partial.auc.correct = FALSE, 
                           partial.auc.focus = "sp")
  return( list(roc = entire_roc, pauc = partial_auc / .15) ) # Return normalized pAUC
}
