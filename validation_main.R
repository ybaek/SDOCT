library(ggplot2)
library(reshape2)
library(scales)
source("prior_distances.R")
source("adj_grid.R")

# Downsampling prior mean
prior_mean <- prior_mean[seq(1, 288, 8), ]
# Graph Laplacian and eigenvalues
Q <- diag(rowSums(adj_grid(8))) - adj_grid(8)
lambdas <- eigen(Q)$values

# Plotting the coefficient
data <- readRDS("data/fit_data.Rds")
pm <- array(0, dim = c(36 * 64, 5))
for (k in 1:5) {
    pm[, k] <- rowMeans(
                    readRDS(paste0("./objects/lmmfit", k, ".rds"))[1:2304, ]
                )
}
pmm <- matrix(rowMeans(pm), 36, 64)
g <- ggplot() +
    geom_tile(data = reshape2::melt(pmm),
              aes(as.numeric(Var1), as.numeric(Var2), fill = value))
# Generating a heatmap
g <- g + labs(x = "cpRNFL", y = "macula", fill = "Deviation") +
    scale_fill_distiller(palette = "RdBu",
        direction = -1, lim = c(-1, 1) * max(abs(pm))) +
    theme(panel.background = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())
# Add bounding box to heatmap for comparison to prior
prior_inds <- which(!!prior_mean, arr.ind = TRUE)
prior_inds_y <- unname(prior_inds[, 2])
prior_inds_x <- as.integer(rownames(prior_inds))
prior_inds_x <- ifelse(prior_inds_x <= 121,
                       prior_inds_x %/% 8 + 1,
                       (prior_inds_x - 480) %/% 8 + 1)
g <- g +
    geom_rect(
        data = NULL,
        aes(xmin = prior_inds_x - .5,
            xmax = prior_inds_x + .5,
            ymin = prior_inds_y - .5,
            ymax = prior_inds_y + .5),
        fill = "transparent", color = "grey50"
    )
ggsave("./images/coefficient_bayes.png", g, "png",
       width = 6, height = 6, units = "in", dpi = 150)

# Plotting ROC curves
# Rule 1. counting low deviations (<0?)
# Rule 2. Construct a different logistic fit and use
# deviations as predictors
validate_roc <- function(pred_dev, test_response, test_X = NULL, thres, ...) {
    if (all(loadedNamespaces() != "pROC")) require(pROC)
    sdm_Y <- apply(pred_dev, 1, 
                 function(x) sum(x < thres, na.rm = TRUE))
    sdm <- sdm_Y
    if (!is.null(test_X)) {
        sdm_Z <- apply(test_X, 1,
                 function(x) sum(x < thres, na.rm = TRUE))
        sdm <- sdm_Y + sdm_Z / (64 + 36)
    }
    full_roc <- pROC::roc(test_response ~ sdm, ...)
    part_roc <- pROC::roc(test_response ~ sdm,
                          partial.auc = c(.85, 1),
                          partial.auc.correct = FALSE,
                          partial.auc.focus = "sp", ...)
    return(list(roc = full_roc, proc = part_roc))
}
# Consider both raw and (CAR-based) smoothed deviations
smoothed_Y_s <- t(
    apply(data$Y_s,
          1, 
          function(row) .01 * solve(.99 * Q + diag(.01, 64)) %*% row)
)
lss_sdm_rocs1 <- lss_sdm_rocs2 <- vector(mode = "list", length = 5L)
raw_sdm_rocs1 <- raw_sdm_rocs2 <- vector(mode = "list", length = 5L)
smoothed_sdm_rocs1 <- smoothed_sdm_rocs2 <- vector(mode = "list", length = 5L)
for (k in 1:5) {
    preds <- data$Z_s[data$test[[k]], ] %*% matrix(pm[, k], 36, 64)
    lss_sdm_rocs1[[k]] <- validate_roc(
        preds, data$labels[data$test[[k]]],
        thres = 0L
    )
    lss_sdm_rocs2[[k]] <- validate_roc(
        preds, data$labels[data$test[[k]]],
        data$Z_s[data$test[[k]], ], 0L
    )
    raw_sdm_rocs1[[k]] <- validate_roc(
        data$Y_s[data$test[[k]], ], data$labels[data$test[[k]]], 
        thres = 0L
    )
    raw_sdm_rocs2[[k]] <- validate_roc(
        data$Y_s[data$test[[k]], ], data$labels[data$test[[k]]], 
        data$Z_s[data$test[[k]], ], 0L
    )
    smoothed_sdm_rocs1[[k]] <- validate_roc(
        smoothed_Y_s[data$test[[k]], ], data$labels[data$test[[k]]], 
        thres = 0L
    )
    smoothed_sdm_rocs2[[k]] <- validate_roc(
        smoothed_Y_s[data$test[[k]], ], data$labels[data$test[[k]]], 
        data$Z_s[data$test[[k]], ], 0L
    )
}
# Bootstrap the sensitivity (at specificity .95) + AUC / partial AUC
sens_ci1 <- sapply(lss_sdm_rocs1, function(x) ci.se(x$roc, .95))
sens_ci2 <- sapply(raw_sdm_rocs1, function(x) ci.se(x$roc, .95))
auc_ci1  <- sapply(lss_sdm_rocs1, function(x) ci.auc(x$roc, .95))
auc_ci2  <- sapply(raw_sdm_rocs1, function(x) ci.auc(x$roc, .95))
pauc_ci1 <- sapply(lss_sdm_rocs1, function(x) ci.auc(x$proc, .95))
pauc_ci2 <- sapply(raw_sdm_rocs1, function(x) ci.auc(x$proc, .95))

# Reduce the list components down to their
# (1 - Specifities), Sensitivities, Thresholds, and AUCs
lss_sdm_rocs1 <- lapply(lss_sdm_rocs1,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
lss_sdm_rocs2 <- lapply(lss_sdm_rocs2,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
raw_sdm_rocs1 <- lapply(raw_sdm_rocs1,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
raw_sdm_rocs2 <- lapply(raw_sdm_rocs2,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
smoothed_sdm_rocs1 <- lapply(smoothed_sdm_rocs1,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
smoothed_sdm_rocs2 <- lapply(smoothed_sdm_rocs2,
    function(x) list(fpr = 1 - x$roc$specificities,
                     sens = x$roc$sensitivities,
                     thres = x$roc$thresholds,
                     auc = x$roc$auc,
                     pauc = x$proc$auc / .15))
# Since ROCs are supported on different number of points
# will make simple linear interpolations
structure_roc_for_plot <- function(l) {
    # Assume object structured from above!
    grid <- sort(Reduce(union, lapply(l, function(x) x$thres)))
    fpr_interp <- sapply(l, function(x) approx(x$thres, x$fpr, grid)$y)
    tpr_interp <- sapply(l, function(x) approx(x$thres, x$sens, grid)$y)
    avg_fpr <- rowMeans(fpr_interp, na.rm = T)
    avg_tpr <- rowMeans(tpr_interp, na.rm = T)
    se_tpr <- apply(tpr_interp, 1, function(x) sd(x, T))
    return(
        list(
            avg_fpr = avg_fpr,
            avg_tpr = avg_tpr,
            avg_lower = avg_tpr - se_tpr,
            avg_upper = avg_tpr + se_tpr
        )
    )
}
lss_for_plot1 <- structure_roc_for_plot(lss_sdm_rocs1)
raw_for_plot1 <- structure_roc_for_plot(raw_sdm_rocs1)
smoothed_for_plot1 <- structure_roc_for_plot(smoothed_sdm_rocs1)
lss_for_plot2 <- structure_roc_for_plot(lss_sdm_rocs2)
raw_for_plot2 <- structure_roc_for_plot(raw_sdm_rocs2)
smoothed_for_plot2 <- structure_roc_for_plot(smoothed_sdm_rocs2)

# Showing ROCs for two different scenarios (Y only / Y + Z)
jpeg("./images/roc1.jpeg", width = 4.5, height = 3.25, units = "in", res = 200)
.env$nice_par()
plot(lss_for_plot1$avg_fpr,
     lss_for_plot1$avg_tpr,
     type = "l", col = 2, lwd = 2,
     xlab = "1-Specificity", ylab = "Sensitivity")
     # main = "Using macula deviation maps")
lines(raw_for_plot1$avg_fpr, raw_for_plot1$avg_tpr, lwd = 2)
# lines(smoothed_for_plot1$avg_fpr,
#      smoothed_for_plot1$avg_tpr, col = "darkgreen", lwd = 2)
polygon(c(lss_for_plot1$avg_fpr, rev(lss_for_plot1$avg_fpr)),
        c(lss_for_plot1$avg_upper, rev(lss_for_plot1$avg_lower)),
        border = scales::alpha(2, .2),
        col = scales::alpha(2, .2))
polygon(c(raw_for_plot1$avg_fpr, rev(raw_for_plot1$avg_fpr)),
        c(raw_for_plot1$avg_upper, rev(raw_for_plot1$avg_lower)),
        border = scales::alpha(1, .2),
        col = scales::alpha(1, .2))
# polygon(c(smoothed_for_plot1$avg_fpr[-2], 
#         rev(smoothed_for_plot1$avg_fpr[-2])),
#         c(smoothed_for_plot1$avg_upper[-2], 
#         rev(smoothed_for_plot1$avg_lower[-2])),
#         border = scales::alpha("darkgreen", .2),
#         col = scales::alpha("darkgreen", .2))
abline(a = 0, b = 1, lty = 3)
abline(v = .15, lty = 4)
legend("bottomright",
       legend = c("Raw (79/39)", 
#                  "Smoothed (72/34)", 
                  "Bayesian (78/47)"),
       col = c(1, 
       #"darkgreen", 
       2), 
       lty = 1, cex = .8)
dev.off()

png("./images/roc2.png", width = 6, height = 6, units = "in", res = 120)
plot(lss_for_plot2$avg_fpr,
     lss_for_plot2$avg_tpr,
     type = "l", col = 2, lwd = 2,
     xlab = "1-Specificity", ylab = "Sensitivity",
     main = "Using macula and ONH deviation maps")
lines(raw_for_plot2$avg_fpr, raw_for_plot2$avg_tpr, lwd = 2)
lines(smoothed_for_plot2$avg_fpr, 
      smoothed_for_plot2$avg_tpr, col = "darkgreen", lwd = 2)
polygon(c(lss_for_plot2$avg_fpr, rev(lss_for_plot2$avg_fpr)),
        c(lss_for_plot2$avg_upper, rev(lss_for_plot2$avg_lower)),
        border = scales::alpha(2, .2),
        col = scales::alpha(2, .2))
polygon(c(raw_for_plot2$avg_fpr, rev(raw_for_plot2$avg_fpr)),
        c(raw_for_plot2$avg_upper, rev(raw_for_plot2$avg_lower)),
        border = scales::alpha(1, .2),
        col = scales::alpha(1, .2))
polygon(c(smoothed_for_plot2$avg_fpr, 
        rev(smoothed_for_plot2$avg_fpr)),
        c(smoothed_for_plot2$avg_upper, 
        rev(smoothed_for_plot2$avg_lower)),
        border = scales::alpha("darkgreen", .2),
        col = scales::alpha("darkgreen", .2))
abline(a = 0, b = 1, lty = 3)
abline(v = .15, lty = 4)
legend("bottomright",
       legend = c("Raw (79/40)",
                  "Smoothed (77/45)",
                  "Smoothed + ONH (78/46)"),
       col = c(1, "darkgreen", 2), 
       lty = 1)
dev.off()

# Bootstrap cross-validation statistics
# (AUC + pAUC)
lss_auc_samples  <- sapply(lss_sdm_rocs1, function(x) x$auc)
lss_pauc_samples <- sapply(lss_sdm_rocs1, function(x) x$pauc)
raw_auc_samples  <- sapply(raw_sdm_rocs1, function(x) x$auc)
raw_pauc_samples <- sapply(raw_sdm_rocs1, function(x) x$pauc)

boots <- matrix(0, 1000, 4)
for (i in 1:1000) {
    boots[i, 1] <- mean(lss_auc_samples[sample(5, 5, TRUE)])
    boots[i, 2] <- mean(lss_pauc_samples[sample(5, 5, TRUE)])
    boots[i, 3] <- mean(raw_auc_samples[sample(5, 5, TRUE)])
    boots[i, 4] <- mean(raw_pauc_samples[sample(5, 5, TRUE)])
}


# Inference: on other possible model parameters
# Why we can diagnose ill model fit
# c_samples <- samples[(2304 + 1 + 1):(2304 + 1 + 1712), ]^.5
# c_densities <- apply(c_samples, 1, function(x) density(x, from = 0))
# plot_xmax <- max(sapply(c_densities, function(y) max(y$x)))
# plot_ymax <- max(sapply(c_densities, function(x) max(x$y)))
# png("./images/densities.png", width = 6, height = 6, units = "in", res = 120)
# plot(c_densities[[1]],
#      xlab = "c",
#      ylab = "density",
#      main = "",
#      xlim = c(0, plot_xmax),
#      ylim = c(0, plot_ymax),
#      col = scales::alpha(1, .2))
# for (i in 2:1712) {
#     lines(c_densities[[i]], col = scales::alpha(1, .2))
# }
# dev.off()
# 
# psiMean <- matrix(0, 36, 36)
# psiMean[lower.tri(psiMean, diag = T)] <- rowMeans(
#     samples[(2304 + 1 + 1712 + 1):(2304 + 1 + 1712 + 666), ]
# )
# psiMean[upper.tri(psiMean)] <- t(psiMean)[upper.tri(psiMean)]
# 
# gPsi <- ggplot(reshape2::melt(solve(psiMean))) +
#     geom_tile(aes(Var1, Var2, fill = value)) +
#     labs(x = "row", y = "col")
# ggsave("./images/psiMat.png", gPsi, "png",
#        width = 6, height = 6, units = "in", dpi = 150)
# 