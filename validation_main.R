library(ggplot2)
library(reshape2)
library(scales)
source("validation.R")
source("prior_distances.R")

# Downsampling prior mean
prior_mean <- prior_mean[seq(1, 288, 8), ]

# Plotting the coefficient
samples <- readRDS("./objects/samples_lmm_cpp.Rds")
pm <- matrix(rowMeans(samples[1:2304, ]), 36)
g <- ggplot() +
    geom_tile(data = reshape2::melt(pm),
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
preds_lss <- as.data.frame(data$Z_s %*% pm)
lss_rocs <- validate_roc(preds_lss, data$labels,
                         data$train, data$test, data$Z_s, data$Y_s)
lss_rocs$pauc # 39%
# Another decision rule
# No logistic regression, just make a decision based on
# negative deviations
# Raw vs. smoothed (w/ cpRNFL info.)
raw_count_roc <- pROC::roc(data$labels[data$test] ~ counts_raw)
raw_count_proc <- pROC::roc(data$labels[data$test] ~ counts_raw,
                            partial.auc = c(.85, 1),
                            partial.auc.correct = F,
                            partial.auc.focus = "sp")
pROC::auc(raw_count_proc,
          partial.auc = c(.85, 1),
          partial.auc.correct = F,
          partial.auc.focus = "sp") / .15 # 25%
lss_count_roc <- pROC::roc(data$labels[data$test] ~ counts_smoothed)
lss_count_proc <- pROC::roc(data$labels[data$test] ~ counts_smoothed,
                            partial.auc = c(.85, 1),
                            partial.auc.correct = F,
                            partial.auc.focus = "sp")
pROC::auc(lss_count_proc,
          partial.auc = c(.85, 1),
          partial.auc.correct = F,
          partial.auc.focus = "sp") / .15 # 23%
# Plot
png("./images/roc.png", width = 6, height = 6, units = "in", res = 120)
plot(lss_rocs$roc, main = "ROC curves", col = 2)
lines(raw_count_roc, col = 1)
lines(lss_count_roc, col = 4)
legend("bottomright",
       legend = c("SDM (71% / 25%)",
                  "Smoothed SDM (76% / 23%)",
                  "Smoothed + Logistic (67% / 39%)"),
       lty = 1, col = c(1, 4, 2))
dev.off()

# Inference: on other possible model parameters
# Why we can diagnose ill model fit
c_samples <- samples[(2304 + 1 + 1):(2304 + 1 + 1712), ]^.5
c_densities <- apply(c_samples, 1, function(x) density(x, from = 0))
plot_xmax <- max(sapply(c_densities, function(y) max(y$x)))
plot_ymax <- max(sapply(c_densities, function(x) max(x$y)))
png("./images/densities.png", width = 6, height = 6, units = "in", res = 120)
plot(c_densities[[1]],
     xlab = "c",
     ylab = "density",
     main = "",
     xlim = c(0, plot_xmax),
     ylim = c(0, plot_ymax),
     col = scales::alpha(1, .2))
for (i in 2:1712) {
    lines(c_densities[[i]], col = scales::alpha(1, .2))
}
dev.off()

psiMean <- matrix(0, 36, 36)
psiMean[lower.tri(psiMean, diag = T)] <- rowMeans(
    samples[(2304 + 1 + 1712 + 1):(2304 + 1 + 1712 + 666), ]
)
psiMean[upper.tri(psiMean)] <- t(psiMean)[upper.tri(psiMean)]

gPsi <- ggplot(reshape2::melt(solve(psiMean))) +
    geom_tile(aes(Var1, Var2, fill = value)) +
    labs(x = "row", y = "col")
ggsave("./images/psiMat.png", gPsi, "png",
       width = 6, height = 6, units = "in", dpi = 150)
