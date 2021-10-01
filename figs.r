# Script for side effects from inference (a.k.a. figures)
source("infer.r")
#
png("images/roc_logistic.png")
plot(roc_model, col = 4)
plot(roc_raw, add = T, col = 2)
abline(v = .85, lty = 2)
legend(
    "bottomright",
    legend = c("Raw", "Model"),
    col    = c(2, 4),
    lty    = 1
)
dev.off()
#
raw_pcs <- prcomp(y_test2)$x # imputed copy for PCA
smoothed_pcs <- prcomp(f_test_full)$x # imputed copy for PCA
png("images/raw_pcs.png")
plot(raw_pcs, col = labels_test + 1, pch = 19)
text(x = raw_pcs[, 1], y = raw_pcs[, 2],
     labels = rownames(y_test), adj = c(1, 1), cex = .7,
     col = labels_test + 1)
dev.off()
png("images/smoothed_pcs.png")
plot(smoothed_pcs, col = labels_test + 1, pch = 19)
text(x = smoothed_pcs[, 1], y = smoothed_pcs[, 2],
     labels = rownames(y_test), adj = c(1, 1), cex = .7,
     col = labels_test + 1)
dev.off()