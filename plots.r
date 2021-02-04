# General plotting functions

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  ## Code for plotting multiple ggplot object print-outs
  ## (Taken from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
  ## (R Graphics CookBook, Winston Chang, CC0 licensed)
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) print(plots[[1]]) 
  else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))  
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# 1. Code for plotting covariance kernel density around the temporal
# for Gaussian process square-exp. kernel on the cpRNFL surface
# library(ggplot2)
# angles <- seq(-pi/2, pi/2, by = pi/64) # angular distances
# k <- c(.05, .5, 1, 2) # Bandwidth parameter
# df <- expand.grid(angles = angles, k = k)
# df$covariance <- exp( -(1 - cos(df$angles)) / (2 * df$k^2) )
# g <- ggplot(df) + 
#         geom_line(aes(angles, covariance)) + 
#         geom_vline(xintercept = 0, col = "red") + 
#         facet_wrap( ~ k, labeller = label_both) + 
#         labs(x = "angle (radian)", y = "correlation", title = "Metric: Cosine distance")    
# g + coord_flip() + scale_x_continuous(position = "top") + scale_y_reverse()

post_coef_plot <- function(coef, samples) {
  # Function for visualizing the coefficient posterior
  # Takes length P true coefficient (vector) and
  # an S x P matrix of coefficient samples (posterior draws)
  # Most of the idea comes from 
  # https://github.com/betanalpha/knitr_case_studies/blob/master/bayes_sparse_regression/plot_utility.R
  coef_ids <- rep(1:length(coef), each = 2)
  # Add/subtract .5 from each ID for histogram-like plotting
  plot_ids <- sapply(1:length(coef_ids), 
                     function(j) if (!(j%%2)) coef_ids[j] + .5 else coef_ids[j] - .5)  
  quants <- sapply(1:length(coef), 
                   function(j) quantile(samples[,j], probs = 1:9 * .1))
  # "Pad" the quantiles for polygon plotting
  pad_truths <- do.call(cbind, lapply(coef_ids, function(x) coef[x]))
  pad_quants <- do.call(cbind, lapply(coef_ids, function(x) quants[,x]))
  
  plot(1, type = "n", main = NULL,
       xlab = "Coefficient Index",
       ylab = "Posterior Draw",
       xlim = c(.5, length(coef)+.5),
       ylim = c(min(c(coef, quants[1,])),
                max(c(coef, quants[9,])))
  )
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[1,], rev(pad_quants[9,])),
          col = "#DCBCBC", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[2,], rev(pad_quants[8,])),
          col = "#C79999", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[3,], rev(pad_quants[7,])),
          col = "#B97C7C", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[4,], rev(pad_quants[6,])),
          col = "#A25050", border = NA)
  lines(plot_ids, pad_quants[5,], col = "#8F2727", lwd = 2)
  points(pch = 4, colMeans(samples), cex = .1) # posterior mean
  abline(h = 0)
  
  lines(plot_ids, pad_truths, col = "white", lwd = 1.5)
  lines(plot_ids, pad_truths, col = "black", lwd = 1.25)
}

post_resid_plot <- function(coef, samples) {
  # Function for visualizing the residual posterior
  # Takes length P true coefficient (vector) and
  # an S x P matrix of coefficient samples (posterior draws)
  # Most of the idea comes from 
  # https://github.com/betanalpha/knitr_case_studies/blob/master/bayes_sparse_regression/plot_utility.R
  samples <- sweep(samples, 2, coef)
  coef_ids <- rep(1:length(coef), each = 2)
  # Add/subtract .5 from each ID for histogram-like plotting
  plot_ids <- sapply(1:length(coef_ids), 
       function(j) if (!(j%%2)) coef_ids[j] + .5 else coef_ids[j] - .5)  
  quants <- sapply(1:length(coef), 
       function(j) quantile(samples[,j], probs = 1:9 * .1))
  # "Pad" the quantiles for polygon plotting
  pad_quants <- do.call(cbind, lapply(coef_ids, function(x) quants[,x]))
  
  plot(1, type = "n", main = NULL,
       xlab = "Coefficient Index",
       ylab = "Posterior Draw Minus Truth",
       xlim = c(.5, length(coef)+.5),
       ylim = c(min(c(coef, quants[1,])),
                max(c(coef, quants[9,])))
       )
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[1,], rev(pad_quants[9,])),
          col = "#DCBCBC", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[2,], rev(pad_quants[8,])),
          col = "#C79999", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[3,], rev(pad_quants[7,])),
          col = "#B97C7C", border = NA)
  polygon(c(plot_ids, rev(plot_ids)), 
          c(pad_quants[4,], rev(pad_quants[6,])),
          col = "#A25050", border = NA)
  lines(plot_ids, pad_quants[5,], col = "#8F2727", lwd = 2)
  points(pch = 4, colMeans(samples), cex = .1) # posterior mean
  abline(h = 0, col = "white", lwd = 1.5)
  abline(h = 0, col = "black", lwd = 1.25)
}
