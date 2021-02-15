#
# Main script for fitting the linear mixed model
# (Script for ARVO submission paper)
# Last Updated: Feb. 10th, 21
#
.utils <- new.env()
source("./utilities.r", local = .utils)
ds <- readRDS("./data/dataset.rds")
stats <- readRDS("./data/stats.rds")
