#
# Main script for fitting the linear mixed model
# (Script for ARVO submission paper)
# Last Updated: Feb. 10th, 21
#
library(rstan)
m <- stan_model(file = 'constrained_nnet.stan')
