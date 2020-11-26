#
# A script collecting user functions for repeated usage
#

### Following are essentially for internal use only
### (no rigorous type checks are made, etc.)

logit <- function(p) return(log(p/(1-p)))

expit <- function(lo) return(1/(1+exp(-lo)))

logsumexp <- function(lp) {
  # Calculates log(sum(exp(lp))) in a stable way
  pv <- max(lp)
  return(pv + log(sum(exp(lp-pv))))
}

dlmvnorm <- function(x, mu, Sigma) {
  # density is on the log scale
  if (length(Sigma)==1) {
  # useful special case for linear regression
    return( -.5*(length(x)*(log(2*pi)+log(Sigma))+sum((x-mu)^2)/Sigma) )
  }
  Sigma <- as.matrix(Sigma)
  if (dim(Sigma)[1] != dim(Sigma)[2]) stop("Sigma must be a real number or a square matrix")
  L <- t(chol(Sigma))
  md <- sum(forwardsolve(L, x-mu)^2)
  return(-.5*(nrow(L)*log(2*pi)+md)-log(det(L)) )
}

forbacksolve <- function(r, x) {
  # Solves for y in equation (r^Tr)y = x
  # Assumes upper triangular r
  # For some reason transpose argument screws things up, so don't use that!
  if (r[2,1]) stop("r is not an upper triangular matrix")
  return(backsolve(r, forwardsolve(t(r), x)))
}

moranI <- function(x, P=8) {
  ## Calculates moran's I statistic on a P x P square grid
  x <- matrix(x, ncol = length(x))
  A <- adj_grid(P)
  return( 
    P^2/sum(A) * rowSums((x - mean(x)) %*% A * (x - mean(x))) / sum((x-mean(x))^2)
  )
}

wendland_c2 <- function(d, rho) {
  # Takes a distance matrix and returns a Wendland covariance function
  # matrix that is twice-differentiable
  # Function is compactly supported on the interval (0,rho)
  return( (1+4 * d/rho) * (1-d/rho)^4 * (d < rho) )
}
