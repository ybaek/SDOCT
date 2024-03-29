### Following are essentially for internal use only
### (no rigorous type checks are made, etc.)

logit <- function(p) return(log(p/(1-p)))

expit <- function(lo) return(1/(1+exp(-lo)))

logsumexp <- function(lp) {
  # Calculates log(sum(exp(lp))) in a stable way
  pv <- max(lp)
  return(pv + log(sum(exp(lp-pv))))
}

form_kernel_matrix <- function(X, phi) {
  # Forms a kernel matrix based on predictors X for nonlinear regression
  # Tuning parameter: phi > 0
    exp(-as.matrix(dist(X))^2 * phi)
}

## Code to map an NxN regular rectangular grid to an adjacency matrix
adj_grid <- function(N, queen = FALSE, order = 1) {
  # Returns an N^2 x N^2 binary adjacency matrix
  # Nodes of a raster are numbered from left to right & scrolls after reaching N multiples
  # N: number of dimensions of a rectangular grid (positive integer)
  # queen: whether rectangles are adjacent when they just share a vertex (default is FALSE; "rook's case")
  # order: to what order are nodes considered adjacent (FOR FUTURE POSSIBLE USE)
  if (N!=as.integer(N) || N<1 ) stop("argument N is not a valid raster dimension")
  result <- matrix(0, nrow = N*N, ncol = N*N)
  for ( i in 1:(N*N) ) {
    i_north <- (i-N) > 0 # Does any node  exist to the north?
    i_south <- (i+N) <= N*N # to the south?
    i_west <- (i-1)%%N != 0 # to the west?
    i_east <- i%%N != 0 # to the east?
    if ( i_north ) result[i,i-N] <- 1 
    if ( i_south ) result[i,i+N] <- 1 
    if ( i_west ) result[i,i-1] <- 1 
    if ( i_east ) result[i,i+1] <- 1 
    if ( i_north & i_west ) result[i,i-(N+1)] <- 1*queen
    if ( i_north & i_east ) result[i,i-(N-1)] <- 1*queen
    if ( i_south & i_west ) result[i,i+(N-1)] <- 1*queen
    if ( i_south & i_east ) result[i,i+(N+1)] <- 1*queen
  }
  return(result)
}