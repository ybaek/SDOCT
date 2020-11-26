# Sketch of DP inference components
dmvnorm <- mvtnorm::dmvnorm
dmvt <- mvtnorm::dmvt
forbacksolve <- function(r, x) backsolve(r, forwardsolve(t(r), x))

# TODO: DP label clustering moves (quite hard one to crack!)

impute_car <- function(TARGET, means, var, rho, mis_inds, adj_l) {
  # Missing (at random) for conditionally autoregressive model
  # Assumes Leroux model -- precision mixture so that pos.def.
  for (j in 1:dim(mis_inds)[1]) {
    # iterates across columns of data
    mis_loc <- mis_inds[j,2]
    neighbors <- adj_l[[j]]
    nj <- length(neighbors)
    for (i in mis_inds[,1]) {
      # iterates across rows of data
      TARGET[i,mis_loc] <- (rho*sum(TARGET[i,neighbors])+(1-rho)*means[i,mis_loc]) / (nj*rho + 1-rho) + 
        (var / (nj*rho + 1-rho))^.5 * rnorm(1)
    }
  }
  return(TARGET)
}
