### Classes for macula GCC and cpRNFL objects
source("adj_grid.R")
macula <- function(x) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  Q <- dim(x)[2]; Qsqrt <- as.integer(Q^.5)
  stopifnot(Q^.5 == Qsqrt)
  class(x) <- c("macula", class(x))
  attr(x, "Qsqrt") <- Qsqrt
  return(x)
}

cp <- function(x, d) {
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  stopifnot(dim(d)[1] == dim(d)[2])
  stopifnot(dim(d)[1] == dim(x)[2])
  class(x) <- c("cp", class(x))
  attr(x, "distances") <- d
  return(x)
}

### Class-specific methods
aggregate.cp <- function(CP, W, offset) {
  # reshapes the cpRNFL object by aggregating data in a window of W
  # distances are re-calculated with +offset from each starting point
  P <- dim(CP)[2]
  stopifnot(offset < W && offset >= 0)
  stopifnot(!(P %% W))
  inds <- seq(1+offset, P, by = W)
  D <- attr(CP, "distances"); CP <- as.matrix(CP)
  CP <- matrix(rowMeans(matrix( c(t(CP)), ncol = W, byrow = TRUE )), 
              ncol = P %/% W, byrow = TRUE)
  return( cp(CP, D[inds,inds]) )
}

plot.macula <- function(m, ind=1, print = FALSE) {
  # Plotting method for macula object
  library("ggplot2", quietly = TRUE)
  N <- attr(m, "Qsqrt")
  coords_m <- data.frame( h = rep(1:N, N), v = rep(1:N, each = N), 
                          y = as.numeric(m[ind,]) )
  g <- ggplot(coords_m) + geom_tile(aes(h, v, fill = y), col = "grey50") + 
    labs(title = "Macula GCC", fill = expression(paste(mu, "m")) ) + 
    theme( axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           axis.title.y = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank(),
           legend.position = "bottom",
           panel.background = element_rect(fill = "white", colour = NA) )
  if (print) print(g); invisible(g)
  return(g)
}

plot.cp <- function(cp, ind=1, print = FALSE) {
  # Plotting method for cpRNFL object
  # constants are determined by the finest scale cpRNFL scan avavilable
  library("ggplot2", quietly = TRUE)
  
  TOTAL_ANGLES_CONST <- 768; REL_ANGLES_CONST <- 288
  FIRST_THRES_IND_CONST <- 128; SECOND_THRES_IND_CONST <- 608
  
  ua <- 360 / TOTAL_ANGLES_CONST # in degrees
  W <- REL_ANGLES_CONST / dim(cp)[2]
  pinds <- c(seq(1, FIRST_THRES_IND_CONST), seq(SECOND_THRES_IND_CONST + 1, 768))
  coords_c <- data.frame( h = pinds*ua, 
                          v = 1, # placeholder for "wrapping" the square grid
                          z = rep( as.numeric(cp[ind,]), each = W) )
  temp_off <- 270/360 * 2*pi # offset for centering the temporal (radians)
  g <- ggplot(coords_c) + geom_tile(aes(h, v, fill = z)) + 
    coord_polar(theta = "x", start = temp_off) + 
    labs(title = "cpRNFL", fill = expression(paste(mu, "m")) ) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(), 
           axis.text.y = element_blank(), 
           axis.ticks.y = element_blank(),
           legend.position = "bottom",
           panel.background = element_rect(fill = "white", colour = NA) )
  if (print) print(g); invisible(g)
  return(g)
}
