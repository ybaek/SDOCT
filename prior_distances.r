# A script that forms a binary connectivity matrix
# That can be used for centering the prior for 
# regression coefficient matrix

# Forming a dictionary of square-sector
d <- rep(
  list(
    c(12, 13), c(12, 13), c(12, 13), c(13, 14),
    c(13, 14), c(13, 14), c(1, 14, 15), c(1, 13, 14),
    4, c(2, 3, 4), c(1, 2), c(1, 2),
    4, c(3, 4), c(3, 4), c(3, 4)
  ), each = 2)
d_map <- rep(d[1:8], 2)
for (i in 2:4) d_map <- append( d_map, rep(d[(8*(i-1)+1):(8*i)], 2) )
names(d_map) <- 1:64

# Forming a named array of angles-sectors
ua <- 360 / 768 # unit angle
a <- c(seq(285/ua+1, 360/ua), seq(1, 60/ua))
names(a) <- c(rep(12, 32), rep(13, 32), rep(14, 32), rep(15, 32), 
              rep(1, 64), rep(2, 32), rep(3, 32), rep(4, 32))

prior_mean <- matrix(0, 64, 288)
rownames(prior_mean) <- as.character(1:64)
colnames(prior_mean) <- as.character(unname(a))

# Forming a connection matrix
for (i in c(1, 2, 3, 4, 12, 13, 14, 15)) {
  cind <- as.character(a[names(a) == i])
  rind <- names(d_map)[unlist(lapply(d_map, function(x) i %in% x))]
  prior_mean[rind, cind] <- 1
}
# Reordering the indices and final transposition
prior_mean <- t(prior_mean[,order(as.integer(colnames(prior_mean)))])

#########

# Forming an angular distance matrix
angles <- as.integer(rownames(prior_mean)) * ua * (2*pi/360) # convert to radians
distMat <- matrix(0, 288, 288) 
for (i in 1:(nrow(distMat)-1)) {
  distMat[i,i] <- 0
  for (j in (i+1):nrow(distMat)) {
    diff <- abs(angles[i] - angles[j])
    distMat[i,j] <- distMat[j,i] <- (diff <= pi) * diff + 
      (diff > pi) * (2*pi-diff)
  }
  distMat[j,j] <- 0
}

rm(list = c("d", "d_map", "a", "angles", "cind", "diff", "i", "j",
            "rind", "ua"))
gc()