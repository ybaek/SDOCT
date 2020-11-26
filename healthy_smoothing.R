# Over-smoothing healthy patients' mean for modeling
# (Clinicians model deviations, which are obtained after subtracting "healthy mean"
# from raw thicknesses. In practice, this prior info. comes from proprietary database.
# We can approximate the procedure by replacing the sample means with some approximate values)

source("./adj_grid.R") # Code for conditionally autoregressive modeling
source("./utilities.R") # Code for a special GP kernel
source("./prior_distances.R") # Code computing the distance matrix
# L_RANGE <- 135/360 * 2*pi # Constant restricting the range of GP bandwidth
#####

# Rehash from data_processing.R
gcl <- read.csv("../data/gcl.csv")
ipl <- read.csv("../data/ipl.csv")
cpRnfl <- read.csv("../data/rnfl.csv")
healthy_mrns_c <- read.csv("../data/healthy_ids.csv", header = FALSE)$V1

idMatches <- intersect(
  intersect(unique(ipl$maskedid), unique(gcl$maskedid)),
  unique(cpRnfl$PatientID)
)
idMatches_ctr <- intersect(idMatches, healthy_mrns_c) 
idMatches_tr <- setdiff(idMatches, healthy_mrns_c) 
incl.t <- paste("RNFLT", c(1:128, 609:768), sep = ".")

#####

# Raw sample averages
healthy_i <- colMeans(  
  apply(ipl[ipl$maskedid %in% idMatches_ctr, grep("^mean", colnames(ipl))], 2, as.numeric ), 
  na.rm = TRUE)
healthy_g <- colMeans(  
  apply(gcl[gcl$maskedid %in% idMatches_ctr, grep("^mean", colnames(gcl))], 2, as.numeric), 
  na.rm = TRUE)
healthy_Y <- healthy_i + healthy_g
healthy_Z <- colMeans(  
  apply(cpRnfl[cpRnfl$PatientID %in% healthy_mrns_c, incl.t], 2, as.numeric), 
  na.rm = TRUE)

# Very simple idea: Why not just take +/- 3SD of pointwise means?
# (3 not 2 seems still okay because there seems to be a considerable 
# heterogeneity across these "databases")
Y_sd <- apply(ipl[ipl$maskedid %in% idMatches_ctr, grep("^mean", colnames(ipl),)], 
      2, function(x) sd(x, na.rm = TRUE)) + 
        apply(gcl[gcl$maskedid %in% idMatches_ctr, grep("^mean", colnames(ipl),)], 
              2, function(x) sd(x, na.rm = TRUE))
Z_sd <- apply(apply(cpRnfl[cpRnfl$PatientID %in% healthy_mrns_c, incl.t], 2, as.numeric), 
      2, function(x) sd(x, na.rm = TRUE))
set.seed(2020-11-21)
healthy_Z2 <- healthy_Z + 1.5 * rnorm(length(healthy_Z)) * Z_sd / sqrt(length(healthy_Z))
healthy_Y2 <- healthy_Y + 1.5 * rnorm(length(healthy_Y)) * Y_sd / sqrt(length(healthy_Y))

saveRDS(healthy_Y2, file = "../data/control_m.Rds")
saveRDS(healthy_Z2, file = "../data/control_cp.Rds")
