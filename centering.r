# A set of "centering" statistics are derived from
# both the samples and the held-out data
gcl <- read.csv("./data/gcl.csv") # includes ILM
ipl <- read.csv("./data/ipl.csv") # includes GCL
cpRnfl <- read.csv("./data/rnfl.csv")
healthy_mrns_c <- read.csv("./data/healthy_ids.csv", header = FALSE)$V1

for (i in 1:8) {
  gcl[gcl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)] <-
  rev(gcl[gcl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)])
  ipl[ipl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)] <-
  rev(ipl[ipl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)])
}

idMatches <- intersect(
  intersect(unique(ipl$maskedid), unique(gcl$maskedid)),
  unique(cpRnfl$PatientID)
) # 506 ID's
idMatches_ctr <- intersect(idMatches, healthy_mrns_c) # Healthy control (64)
incl_t <- paste("RNFLT", 1:768, sep = ".")
incl_col1 <- colnames(ipl)[grep("^mean", colnames(ipl))]
incl_col2 <- colnames(gcl)[grep("^mean", colnames(gcl))]
# Consider the sample average / quantiles
m_ctr_avg <- apply(ipl[ipl$maskedid %in% idMatches_ctr, incl_col1], 2, function(x) mean(as.numeric(x), na.rm = TRUE)) +
  apply(gcl[gcl$maskedid %in% idMatches_ctr, incl_col2], 2, function(x) mean(as.numeric(x), na.rm = TRUE))
cp_ctr_avg <- apply(cpRnfl[cpRnfl$PatientID %in% idMatches_ctr, incl_t], 2, function(x) mean(as.numeric(x), na.rm = TRUE))

m_ctr_q5 <- apply(ipl[ipl$maskedid %in% idMatches_ctr, incl_col1], 2, function(x) quantile(as.numeric(x), prob = .05, na.rm = TRUE)) +
  apply(gcl[gcl$maskedid %in% idMatches_ctr, incl_col2], 2, function(x) quantile(as.numeric(x), prob = .05, na.rm = TRUE))
cp_ctr_q5 <- apply(cpRnfl[cpRnfl$PatientID %in% idMatches_ctr, incl_t], 2, function(x) quantile(as.numeric(x), prob = .05, na.rm = TRUE))

m_stats <- rbind(m_ctr_avg, m_ctr_q5)
cp_stats <- rbind(cp_ctr_avg, cp_ctr_q5)
rownames(m_stats) <- c("mean", "q5")
colnames(m_stats) <- NULL
rownames(cp_stats) <- c("mean", "q5")
colnames(cp_stats) <- NULL
saveRDS(list(cp = cp_stats, macula = m_stats), file = "data/stats.rds")
