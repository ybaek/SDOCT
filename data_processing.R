#
# Preprocessing steps: LAST UPDATED Jan. 31, 2021
#
# Remember Alessandro's comment:
# "There is sth to be said about either GCL ONLY or GCC+
# (Spectralis has higher resolution but we go with previous conventions)"
#
# Sam's idea:
# "We want to keep this strictly cross-sectional--
# for now, just single image per patient, the most recent pair of macula + RNFL"
#
# Dr. Medeiros' point:
# "We want to keep both right and left eye. No reason NOT to!!!"
# (Prob. Sam mentioned this too at some point but it got missed on me)
#
gcl <- read.csv("./data/gcl.csv") # includes ILM
ipl <- read.csv("./data/ipl.csv") # includes GCL
cpRnfl <- read.csv("./data/rnfl.csv")
healthy_mrns_c <- read.csv("./data/healthy_ids.csv", header = FALSE)$V1
# Macula: need to flip the grid if the eye is left!
for (i in 1:8) {
  gcl[gcl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)] <-
  rev(gcl[gcl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)])
  ipl[ipl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)] <-
  rev(ipl[ipl$eye == "L", (12 + 8 * (i - 1) + 1):(12 + 8 * i)])
}

# Joining the tables by Subject ID
idMatches <- intersect(
  intersect(unique(ipl$maskedid), unique(gcl$maskedid)),
  unique(cpRnfl$PatientID)
) # 506 ID's
idMatches_ctr <- intersect(idMatches, healthy_mrns_c) # Healthy control (64)
idMatches_tr <- setdiff(idMatches, healthy_mrns_c) # Treatment group (442)

# Pre-processing cpRNFL and GCL
# (Indices included: cf. Fujino et al. 2018)
# UPDATE: We may include some points in the superior quadrant too,
# so including 288 points may not be the best practice!
incl_t <- paste("RNFLT", 1:768, sep = ".")
incl_col1 <- colnames(ipl)[grep("^mean", colnames(ipl))]
incl_col2 <- colnames(gcl)[grep("^mean", colnames(gcl))]

# Raw sample averages
# UPDATE: instead of averaging, take 95% quartiles?? (Dr. Medeiros' idea)
# Rationale: We want SIGNIFICANT DEVIATIONS
ipl_sample    <- subset(ipl, maskedid %in% idMatches)
gcl_sample    <- subset(gcl, maskedid %in% idMatches)
cpRnfl_sample <- subset(cpRnfl, PatientID %in% idMatches)
ipl_ctr    <- subset(ipl_sample, maskedid %in% idMatches)[, incl_col1]
gcl_ctr    <- subset(gcl_sample, maskedid %in% idMatches)[, incl_col2]
cpRnfl_ctr <- subset(cpRnfl_sample, PatientID %in% idMatches)[, incl_t]
rm(list = c("ipl", "gcl", "cpRnfl"))

# An array of statistics that may be considered as "centering"
# from thicknesses to deviations (arbitrary!)
# Will need to separate this code -- for now, leaving this as blank
cpRnfl_sample <- cpRnfl_sample[, c("PatientID", "ExamDate", "Eye", incl_t)]
cpRnfl_sample[, incl_t] <- apply(cpRnfl_sample[, incl_t], 2, as.numeric)
ipl_sample <- ipl_sample[, c("maskedid", "examdate", "eye", incl_col1)]
ipl_sample[, incl_col1] <- apply(ipl_sample[, incl_col1], 2, as.numeric)
gcl_sample <- gcl_sample[, c("maskedid", "examdate", "eye", incl_col2)]
gcl_sample[, incl_col2] <- apply(gcl_sample[, incl_col2], 2, as.numeric)

# Let's change the column names with no weird characters
# also match conventions between RNFL and GCL+IPL
colnames(cpRnfl_sample)[1:3] <- c("patientid", "examdate", "eye")
colnames(gcl_sample)[1] <- c("patientid")
colnames(ipl_sample)[1] <- c("patientid")
incl_t_mod <- unlist(lapply(strsplit(incl_t, "T\\."),
  function(str) paste(str, collapse = "")))
# Weird character UTF-8: raw hex c3 a2
incl_col1_mod <- unlist(lapply(strsplit(incl_col1,
  rawToChar(charToRaw("â"))), function(str) str[1]))
incl_col2_mod <- unlist(lapply(strsplit(incl_col2,
  rawToChar(charToRaw("â"))), function(str) str[1]))
colnames(cpRnfl_sample)[-c(1:3)] <- tolower(incl_t_mod)
colnames(gcl_sample)[-c(1:3)] <- incl_col1_mod
colnames(ipl_sample)[-c(1:3)] <- incl_col2_mod

# Change date/times to POSIX-compliant formats (EST)
# Get rid of dates older than the last RNFL measurement
# Since ipl and gcl are exported from the same file,
# they share all patient ID's and dates
cpRnfl_sample$examdate <- as.Date(cpRnfl_sample$examdate, try = "%m/%d/%Y", tz = "EST")
ipl_sample$examdate <- as.Date(ipl_sample$examdate, try = "%m/%d/%Y", tz = "EST")
gcl_sample$examdate <- as.Date(gcl_sample$examdate, try = "%m/%d/%Y", tz = "EST")

# New dataset: most recent pair from each patient (Left AND Right)
gcl_sample[, -c(1:3)] <- gcl_sample[, -c(1:3)] + ipl_sample[, -c(1:3)]
cpRnfl_dates <- aggregate(cpRnfl_sample$examdate,
  by = list(patientid = cpRnfl_sample$patientid), max)
gcl_dates <- aggregate(gcl_sample$examdate,
  by = list(patientid = gcl_sample$patientid), max)
# Take images from the most recent dates (can be multiple)
colnames(cpRnfl_dates)[2] <- "examdate"
colnames(gcl_dates)[2] <- "examdate"
# No pairing if more than 6 months apart
cpRnfl_recents <- merge(cpRnfl_dates, cpRnfl_sample)
gcl_recents <- merge(gcl_dates, gcl_sample)
my_mean <- function(x) {
  if(!length(na.omit(x))) return(NA)
  return(mean(x, na.rm = T))
}
cp_cs <- aggregate(cpRnfl_recents[, -c(1:3)],
  list(patientid = cpRnfl_recents$patientid,
       examdate  = cpRnfl_recents$examdate,
       eye       = cpRnfl_recents$eye),
  my_mean)
gcl_cs <- aggregate(gcl_recents[, -c(1:3)],
  list(patientid = gcl_recents$patientid,
       examdate  = gcl_recents$examdate,
       eye       = gcl_recents$eye),
  my_mean)
cs_data <- merge(cp_cs, gcl_cs, by = c("patientid", "eye"))
cs_data <- subset(cs_data, examdate.y - examdate.x <= 180)

## 731 eyes (=images by our design)
## from 393 patients
## 11 patient missing the left eye (subunit)
## 615 eyes glaucoma, 116 eyes healthy
## (333 patients glaucoma, 60 healthy)
cs_data$group <- 0L
cs_data$group[cs_data$patientid %in% idMatches_tr] <- 1L
colnames(cs_data)[3] <- "examdate"
cs_data <- cs_data[, -which(colnames(cs_data) == "examdate.y")]
# Save in form of several objects
dataset <- list(
  Z = as.matrix(cs_data[, grep("rnfl", colnames(cs_data))]),
  Y = as.matrix(cs_data[, grep("mean", colnames(cs_data))]),
  id = as.integer(as.factor(cs_data$patientid)),
  group = cs_data$group
)
saveRDS(dataset, file = "./data/dataset.Rds")
