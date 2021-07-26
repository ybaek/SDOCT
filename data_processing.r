#
# Preprocessing steps: LAST UPDATED Jun. 28th, 2021
#
gcl <- read.csv("./data/gcl.csv") # includes ILM
ipl <- read.csv("./data/ipl.csv") # includes GCL
cpRnfl <- read.csv("./data/rnfl.csv")
demo <- read.csv("./data/20170322_Demographics.csv") # demographics
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
incl_t <- paste("RNFLT", 1:768, sep = ".")
incl_col1 <- colnames(ipl)[grep("^mean", colnames(ipl))]
incl_col2 <- colnames(gcl)[grep("^mean", colnames(gcl))]
# Raw sample 5% quartiles
ipl_sample    <- subset(ipl, maskedid %in% idMatches)
gcl_sample    <- subset(gcl, maskedid %in% idMatches)
cpRnfl_sample <- subset(cpRnfl, PatientID %in% idMatches)
ipl_ctr    <- subset(ipl_sample, maskedid %in% idMatches)[, incl_col1]
gcl_ctr    <- subset(gcl_sample, maskedid %in% idMatches)[, incl_col2]
cpRnfl_ctr <- subset(cpRnfl_sample, PatientID %in% idMatches)[, incl_t]
rm(list = c("ipl", "gcl", "cpRnfl"))

cpRnfl_sample <- cpRnfl_sample[, c("PatientID", "ExamDate", "Eye", incl_t)]
cpRnfl_sample[, incl_t] <- apply(cpRnfl_sample[, incl_t], 2, as.numeric)
ipl_sample <- ipl_sample[, c("maskedid", "examdate", "eye", incl_col1)]
ipl_sample[, incl_col1] <- apply(ipl_sample[, incl_col1], 2, as.numeric)
gcl_sample <- gcl_sample[, c("maskedid", "examdate", "eye", incl_col2)]
gcl_sample[, incl_col2] <- apply(gcl_sample[, incl_col2], 2, as.numeric)
# Let's change the column names with no weird characters
# and match conventions between RNFL and GCL+IPL
colnames(cpRnfl_sample)[1:3] <- c("maskedid", "examdate", "eye")
colnames(gcl_sample)[1] <- c("maskedid")
colnames(ipl_sample)[1] <- c("maskedid")
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
cpRnfl_sample$examdate <- as.Date(cpRnfl_sample$examdate, try = "%m/%d/%Y", tz = "EST")
ipl_sample$examdate <- as.Date(ipl_sample$examdate, try = "%m/%d/%Y", tz = "EST")
gcl_sample$examdate <- as.Date(gcl_sample$examdate, try = "%m/%d/%Y", tz = "EST")
# New dataset: most recent pair from each patient
gcl_sample[, -c(1:3)] <- gcl_sample[, -c(1:3)] + ipl_sample[, -c(1:3)]
cpRnfl_dates <- aggregate(cpRnfl_sample$examdate,
  by = list(maskedid = cpRnfl_sample$maskedid), max)
gcl_dates <- aggregate(gcl_sample$examdate,
  by = list(maskedid = gcl_sample$maskedid), max)
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
  list(maskedid = cpRnfl_recents$maskedid,
       examdate  = cpRnfl_recents$examdate,
       eye       = cpRnfl_recents$eye),
  my_mean)
gcl_cs <- aggregate(gcl_recents[, -c(1:3)],
  list(maskedid = gcl_recents$maskedid,
       examdate  = gcl_recents$examdate,
       eye       = gcl_recents$eye),
  my_mean)
cs_data <- merge(cp_cs, gcl_cs, by = c("maskedid", "eye"))
cs_data <- subset(cs_data, examdate.y - examdate.x <= 180)
# Group variable assignment: healthy vs. glaucomatose
cs_data$group <- 0L
cs_data$group[cs_data$maskedid %in% idMatches_tr] <- 1L
colnames(cs_data)[3] <- "examdate"
cs_data <- cs_data[, -which(colnames(cs_data) == "examdate.y")]
# Subsetting demographic variables to unique patient ID's
demo <- demo[demo$maskedid %in% unique(cs_data$maskedid), ]
demo <- demo[, c("maskedid", "dob", "gender", "race_nih", "race_primary")]
cs_data_demo <- merge(demo, cs_data, by = "maskedid")
cs_data_demo$age <- as.integer(format(cs_data_demo$examdate, "%Y")) -
  as.integer(
    format(as.Date(cs_data_demo$dob, try = "%m/%d/%Y", tz = "EST"), "%Y")
  )
# Save in form of several objects
dataset <- list(
  Z = as.matrix(cs_data[, grep("rnfl", colnames(cs_data))]),
  Y = as.matrix(cs_data[, grep("mean", colnames(cs_data))]),
  id = as.integer(as.factor(cs_data$maskedid)),
  group = cs_data$group
)
demo_df <- cs_data_demo[, c("age", "gender", "race_nih", "race_primary")]
saveRDS(dataset, file = "./data/dataset.rds")
saveRDS(demo_df, file = "./data/demo.rds")