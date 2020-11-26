#
# Preprocessing steps: LAST UPDATED Nov. 21, 2020
#
# Remember Alessandro's comment:
# 1. Is it sensible to average over the same year cpRNFL scans? (we have so few m scans that they aren't affected)
# 2. There is sth to be said about either GCL ONLY or GCC+ (Spectralis has higher resolution but we go with previous conventions)
#
library(lubridate)
gcl <- read.csv("../data/gcl.csv") # includes ILM
ipl <- read.csv("../data/ipl.csv") # includes GCL
cpRnfl <- read.csv("../data/rnfl.csv")
healthy_mrns_c <- read.csv("../data/healthy_ids.csv", header = FALSE)$V1
# Macula: need to flip the grid if the eye is left!
for (i in 1:8) {
  gcl[gcl$eye=="L", (12 + 8*(i-1) + 1):(12 + 8*i)] <-
    rev(gcl[gcl$eye=="L", (12 + 8*(i-1) + 1):(12 + 8*i)])
  ipl[ipl$eye=="L", (12 + 8*(i-1) + 1):(12 + 8*i)] <-
    rev(ipl[ipl$eye=="L", (12 + 8*(i-1) + 1):(12 + 8*i)])
}

# Joining the tables by Subject ID
idMatches <- intersect(
  intersect(unique(ipl$maskedid), unique(gcl$maskedid)),
  unique(cpRnfl$PatientID)
) # 506 ID's
idMatches_ctr <- intersect(idMatches, healthy_mrns_c) # Healthy control (64)
idMatches_tr <- setdiff(idMatches, healthy_mrns_c) # Treatment group (442)

# Pre-processing cpRNFL
# (Indices included: cf. Fujino et al. 2018)
incl.t <- paste("RNFLT", c(1:128, 609:768), sep = ".")
Z_sample <- subset(cpRnfl, cpRnfl$PatientID %in% idMatches)
Z_sample <- Z_sample[,c("PatientID", "ExamDate", incl.t)]
#nas <- sum(apply(Z_sample[,incl.t], 2, function(x) sum(x=="n/a")))
Z_sample[,incl.t] <- apply(Z_sample[,incl.t], 2, as.numeric) # n/a coerced into NA
#nas == sum(is.na(Z_sample))

# Pre-processing macula GCC
incl.col1 <- colnames(ipl)[grep("^mean", colnames(ipl))]
incl.col2 <- colnames(gcl)[grep("^mean", colnames(gcl))]
#
ipl_matched <- subset(ipl, ipl$maskedid %in% idMatches)
ipl_matched <- ipl_matched[, c("maskedid", "examdate", incl.col1)]
nas <- sum(apply(ipl_matched[,incl.col1], 2, function(x) sum(x=="n/a")))
ipl_matched[,incl.col1] <- apply(ipl_matched[,incl.col1], 2, as.numeric) # n/a coerced into NA
#nas == sum(is.na(ipl_matched))
#
gcl_matched <- subset(gcl, gcl$maskedid %in% idMatches)
gcl_matched <- gcl_matched[, c("maskedid", "examdate", incl.col2)]
nas <- sum(apply(gcl_matched[,incl.col2], 2, function(x) sum(x=="n/a")))
gcl_matched[,incl.col2] <- apply(gcl_matched[,incl.col2], 2, as.numeric) # n/a coerced into NA
#nas == sum(is.na(gcl_matched))

# Let's change the column names with no weird characters
# also match conventions between RNFL and GCL+IPL
colnames(Z_sample)[1:2] <- c("patientid", "examdate")
colnames(gcl_matched)[1] <- c("patientid")
colnames(ipl_matched)[1] <- c("patientid")
incl.t.mod <- unlist(lapply(strsplit(incl.t, "T\\."), function(str) paste(str, collapse = "")))
# Weird character UTF-8: hex C383 	
incl.col1.mod <- unlist(lapply(strsplit(incl.col1, "Ã"), function(str) str[1]))
incl.col2.mod <- unlist(lapply(strsplit(incl.col2, "Ã"), function(str) str[1]))
colnames(Z_sample)[-c(1:2)] <- tolower(incl.t.mod)
colnames(gcl_matched)[-c(1:2)] <- incl.col1.mod
colnames(ipl_matched)[-c(1:2)] <- incl.col2.mod

# Change date/times to POSIX-compliant formats (EST)
# Get rid of dates older than the last RNFL measurement
# Since ipl and gcl are exported from the same file, they share all patient ID's and dates
Z_sample$examdate <- as.POSIXct(Z_sample$examdate, try = "%m/%d/%Y", tz = "EST")
ipl_matched$examdate <- as.POSIXct(ipl_matched$examdate, try = "%m/%d/%Y", tz = "EST")
gcl_matched$examdate <- as.POSIXct(gcl_matched$examdate, try = "%m/%d/%Y", tz = "EST")
Z_sample <- subset(Z_sample, Z_sample$examdate >= min(gcl_matched$examdate))
ipl_matched <- subset(ipl_matched, ipl_matched$examdate <= max(Z_sample$examdate))
gcl_matched <- subset(gcl_matched, gcl_matched$examdate <= max(Z_sample$examdate))

####
# Averaging by CY: We have so few macula measurements that it doesn't matter for them
# Lots of lowering resolution, OTOH, for cpRNFL measurements

my_mean <- function(x) { if (!length(na.omit(x))) return(NA); mean(x, na.rm = TRUE)}
Z_byDate <- aggregate(Z_sample[,-c(1,2)], 
                      by = list(patientid = Z_sample$patientid, examyear = lubridate::year(Z_sample$examdate)), 
                      FUN = my_mean)
i_byDate <- aggregate(ipl_matched[,-c(1,2)],
                      by = list(patientid = ipl_matched$patientid, examyear = lubridate::year(ipl_matched$examdate)),
                      FUN = my_mean)
g_byDate <- aggregate(gcl_matched[,-c(1,2)],
                      by = list(patientid = gcl_matched$patientid, examyear = lubridate::year(gcl_matched$examdate)),
                      FUN = my_mean)
# Omit all NA's (No observation for that ID and CY)
i_byDate <- i_byDate[!apply(i_byDate, 1, function(x) sum(is.na(x))== 64),]
g_byDate <- g_byDate[!apply(g_byDate, 1, function(x) sum(is.na(x))== 64),]
# all(i_byDate$examyear == g_byDate$examyear & i_byDate$patientid == g_byDate$patientid)
# Matching the datasets with a natural (inner) join
# For IPL and GCL, all you do is just summing up the figures
m_byDate <- i_byDate
m_byDate[,-c(1,2)] <- i_byDate[,-c(1,2)] + g_byDate[,-c(1,2)] # There is 1 for which there is NA only for IPL
agg_data <- merge(Z_byDate, m_byDate, by = c("patientid", "examyear"))

# Now we code the groups: control vs. treatment
# Coding: 0 == control, 1 == treatment
agg_data$group <- 0
agg_data$group[agg_data$patientid %in% idMatches_tr] <- 1
## In total 625 eyes for unhealthy, 93 eyes for healthy
saveRDS(agg_data, file = "../data/macula_cross14-17.Rds")

#######
