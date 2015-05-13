# March 24, 2015
# Read in SIPP .dta files (merged A, B, and D sets from CEPR)
# Select relevant variables
# Save files into ./Data directory as RData format
library(foreign)
library(dplyr)

# specify which variables to keep from CEPR sets A, B, and D
keepVars <- c("age", "educ", "female", "id", "race", "month", 
              "wpfinwgt", "year", "earnm", "occ", "job", "esr",
              "occ14", "ind23")

# input: full path to dta file with merged CEPR sets A, B, and D
# output: data frame with variables specified in keepVars and dates in R format
readCEPR <- function(filePath) {
        dta <- read.dta(filePath, convert.factors = FALSE)
        panel <- dta %>%
                select(one_of(keepVars)) %>%
                mutate(date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
                mutate(qtrdate = as.Date(paste(3*ceiling(month/3)-2, "/1/", year, sep =""), "%m/%d/%Y")) %>%
                select(-year, -month)
        return(panel)
}

# 1996
sipp96 <- readCEPR("G:/Research_Analyst/Eubanks/SIPP/sippsets96/sippsets96ABD.dta")
saveRDS(sipp96, file("./Data/sipp96.RData"))
rm(sipp96)

# 2001
sipp01 <- readCEPR("G:/Research_Analyst/Eubanks/SIPP/sippsets01/sippsets01ABD.dta")
saveRDS(sipp01, file("./Data/sipp01.RData"))
rm(sipp01)

# 2004
sipp04 <- readCEPR("G:/Research_Analyst/Eubanks/SIPP/sippsets04/sippsets04ABD.dta")
saveRDS(sipp04, file("./Data/sipp04.RData"))
rm(sipp04)

# 2008
sipp08 <- readCEPR("G:/Research_Analyst/Eubanks/SIPP/sippsets08/sippsets08ABD.dta")
saveRDS(sipp08, file("./Data/sipp08.RData"))
rm(sipp08)