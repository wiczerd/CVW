# March 24, 2015
# Read in SIPP .dta files (merged A, B, and D sets from CEPR)
# Select relevant variables
# Save files into ./Data directory as RData format
# Clear memory, repeat
library(foreign)
library(dplyr)

keepVars <- c("age", "educ", "female", "id", "race", "month", 
              "wpfinwgt", "year", "earnm", "occ", "job", "esr")

# 1996
sipp96ABD <- read.dta("G:/Research_Analyst/Eubanks/SIPP/sippsets96/sippsets96ABD.dta", convert.factors = FALSE)
sipp96 <- sipp96ABD %>%
        select(one_of(keepVars)) %>%
        mutate(date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        mutate(qtrdate = as.Date(paste(3*ceiling(month/3)-2, "/1/", year, sep =""), "%m/%d/%Y")) %>%
        select(-year, -month)
saveRDS(sipp96, file("./Data/sipp96.RData"))
rm(list = setdiff(ls(), "keepVars"))

# 2001
sipp01ABD <- read.dta("G:/Research_Analyst/Eubanks/SIPP/sippsets01/sippsets01ABD.dta", convert.factors = FALSE)
sipp01 <- sipp01ABD %>%
        select(one_of(keepVars)) %>%
        mutate(date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        mutate(qtrdate = as.Date(paste(3*ceiling(month/3)-2, "/1/", year, sep =""), "%m/%d/%Y")) %>%
        select(-year, -month)
saveRDS(sipp01, file("./Data/sipp01.RData"))
rm(list = setdiff(ls(), "keepVars"))

# 2004
sipp04ABD <- read.dta("G:/Research_Analyst/Eubanks/SIPP/sippsets04/sippsets04ABD.dta", convert.factors = FALSE)
sipp04 <- sipp04ABD %>%
        select(one_of(keepVars)) %>%
        mutate(date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        mutate(qtrdate = as.Date(paste(3*ceiling(month/3)-2, "/1/", year, sep =""), "%m/%d/%Y")) %>%
        select(-year, -month)
saveRDS(sipp04, file("./Data/sipp04.RData"))
rm(list = setdiff(ls(), "keepVars"))

# 2008
sipp08ABD <- read.dta("G:/Research_Analyst/Eubanks/SIPP/sippsets08/sippsets08ABD.dta", convert.factors = FALSE)
sipp08 <- sipp08ABD %>%
        select(one_of(keepVars)) %>%
        mutate(date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        mutate(qtrdate = as.Date(paste(3*ceiling(month/3)-2, "/1/", year, sep =""), "%m/%d/%Y")) %>%
        select(-year, -month)
saveRDS(sipp08, file("./Data/sipp08.RData"))
rm(list = ls())