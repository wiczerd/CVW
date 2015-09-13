# April 27, 2015
# Modify calculateResidWages and inflation adjusts to create quarterly total earnings and use those
# to calculate wage changes.
# Save analytic data in ./Data/ directory.
# Precondition: processData.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)
library(xlsx)

setwd("~/workspace/CVW/R")
#setwd("G:/Research_Analyst/Eubanks/Occupation Switching Lab/")
source("./RevArt/functions.R")

useRegResid <- TRUE
# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- TRUE

# import PCE deflator
PCE <- read.csv("./Data/PCE.csv")
PCE$date <- as.Date(PCE$date)

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
				   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
	mutate(month = format(date, "%m"),
		   year = format(date, "%Y"),
		   date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
	select(-year, -month)

setwd("./Data")

detach(package:xlsx)

# Setup some arrays ---------------------------------------------------------------


regressors <- c("age", 
                "educ", 
                "female", 
                "race", 
                "yearsSchool",
                "experience",
                "black", 
                "hispanic", 
                "year", 
                "earnm", 
                "logEarnm")


toKeep <- c("id",
	    "wpfinwgt", 
	    "switchedOcc",
	    "soc2d", 
	    "occ",
	    "nextOcc",
	    "EE", 
	    "UE",
	    "EU",
		"Young",
		"HSCol",
	    "wageChange", 
	    "wageChange_EUE", 
	    "wageChange_all",
	    "lfStat", 
	    "date",
	    "occWage",
	    "occWageChange",
	    "useWage",
	    "nextWage",
	    "nextOccWage",
	    "recIndic",
	    "waveRec")

# Combine panels  --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed9608 <- readRDS("processed96.RData")

processed01 <- readRDS("processed01.RData")
processed9608 <- bind_rows(processed01, processed9608)
rm(processed01)

processed04  <- readRDS("processed04.RData")
processed9608 <- bind_rows(processed04, processed9608)
rm(processed04)

processed08  <- readRDS("processed08.RData")
processed9608 <- bind_rows(processed08, processed9608)
rm(processed08)

processed9608 <- left_join(processed9608, PCE, by = "date")
rm(PCE)

# Do the regressions -------------------------------------------

# Generate regressor variables
analytic9608 <- genRegressors(processed9608)
rm(processed9608)

# Find average log wage for 1996 panel
# avg19962008 <- weighted.mean(analytic9608$logEarnm, analytic9608$wpfinwgt, na.rm = TRUE)
logEarnmXweights <- analytic9608$logEarnm * analytic9608$wpfinwgt
weights <- analytic9608$wpfinwgt
avg19962008 <- sum(logEarnmXweights, na.rm = TRUE)/sum(weights, na.rm = TRUE)
rm(list = c("weights", "logEarnmXweights"))

# Run regression for whole sample, remove regressors
analytic9608 <- calculateUseWage(analytic9608, avg19962008)
analytic9608 <- calculateOccWage(analytic9608, 0.)
analytic9608$Young <- (analytic9608$age <30)
analytic9608$Young <- ifelse(is.na(analytic9608$age),NA,analytic9608$Young)
analytic9608$HSCol <- (analytic9608$educ>=4) + (analytic9608$educ>=2)
analytic9608$HSCol <- ifelse(is.na(analytic9608$educ),NA,analytic9608$HSCol)
analytic9608 <- select(analytic9608, -one_of(regressors))
analytic9608 <- analytic9608 %>%
	mutate(useWageLevel = 0.5*(exp(useWage) - exp(-useWage))) %>%
	mutate(useWageLevel = ifelse(lfStat > 1 , 0. , useWageLevel)) #wage if not working is 0 in levels
	# mutate(useWage  = ifelse(lfStat > 1 , 1. , useWage)) #wage if not working is inv hyperbolic sine(0)

# Create quarterly residual wage and quarterly lfStat
analytic9608 <- analytic9608 %>%
	filter(!is.na(lfStat)) %>%
	group_by(id, qtrdate) %>%
	mutate(useWageQtr = sum(useWageLevel)) %>%
	mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
	mutate(useWageQtr = log(useWageQtr + sqrt(useWageQtr^2 + 1)))

# Calculate last residual wage, fill up
analytic9608 <- fillUpWage(analytic9608)

#fix SwitchedOcc for UE
# Occupation is already filled up; these commands do not change switchedOcc
#analytic9608$switchedOcc <- ifelse(analytic9608$UE, (analytic9608$occ != analytic9608$nextOcc), analytic9608$switchedOcc)
#analytic9608$switchedOcc <- ifelse(analytic9608$EU, (analytic9608$occ != analytic9608$nextOcc), analytic9608$switchedOcc)

# Calculate residual wage change
analytic9608 <- calculateWageChange(analytic9608)

# Save data, remove from environment
if(useRegResid) {
	setwd("./RegResid")
}else{
	setwd("./Raw")
}

saveRDS(analytic9608, "analytic9608.RData")

# Set up wage changes data ----------------------------------------

wageChanges <- analytic9608 %>%
	select(one_of(toKeep))
rm(analytic9608)

#merge in unemployment
wageChanges <- left_join(wageChanges, haver, by = "date")

# throw out infinity and missing values
wageChanges <- wageChanges %>%
	# drop the ones with no wage change (i.e missing values)
	filter(!is.infinite(wageChange) & !is.na(wageChange))

# store full set
saveRDS(wageChanges, "wageChanges.RData")


