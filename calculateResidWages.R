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

useRegResid <- T
# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T

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

# Functions ---------------------------------------------------------------
# Create function to calculate wage to use for analsysis
# add a constant (previously was constant from regression, now is 1996 avg)
# useWage is wage to use (previously always called resid, even if not a residual)
calculateUseWage <- function(df, const = 0) {
        if(useRegResid) {
                model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
                                    female + black + hispanic + factor(soc2d), data = df,
                            na.action = na.exclude, weights = wpfinwgt)
                useWage <- residuals(model) + const
                result <- data.frame(df, useWage)
        } else {
                result <- df %>%
                        mutate(useWage = logEarnm) 
        }
        return(result)
}

calculateOccWage <- function(df, const =0 ){
	if(useSoc2d){
		df <- group_by(df, occ)
	}else{
		df <- group_by(df, soc2d)
	}
	model <- lm( logEarnm ~ experience + I(experience^2) + factor(educ) 
				+ female + black + hispanic, 
				data= df, na.action = na.exclude, weights = wpfinwgt)
	occWage <- fitted(model) + const
	result <- data.frame(df,occWage)
	result <- result %>% 
		ungroup
	return(result)
}


# Create function to calculate last observed useWage
fillDownWage <- function(df) {
        result <- df %>%
                group_by(id) %>%
                arrange(id, date) %>%
                mutate(lastWage = as.numeric(ifelse(switchedJob & job != 0, lag(useWage), NA))) %>%
                mutate(lastWage = na.locf(lastWage, na.rm = FALSE)) %>%
                mutate(lastWageQtr = as.numeric(ifelse(switchedJob & job != 0, lag(useWageQtr, 3), NA))) %>%
                mutate(lastWageQtr = na.locf(lastWageQtr, na.rm = FALSE))
        result <- result %>%
                group_by(id) %>%
                arrange(id, date) %>%
                mutate(lastOcc = as.integer(ifelse(switchedJob & job != 0, lag(occ), NA))) %>%
                mutate(lastOcc = na.locf(lastOcc, na.rm = FALSE)) %>%
                mutate(lastOccWage = as.numeric(ifelse(switchedJob & job != 0, lag(occWage), NA))) %>%
                mutate(lastOccWage = na.locf(lastOccWage, na.rm = FALSE))  %>%
                ungroup
        return(result)
}

fillUpWage <- function(df) {
        result <- df %>%
                group_by(id) %>%
                arrange(id, date) %>%
                mutate(nextWage = as.numeric(ifelse(switchedJob & job != 0, lead(useWage), NA))) %>%
                mutate(nextWage = na.locf(nextWage, na.rm = FALSE, fromLast = TRUE)) %>%
                mutate(nextWageQtr = as.numeric(ifelse(switchedJob & job != 0, lead(useWageQtr, 3), NA))) %>%
                mutate(nextWageQtr = na.locf(nextWageQtr, na.rm = FALSE, fromLast = TRUE))
        result <- result %>%
                group_by(id) %>%
                arrange(date) %>%
                mutate(nextOccWage = as.numeric(ifelse(switchedJob & job != 0, lead(occWage), NA))) %>%
                mutate(nextOccWage = na.locf(nextOccWage, na.rm = FALSE, fromLast = TRUE)) %>%
                ungroup
        return(result)
}

# Create function to generate regressor variables and inflation adjusts
genRegressors <- function(df) {
        # import PCE data
        df <- left_join(df,PCE, by="date")
        result <- df %>%
                mutate(nomearnm = earnm) %>%
                mutate(earnm = earnm/PCEPI*100) %>%
                mutate(logEarnm = log(earnm)) %>%
                mutate(yearsSchool = as.integer(ifelse(educ == 1, 9, NA)),
                       yearsSchool = as.integer(ifelse(educ == 2, 12, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 3, 14, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 4, 16, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 5, 18, yearsSchool))) %>%
                mutate(experience = age - yearsSchool) %>%
                mutate(black = (race == 2)) %>%
                mutate(hispanic = (race == 3)) %>%
                mutate(year = as.numeric(format(date, "%Y")))
        return(result)
}

calculateWageChange <- function(df) {
        result <- df %>%
                mutate(wageChange = nextWage - lag(useWage)) %>%
                mutate(wageChangeQtr = nextWageQtr - lag(useWageQtr)) %>%
                mutate(occWageChange = nextOccWage - lag(occWage)) %>%
                mutate(wageChange_stayer = as.numeric(ifelse(!switchedJob, nextWage - lag(useWage), NA) )) %>%
                mutate(wageChange_all = as.numeric(ifelse(!switchedJob, wageChange_stayer, wageChange)))
#                 mutate(wageChange = lead(useWage) - lastWage) %>%
#                 mutate(wageChangeQtr = lead(useWageQtr) - lastWageQtr) %>%
#                 mutate(occWageChange = lead(occWage) - lastOccWage)) %>%
#                 mutate(wageChange_stayer = as.numeric(ifelse(!switchedJob, lead(useWage) - lastWage, NA) )) %>%
#                 mutate(wageChange_all = as.numeric(ifelse(!switchedJob, wageChange_stayer, wageChange)))
        return(result)
}

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
			"EE", 
			"UE",  
			"wageChange", 
			"wageChange_wU", 
			"wageChange_stayer",
			"lfStat", 
			"date",
			"wageChangeQtr", 
			"wageChangeQtr_wU",
			"occWage",
			"occWageChange",
			"useWage",
			"nextWage",
			"nextOccWage",
			#             "lastWage",
			#             "lastOccWage",
			"recIndic",
			"waveRec")

# Combine panels  --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed9608 <- readRDS("processed96.RData")
processed01  <- readRDS("processed01.RData")
processed9608<-bind_rows(processed01,processed9608)
rm(processed01)
processed04  <- readRDS("processed04.RData")
processed9608<-bind_rows(processed04,processed9608)
rm(processed04)
processed08  <- readRDS("processed08.RData")
processed9608<-bind_rows(processed08,processed9608)
rm(processed08)


# Do the regressions -------------------------------------------

# Generate regressor variables
analytic9608 <- genRegressors(processed9608)
rm(processed9608)
# Find average log wage for 1996 panel
avg19962008 <- weighted.mean(analytic9608$logEarnm, analytic9608$wpfinwgt, na.rm = TRUE)
# Run regression within each year, remove regressors
analytic9608 <- calculateUseWage(analytic9608, avg19962008)
analytic9608 <- calculateOccWage(analytic9608, 0.)
analytic9608 <- select(analytic9608, -one_of(regressors))
analytic9608 <- analytic9608 %>%
        mutate(useWageLevel = exp(useWage)) %>%
        mutate(useWageLevel = as.numeric(ifelse(lfStat > 1 , 0. , useWageLevel)))

# Create quarterly residual wage and quarterly lfStat
analytic9608 <- analytic9608 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(useWageQtr = sum(useWageLevel)) %>%
        mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
        mutate(useWageQtr = log(useWageQtr))

# Calculate last residual wage, fill down
#analytic96 <- fillDownWage(analytic96)
analytic9608 <- fillUpWage(analytic9608)

# Calculate residual wage change
analytic9608 <- calculateWageChange(analytic9608)

# job-job changes, including the separations. If separated,  pct change is -1
analytic9608 <- analytic9608 %>%
        mutate(wageChange_wU = as.numeric(ifelse(lfStat == 1 , wageChange , -1.))) %>%
        mutate(wageChangeQtr_wU = as.numeric(ifelse(lfStatQtr == 1 , wageChangeQtr , -1.)))

# Save data, remove from environment
if(useRegResid) {
        setwd("./RegResid")
} else {
        setwd("./Raw")
}

saveRDS(analytic9608, "analytic9608.RData")
rm(analytic9608)

setwd("../../")

# Set up wage changes data ----------------------------------------

setwd("./Data")

if(useSoc2d) {
	setwd("./soc2d")
} else {
	setwd("./occ")
}

if(useRegResid) {
	setwd("./RegResid")
} else {
	setwd("./Raw")
}

wageChanges <- analytic9608 %>%
	select(one_of(toKeep))
rm(analytic9608)

#merge in unemployment
wageChanges <- left_join(wageChanges, haver, by = "date")

# throw out infinity and missing values
wageChanges <- wageChanges %>%
	# drop the ones with no wage change (i.e missing values)
	filter(!is.infinite(wageChange) & !is.na(wageChange_wU) & !is.nan(wageChange_wU) )

# store full set
saveRDS(wageChanges, "wageChanges.RData")


