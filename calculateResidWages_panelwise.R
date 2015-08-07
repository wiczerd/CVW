# April 27, 2015
# Modify calculateResidWages and inflation adjusts to create quarterly total earnings and use those
# to calculate wage changes.
# Save analytic data in ./Data/ directory.
# Precondition: processData.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)

setwd("~/workspace/CVW/R")

useRegResid <- T
# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T

# import PCE deflator
PCE <- read.csv("./Data/PCE.csv")
PCE$date <- as.Date(PCE$date)

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

# 1996 Panel --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed96 <- readRDS("processed96.RData")

# Generate regressor variables
analytic96 <- genRegressors(processed96)

# Find average log wage for 1996 panel
avg1996 <- weighted.mean(analytic96$logEarnm, analytic96$wpfinwgt, na.rm = TRUE)
# Run regression within each year, remove regressors
analytic96 <- calculateUseWage(analytic96, avg1996)
analytic96 <- calculateOccWage(analytic96, 0.)
analytic96 <- select(analytic96, -one_of(regressors))
analytic96 <- analytic96 %>%
        mutate(useWageLevel = exp(useWage)) %>%
        mutate(useWageLevel = as.numeric(ifelse(lfStat > 1 , 0. , useWageLevel)))

# Create quarterly residual wage and quarterly lfStat
analytic96 <- analytic96 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(useWageQtr = sum(useWageLevel)) %>%
        mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
        mutate(useWageQtr = log(useWageQtr))

# Calculate last residual wage, fill down
#analytic96 <- fillDownWage(analytic96)
analytic96 <- fillUpWage(analytic96)

# Calculate residual wage change
analytic96 <- calculateWageChange(analytic96)

# job-job changes, including the separations. If separated,  pct change is -1
analytic96 <- analytic96 %>%
        mutate(wageChange_wU = as.numeric(ifelse(lfStat == 1 , wageChange , -1.))) %>%
        mutate(wageChangeQtr_wU = as.numeric(ifelse(lfStatQtr == 1 , wageChangeQtr , -1.)))

# Save data, remove from environment
if(useRegResid) {
        setwd("./RegResid")
} else {
        setwd("./Raw")
}

saveRDS(analytic96, "analytic96.RData")
rm(list = c("processed96", "analytic96"))

setwd("../../")

# 2001 Panel --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed01 <- readRDS("processed01.RData")

# Generate regressor variables
analytic01 <- genRegressors(processed01)

# Run regression within each year, remove regressors
analytic01 <- calculateUseWage(analytic01, avg1996)
analytic01 <- calculateOccWage(analytic01,0.)
analytic01 <- select(analytic01, -one_of(regressors))

analytic01 <- analytic01 %>%
        mutate(useWageLevel = exp(useWage)) %>%
        mutate(useWageLevel = as.numeric(ifelse(lfStat > 1 , 0. , useWageLevel)))

# Create quarterly residual wage and quarterly lfStat
analytic01 <- analytic01 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(useWageQtr = sum(useWageLevel)) %>%
        mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
        mutate(useWageQtr = log(useWageQtr))

# Calculate last residual wage, fill down
# analytic01 <- fillDownWage(analytic01)
analytic01 <- fillUpWage(analytic01)

# Calculate residual wage change
analytic01 <- calculateWageChange(analytic01)

# job-job changes, including the separations. If separated,  pct change is -1
analytic01 <- analytic01 %>%
        mutate(wageChange_wU = as.numeric(ifelse(lfStat == 1 ,wageChange , -1.))) %>%
        mutate(wageChangeQtr_wU = as.numeric(ifelse(lfStatQtr == 1 ,wageChangeQtr , -1.)))

# Save data, remove from environment
if(useRegResid) {
        setwd("./RegResid")
} else {
        setwd("./Raw")
}

saveRDS(analytic01, "analytic01.RData")
rm(list = c("processed01", "analytic01"))

setwd("../../")

# 2004 Panel --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed04 <- readRDS("processed04.RData")

# Generate regressor variables
analytic04 <- genRegressors(processed04)

# Run regression within each year, remove regressors
analytic04 <- calculateUseWage(analytic04, avg1996)
analytic04 <- calculateOccWage(analytic04,0.)
analytic04 <- select(analytic04, -one_of(regressors))

analytic04 <- analytic04 %>%
        mutate(useWageLevel = exp(useWage)) %>%
        mutate(useWageLevel = as.numeric(ifelse(lfStat > 1 , 0. , useWageLevel)))

# Create quarterly residual wage and quarterly lfStat
analytic04 <- analytic04 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(useWageQtr = sum(useWageLevel)) %>%
        mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
        mutate(useWageQtr = log(useWageQtr))

# Calculate last residual wage, fill down
# analytic04 <- fillDownWage(analytic04)
analytic04 <- fillUpWage(analytic04)

# Calculate residual wage change
analytic04 <- calculateWageChange(analytic04)

# job-job changes, including the separations. If separated,  pct change is -1
analytic04 <- analytic04 %>%
        mutate(wageChange_wU = as.numeric(ifelse(lfStat == 1 ,wageChange , -1.))) %>%
        mutate(wageChangeQtr_wU = as.numeric(ifelse(lfStatQtr == 1 ,wageChangeQtr , -1.)))

# Save data, remove from environment
if(useRegResid) {
        setwd("./RegResid")
} else {
        setwd("./Raw")
}

saveRDS(analytic04, "analytic04.RData")
rm(list = c("processed04", "analytic04"))

setwd("../../")

# 2008 Panel --------------------------------------------------------------

if(useSoc2d) {
        setwd("./soc2d")
}else {
        setwd("./occ")
}

processed08 <- readRDS("processed08.RData")

# Generate regressor variables
analytic08 <- genRegressors(processed08)

# Run regression within each year, remove regressors
analytic08 <- calculateUseWage(analytic08, avg1996)
analytic08 <- calculateOccWage(analytic08,0.)
analytic08 <- select(analytic08, -one_of(regressors))

analytic08 <- analytic08 %>%
        mutate(useWageLevel = exp(useWage)) %>%
        mutate(useWageLevel = as.numeric(ifelse(lfStat > 1 , 0. , useWageLevel)))

# Create quarterly residual wage and quarterly lfStat
analytic08 <- analytic08 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(useWageQtr = sum(useWageLevel)) %>%
        mutate(lfStatQtr = as.integer(ifelse(useWageQtr > 0, 1, 2))) %>%
        mutate(useWageQtr = log(useWageQtr))

# Calculate last residual wage, fill down
# analytic08 <- fillDownWage(analytic08)
analytic08 <- fillUpWage(analytic08)

# Calculate residual wage change
analytic08 <- calculateWageChange(analytic08)

# job-job changes, including the separations. If separated,  pct change is -1
analytic08 <- analytic08 %>%
        mutate(wageChange_wU = as.numeric(ifelse(lfStat == 1 ,wageChange , -1.))) %>%
        mutate(wageChangeQtr_wU = as.numeric(ifelse(lfStatQtr == 1 ,wageChangeQtr , -1.)))

# Save data, remove from environment
if(useRegResid) {
        setwd("./RegResid")
} else {
        setwd("./Raw")
}

saveRDS(analytic08, "analytic08.RData")
rm(list = c("processed08", "analytic08"))

setwd("../../")
