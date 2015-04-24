# April 2, 2015
# Calculate residual log wage and mean residual log wage per 
# person per occupation.
# Save analytic data in ./Data/ directory.
# Precondition: processData.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)

setwd("~/workspace/CVW/R")

# Create residuals function to apply with do()
calculateResiduals <- function(df) {
	model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
					female + black + hispanic + factor(soc2d), data = df,
				na.action = na.exclude, weights = wpfinwgt)
	resid <- residuals(model)
	data.frame(df, resid)
}

# Create function to generate regressor variables
genRegressors <- function(df) {
	mutate(df, logEarnm = log(earnm),
		   yearsSchool = as.integer(ifelse(educ == 1, 9, NA)),
		   yearsSchool = as.integer(ifelse(educ == 2, 12, yearsSchool)),
		   yearsSchool = as.integer(ifelse(educ == 3, 14, yearsSchool)),
		   yearsSchool = as.integer(ifelse(educ == 4, 16, yearsSchool)),
		   yearsSchool = as.integer(ifelse(educ == 5, 18, yearsSchool)),
		   experience = age - yearsSchool,
		   black = (race == 2),
		   hispanic = (race == 3),
		   year = as.numeric(format(date, "%Y")))
}

regressors <- c("age", "educ", "female", "race", "yearsSchool",
				"experience", "black", "hispanic", "year", 
				"earnm", "logEarnm")

# 1996 Panel --------------------------------------------------------------
processed96 <- readRDS("./Data/processed96.RData")

# Generate regressor variables
analytic96 <- genRegressors(processed96)

# Run regression within each year, remove regressors
analytic96 <- analytic96 %>%
	group_by(year) %>%
	do(calculateResiduals(.)) %>%
	ungroup %>%
	select(-one_of(regressors))

# Calculate last residual wage and find residual wage change
analytic96 <- analytic96 %>%
	filter(!is.na(lfStat)) %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0, lag(resid), NA))) %>%
	mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE)) 
analytic96 <- analytic96 %>%
	mutate(residWageChange = lead(resid) - lastResidWage) 
# job-job changes, including the separations. If separated,  pct change is -1
analytic96 <- analytic96 %>%
	mutate(residWageChange_wU = ifelse(lfStat == 1 ,residWageChange , -1.) )
# all changes, including the separations and continuing guys
analytic96 <- analytic96 %>%
	mutate(residWageChange_wA = (lead(resid)) - resid ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lead(lfStat) == 2 | lead(lfStat) == 3, -1. , residWageChange_wA)) ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lfStat == 1, residWageChange_wA,NA ) ) )

# Save data, remove from environment
saveRDS(analytic96, "./Data/analytic96.RData")
rm(list = c("processed96", "analytic96"))

# 2001 Panel --------------------------------------------------------------
processed01 <- readRDS("./Data/processed01.RData")

# Generate regressor variables
analytic01 <- genRegressors(processed01)

# Run regression within each year, remove regressors
analytic01 <- analytic01 %>%
	group_by(year) %>%
	do(calculateResiduals(.)) %>%
	ungroup %>%
	select(-one_of(regressors))

# Calculate mean residual wage per person per occupation
analytic01 <- analytic01 %>%
	filter(!is.na(lfStat)) %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0, lag(resid), NA))) %>%
	mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE))
analytic01 <- analytic01 %>%
	mutate(residWageChange = lead(resid) - lastResidWage) 
analytic01 <- analytic01 %>%
	mutate(residWageChange_wU = ifelse(lfStat == 1 ,residWageChange , -1.) )
analytic01 <- analytic01 %>%
	mutate(residWageChange_wA = (lead(resid)) - resid ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lead(lfStat) == 2 | lead(lfStat) == 3, -1. , residWageChange_wA)) ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lfStat == 1, residWageChange_wA,NA ) ) )

# Save data, remove from environment
saveRDS(analytic01, "./Data/analytic01.RData")
rm(list = c("processed01", "analytic01"))

# 2004 Panel --------------------------------------------------------------
processed04 <- readRDS("./Data/processed04.RData")

# Generate regressor variables
analytic04 <- genRegressors(processed04)

# Run regression within each year, remove regressors after
analytic04 <- analytic04 %>%
	group_by(year) %>%
	do(calculateResiduals(.)) %>%
	ungroup %>%
	select(-one_of(regressors))

# Calculate mean residual wage per person per occupation
analytic04 <- analytic04 %>%
	filter(!is.na(lfStat)) %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0, lag(resid), NA))) %>%
	mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE)) 
analytic04 <- analytic04 %>%
	mutate(residWageChange = lead(resid) - lastResidWage)
analytic04 <- analytic04 %>%
	mutate(residWageChange_wU = ifelse(lfStat == 1 ,residWageChange , -1.) )
analytic04 <- analytic04 %>%
	mutate(residWageChange_wA = (lead(resid)) - resid ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lead(lfStat) == 2 | lead(lfStat) == 3, -1. , residWageChange_wA)) ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lfStat == 1, residWageChange_wA,NA ) ) )

# Save data, remove from environment
saveRDS(analytic04, "./Data/analytic04.RData")
rm(list = c("processed04", "analytic04"))

# 2008 Panel --------------------------------------------------------------
processed08 <- readRDS("./Data/processed08.RData")

# Generate regressor variables
analytic08 <- genRegressors(processed08)

# Run regression within each year, remove regressors after
analytic08 <- analytic08 %>%
	group_by(year) %>%
	filter(!is.na(year)) %>%
	do(calculateResiduals(.)) %>%
	ungroup %>%
	select(-one_of(regressors))

# Calculate mean residual wage per person per occupation
analytic08 <- analytic08 %>%
	filter(!is.na(lfStat)) %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0, lag(resid), NA))) %>%
	mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE))
analytic08 <- analytic08 %>%
	mutate(residWageChange = lead(resid) - lastResidWage)
analytic08 <- analytic08 %>%
	mutate(residWageChange_wU = ifelse(lfStat == 1 ,residWageChange , -1.) )
analytic08 <- analytic08 %>%
	mutate(residWageChange_wA = (lead(resid)) - resid ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lead(lfStat) == 2 | lead(lfStat) == 3, -1. , residWageChange_wA)) ) %>%
	mutate(residWageChange_wA = as.numeric(ifelse( lfStat == 1, residWageChange_wA,NA ) ) )

# Save data, remove from environment
saveRDS(analytic08, "./Data/analytic08.RData")
rm(list = c("processed08", "analytic08"))

