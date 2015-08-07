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
# Functions ---------------------------------------------------------------
# Create function to calculate residuals
# add a constant (previously was constant from regression, now is 1996 avg)
calculateResiduals <- function(df, const = 0) {
        if(useRegResid){
			# group data by year - this doesn't make sense if we're doing business cycles
	        df <- ungroup(df)
	        # regression within each year
	        model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
	                             female + black + hispanic + factor(soc2d), data = df,
	                     na.action = na.exclude, weights = wpfinwgt)
	        # calculate residuals
	        resid <- residuals(model) + const
	        # append residuals to df
	        result <- data.frame(df, resid)
	        # ungroup and remove regressions
	        result <- result %>% 
	        	ungroup
        }else{
        	result <- df %>%
        		mutate(resid = logEarnm) 
        }
        return(result)
}

calculateOccWage <- function(df, const =0 ){
	df <- group_by(df, occ)
	model <- lm( logEarnm ~ experience + I(experience^2) + factor(educ) 
				 + female + black + hispanic, 
				 data= df, na.action = na.exclude, weights = wpfinwgt)
	occWage <- fitted(model) + const
	result <- data.frame(df,occWage)
	result <- result %>% 
		ungroup
	
}


# Create function to calculate last residual observed wage
fillDownResidual <- function(df) {
        result <- df %>%
                group_by(id) %>%
                arrange(id, date) %>%
                mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0 & lag(resid) != NA , lag(resid), NA))) %>%
                mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE)) %>%
                mutate(lastResidWage_q = as.numeric(ifelse(switchedJob & job != 0, lag(resid_q, 3), NA))) %>%
                mutate(lastResidWage_q = na.locf(lastResidWage_q, na.rm = FALSE))
        result <- result %>%
        	group_by(id) %>%
        	arrange(id, date) %>%
        	mutate(lastOcc = as.integer(ifelse(switchedJob & job != 0, lag(occ), NA))) %>%
        	mutate(lastOcc = na.locf(lastOcc, na.rm = FALSE)) %>%
        	mutate(lastOccWage = as.numeric(ifelse(switchedJob & job != 0, lag(occWage), NA))) %>%
        	mutate(lastOccWage = na.locf(lastOccWage, na.rm = FALSE)) 
        return(result)
}

# Create function to generate regressor variables and inflation adjusts
genRegressors <- function(df) {
  # import PCE data
  df <- left_join(df,PCE, by="date")
  result <- df %>%
    mutate(nomwage= wage) %>%
  	mutate(wage = wage/PCEPI*100) %>%
    mutate(logWage = log(wage)) %>%
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
  if(useSoc2d){
  	result <- mutate(result,soc2d = occ)
  }
  return(result)
}

regressors <- c("age", "educ", "female", "race", "yearsSchool",
                "experience", "black", "hispanic", "year", 
                "earnm", "logEarnm")

# 1996 Panel --------------------------------------------------------------

if(useSoc2d) {
        processed96 <- readRDS("./Data/processed96soc2d.RData")
}else {
        processed96 <- readRDS("./Data/processed96.RData")
}

# Generate regressor variables
analytic96 <- genRegressors(processed96)

# Find average log wage for 1996 panel
avg1996 <- weighted.mean(analytic96$logEarnm, analytic96$wpfinwgt, na.rm = TRUE)
# Run regression within each year, remove regressors
analytic96 <- calculateResiduals(analytic96, avg1996)
analytic96 <- calculateOccWage(analytic96,0.)
analytic96 <- select(analytic96, -one_of(regressors))
analytic96 <- analytic96 %>%
        mutate(resid_lev = exp(resid)) %>%
        mutate(resid_lev = as.numeric(ifelse(lfStat > 1 , 0. , resid_lev)))

# Create quarterly residual wage and quarterly lfStat
analytic96 <- analytic96 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(resid_q = sum(resid_lev)) %>%
        mutate(lfStat_q = as.integer(ifelse(resid_q > 0, 1, 2))) %>%
        mutate(resid_q = log(resid_q))

# Calculate last residual wage, fill down
analytic96 <- fillDownResidual(analytic96)

# Calculate residual wage change
analytic96 <- analytic96 %>%
        mutate(residWageChange = lead(resid) - lastResidWage) %>%
        mutate(residWageChange_q = resid_q - lastResidWage_q) %>%
		mutate(occWageChange = lead(occWage) - lastOccWage) %>%
		mutate(residWageChange_stayer = as.numeric(ifelse(!switchedJob, lead(resid) - resid,NA) ))

# job-job changes, including the separations. If separated,  pct change is -1
analytic96 <- analytic96 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
		mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
if(useSoc2d & useRegResid) {
	saveRDS(analytic96, "./Data/analytic96soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	saveRDS(analytic96, "./Data/analytic96soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	saveRDS(analytic96, "./Data/analytic96Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	saveRDS(analytic96, "./Data/analytic96Raw.RData")   
} else{
	saveRDS(analytic96, "./Data/analytic96.RData")   
}
rm(list = c("processed96", "analytic96"))

# 2001 Panel --------------------------------------------------------------

if(useSoc2d) {
        processed01 <- readRDS("./Data/processed01soc2d.RData")
} else {
        processed01 <- readRDS("./Data/processed01.RData")
}

# Generate regressor variables
analytic01 <- genRegressors(processed01)

# Run regression within each year, remove regressors
analytic01 <- calculateResiduals(analytic01, avg1996)
analytic01 <- calculateOccWage(analytic01,0.)
analytic01 <- select(analytic01, -one_of(regressors))

analytic01 <- analytic01 %>%
        mutate(resid_lev = exp(resid)) %>%
		mutate(resid_lev = as.numeric(ifelse(lfStat > 1 , 0. , resid_lev)))

# Create quarterly residual wage and quarterly lfStat
analytic01 <- analytic01 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(resid_q = sum(resid_lev)) %>%
        mutate(lfStat_q = as.integer(ifelse(resid_q > 0, 1, 2))) %>%
        mutate(resid_q = log(resid_q))

# Calculate last residual wage, fill down
analytic01 <- fillDownResidual(analytic01)

# Calculate residual wage change
analytic01 <- analytic01 %>%
        mutate(residWageChange = lead(resid) - lastResidWage) %>%
        mutate(residWageChange_q = resid_q - lastResidWage_q) %>%
		mutate(occWageChange = lead(occWage) - lastOccWage) %>%
		mutate(residWageChange_stayer = as.numeric(ifelse(!switchedJob, lead(resid) - resid,NA) ))


# job-job changes, including the separations. If separated,  pct change is -1
analytic01 <- analytic01 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
if(useSoc2d & useRegResid) {
	saveRDS(analytic01, "./Data/analytic01soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	saveRDS(analytic01, "./Data/analytic01soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	saveRDS(analytic01, "./Data/analytic01Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	saveRDS(analytic01, "./Data/analytic01Raw.RData")   
} else{
	saveRDS(analytic01, "./Data/analytic01.RData")   
}
rm(list = c("processed01", "analytic01"))


# 2004 Panel --------------------------------------------------------------

if(useSoc2d) {
        processed04 <- readRDS("./Data/processed04soc2d.RData")
} else {
        processed04 <- readRDS("./Data/processed04.RData")
}

# Generate regressor variables
analytic04 <- genRegressors(processed04)

# Run regression within each year, remove regressors
analytic04 <- calculateResiduals(analytic04, avg1996)
analytic04 <- calculateOccWage(analytic04,0.)
analytic04 <- select(analytic04, -one_of(regressors))

analytic04 <- analytic04 %>%
        mutate(resid_lev = exp(resid)) %>%
		mutate(resid_lev = as.numeric(ifelse(lfStat > 1 , 0. , resid_lev)))

# Create quarterly residual wage and quarterly lfStat
analytic04 <- analytic04 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(resid_q = sum(resid_lev)) %>%
        mutate(lfStat_q = as.integer(ifelse(resid_q > 0, 1, 2))) %>%
        mutate(resid_q = log(resid_q))

# Calculate last residual wage, fill down
analytic04 <- fillDownResidual(analytic04)

# Calculate residual wage change
analytic04 <- analytic04 %>%
		mutate(residWageChange = lead(resid) - lastResidWage) %>%
		mutate(residWageChange_q = resid_q - lastResidWage_q) %>%
		mutate(occWageChange = lead(occWage) - lastOccWage) %>%
		mutate(residWageChange_stayer = as.numeric(ifelse(!switchedJob, lead(resid) - resid,NA) ))

# job-job changes, including the separations. If separated,  pct change is -1
analytic04 <- analytic04 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
if(useSoc2d & useRegResid) {
	saveRDS(analytic04, "./Data/analytic04soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	saveRDS(analytic04, "./Data/analytic04soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	saveRDS(analytic04, "./Data/analytic04Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	saveRDS(analytic04, "./Data/analytic04Raw.RData")   
} else{
	saveRDS(analytic04, "./Data/analytic04.RData")   
}
rm(list = c("processed04", "analytic04"))

# 2008 Panel --------------------------------------------------------------

if(useSoc2d) {
        processed08 <- readRDS("./Data/processed08soc2d.RData")
} else {
        processed08 <- readRDS("./Data/processed08.RData")
}

# Generate regressor variables
analytic08 <- genRegressors(processed08)

# Run regression within each year, remove regressors
analytic08 <- calculateResiduals(analytic08, avg1996)
analytic08 <- calculateOccWage(analytic08,0.)
analytic08 <- select(analytic08, -one_of(regressors))

analytic08 <- analytic08 %>%
		mutate(resid_lev = exp(resid)) %>%
		mutate(resid_lev = as.numeric(ifelse(lfStat > 1 , 0. , resid_lev)))

# Create quarterly residual wage and quarterly lfStat
analytic08 <- analytic08 %>%
        filter(!is.na(lfStat)) %>%
        group_by(id, qtrdate) %>%
        mutate(resid_q = sum(resid_lev)) %>%
        mutate(lfStat_q = as.integer(ifelse(resid_q > 0, 1, 2))) %>%
        mutate(resid_q = log(resid_q))

# Calculate last residual wage, fill down
analytic08 <- fillDownResidual(analytic08)

# Calculate residual wage change
analytic08 <- analytic08 %>%
        mutate(residWageChange = lead(resid) - lastResidWage) %>%
        mutate(residWageChange_q = resid_q - lastResidWage_q) %>%
		mutate(occWageChange = lead(occWage) - lastOccWage) %>%
		mutate(residWageChange_stayer = as.numeric(ifelse(!switchedJob, lead(resid) - resid,NA) ))

# job-job changes, including the separations. If separated,  pct change is -1
analytic08 <- analytic08 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
if(useSoc2d & useRegResid) {
        saveRDS(analytic08, "./Data/analytic08soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
		saveRDS(analytic08, "./Data/analytic08soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	saveRDS(analytic08, "./Data/analytic08Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	saveRDS(analytic08, "./Data/analytic08Raw.RData")   
} else{
        saveRDS(analytic08, "./Data/analytic08.RData")   
}
rm(list = c("processed08", "analytic08"))
