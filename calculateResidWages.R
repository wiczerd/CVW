# April 27, 2015
# Modify calculateResidWages to create quarterly total earnings and use those
# to calculate wage changes.
# Save analytic data in ./Data/ directory.
# Precondition: processData.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)

setwd("~/workspace/CVW/R")

# Functions ---------------------------------------------------------------
# Create function to calculate residuals
calculateResiduals <- function(df) {
        model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
                            female + black + hispanic + factor(soc2d), data = df,
                    na.action = na.exclude, weights = wpfinwgt)
        resid <- residuals(model)
        resid <- resid + coef(model)[1]
        data.frame(df, resid)
}

# Create function calculate last residual observed wage
fillDownResidual <- function(df) {
        df %>%
                group_by(id) %>%
                arrange(id, date) %>%
                mutate(lastResidWage = as.numeric(ifelse(switchedJob & job != 0, lag(resid), NA))) %>%
                mutate(lastResidWage = na.locf(lastResidWage, na.rm = FALSE)) %>%
                mutate(lastResidWage_q = as.numeric(ifelse(switchedJob & job != 0, lag(resid_q, 3), NA))) %>%
                mutate(lastResidWage_q = na.locf(lastResidWage_q, na.rm = FALSE)) 
}

# Create function to generate regressor variables
genRegressors <- function(df) {
        df %>%
                mutate(logEarnm = log(earnm)) %>%
                mutate(yearsSchool = as.integer(ifelse(educ == 1, 9, NA)),
                       yearsSchool = as.integer(ifelse(educ == 2, 12, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 3, 14, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 4, 16, yearsSchool)),
                       yearsSchool = as.integer(ifelse(educ == 5, 18, yearsSchool))) %>%
                mutate(experience = age - yearsSchool,
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
        calculateResiduals(.) %>%
        ungroup %>%
        select(-one_of(regressors))

analytic96 <- analytic96 %>%
        mutate(resid_lev = exp(resid)) %>%
        mutate(resid_lev = ifelse(lfStat > 1 , 0. , resid_lev))

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
        mutate(residWageChange_q = resid_q - lastResidWage_q)

# job-job changes, including the separations. If separated,  pct change is -1
analytic96 <- analytic96 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

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
        calculateResiduals(.) %>%
        ungroup %>%
        select(-one_of(regressors))

analytic01 <- analytic01 %>%
        mutate(resid_lev = exp(resid)) %>%
        mutate(resid_lev = ifelse(lfStat > 1 , 0. , resid_lev))

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
        mutate(residWageChange_q = resid_q - lastResidWage_q)

# job-job changes, including the separations. If separated,  pct change is -1
analytic01 <- analytic01 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
saveRDS(analytic01, "./Data/analytic01.RData")
rm(list = c("processed01", "analytic01"))


# 2004 Panel --------------------------------------------------------------
processed04 <- readRDS("./Data/processed04.RData")

# Generate regressor variables
analytic04 <- genRegressors(processed04)

# Run regression within each year, remove regressors
analytic04 <- analytic04 %>%
        group_by(year) %>%
        calculateResiduals(.) %>%
        ungroup %>%
        select(-one_of(regressors))

analytic04 <- analytic04 %>%
        mutate(resid_lev = exp(resid)) %>%
        mutate(resid_lev = ifelse(lfStat > 1 , 0. , resid_lev))

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
        mutate(residWageChange_q = resid_q - lastResidWage_q)

# job-job changes, including the separations. If separated,  pct change is -1
analytic04 <- analytic04 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
saveRDS(analytic04, "./Data/analytic04.RData")
rm(list = c("processed04", "analytic04"))

# 2008 Panel --------------------------------------------------------------
processed08 <- readRDS("./Data/processed08.RData")

# Generate regressor variables
analytic08 <- genRegressors(processed08)

# Run regression within each year, remove regressors
analytic08 <- analytic08 %>%
        group_by(year) %>%
        calculateResiduals(.) %>%
        ungroup %>%
        select(-one_of(regressors))

analytic08 <- analytic08 %>%
        mutate(resid_lev = exp(resid)) %>%
        mutate(resid_lev = ifelse(lfStat > 1 , 0. , resid_lev))

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
        mutate(residWageChange_q = resid_q - lastResidWage_q)

# job-job changes, including the separations. If separated,  pct change is -1
analytic08 <- analytic08 %>%
        mutate(residWageChange_wU = as.numeric(ifelse(lfStat == 1 ,residWageChange , -1.))) %>%
        mutate(residWageChange_q_wU = as.numeric(ifelse(lfStat_q == 1 ,residWageChange_q , -1.)))

# Save data, remove from environment
saveRDS(analytic08, "./Data/analytic08.RData")
rm(list = c("processed08", "analytic08"))
