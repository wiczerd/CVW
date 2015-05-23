# April 2, 2015
# Calculate wage percentiles, probability of switching by wage
# percentile, and mean income changes after switching
# Plot results, save to ./Figures/ directory.
# Precondition: calculateResidualWages.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)
library(ggplot2)

setwd("~/workspace/CVW/R")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- F
useRegResid <- F

# Function for weighted moving average
# x, wts, and new.name are strings
movingAverage <- function(df, x, wts, new.name) {
        sumProduct <- rollsum(df[x] * df[wts], 13, align = "center", fill = NA)
        sumWts <- rollsum(df[wts], 13, align = "center", fill = NA)
        movingAvg <- sumProduct/sumWts
        df[new.name] <- movingAvg
        return(df)
}

# Create cumuluative distribution function
calculateCumul <- function(df) {
        cumulFunc <- ecdf(df$resid)
        percentile <- round(cumulFunc(df$resid)*100)
        data.frame(df, percentile)
}

rec_dates   <- as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2009-06-01"))

# 1996 Panel --------------------------------------------------------------

if(useSoc2d & useRegResid) {
	analytic96 <- readRDS( "./Data/analytic96soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic96 <- readRDS("./Data/analytic96soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic96 <- readRDS("./Data/analytic96Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic96 <- readRDS("./Data/analytic96Raw.RData")   
} else{
	analytic96 <- readRDS("./Data/analytic96.RData")   
}

analytic96$Rec <- with(analytic96, (date>rec_dates[1] & date<rec_dates[2]) | 
                                            (date>rec_dates[3] & date<rec_dates[4]))

# Calculate income percentile within two-digit SOC code
incomePercentiles <- analytic96 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE, date)

# Calculate mean residual income changes by labor force flow
incomeChanges <- analytic96 %>%
        filter(switchedOcc) %>%
        group_by(date) %>%
        summarize(meanIncomeChange = weighted.mean(residWageChange[EE | UE],
                                                   wpfinwgt[EE | UE], na.rm = TRUE),
                  meanIncomeChangeEE = weighted.mean(residWageChange[EE], 
                                                     wpfinwgt[EE], na.rm = TRUE),
                  meanIncomeChangeUE = weighted.mean(residWageChange[UE], 
                                                     wpfinwgt[UE], na.rm = TRUE),
                  numPooled = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  numEE = sum(wpfinwgt[EE], na.rm = TRUE),
                  numUE = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = "1996")

# Remove from environment
rm(analytic96)

# 2001 Panel --------------------------------------------------------------


if(useSoc2d & useRegResid) {
	analytic01 <- readRDS( "./Data/analytic01soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic01 <- readRDS("./Data/analytic01soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic01 <- readRDS("./Data/analytic01Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic01 <- readRDS("./Data/analytic01Raw.RData")   
} else{
	analytic01 <- readRDS("./Data/analytic01.RData")   
}

analytic01$Rec <- with(analytic01, (date>rec_dates[1] & date<rec_dates[2]) | 
                               (date>rec_dates[3] & date<rec_dates[4]))

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic01 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE, date) %>%
        bind_rows(incomePercentiles)

# Calculate mean residual income changes by labor force flow, append
incomeChanges <- analytic01 %>%
        filter(switchedOcc) %>%
        group_by(date) %>%
        summarize(meanIncomeChange = weighted.mean(residWageChange[EE | UE],
                                                   wpfinwgt[EE | UE], na.rm = TRUE),
                  meanIncomeChangeEE = weighted.mean(residWageChange[EE], 
                                                     wpfinwgt[EE], na.rm = TRUE),
                  meanIncomeChangeUE = weighted.mean(residWageChange[UE], 
                                                     wpfinwgt[UE], na.rm = TRUE),
                  numPooled = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  numEE = sum(wpfinwgt[EE], na.rm = TRUE),
                  numUE = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = "2001") %>%
        bind_rows(incomeChanges)

# Remove from environment
rm(analytic01)

# 2004 Panel --------------------------------------------------------------

if(useSoc2d & useRegResid) {
	analytic04 <- readRDS( "./Data/analytic04soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic04 <- readRDS("./Data/analytic04soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic04 <- readRDS("./Data/analytic04Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic04 <- readRDS("./Data/analytic04Raw.RData")   
} else{
	analytic04 <- readRDS("./Data/analytic04.RData")   
}

analytic04$Rec <- with(analytic04, (date>rec_dates[1] & date<rec_dates[2]) | 
                               (date>rec_dates[3] & date<rec_dates[4]))

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic04 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE, date) %>%
        bind_rows(incomePercentiles)

# Calculate mean residual income changes by labor force flow, append
incomeChanges <- analytic04 %>%
        filter(switchedOcc) %>%
        group_by(date) %>%
        summarize(meanIncomeChange = weighted.mean(residWageChange[EE | UE],
                                                   wpfinwgt[EE | UE], na.rm = TRUE),
                  meanIncomeChangeEE = weighted.mean(residWageChange[EE], 
                                                     wpfinwgt[EE], na.rm = TRUE),
                  meanIncomeChangeUE = weighted.mean(residWageChange[UE], 
                                                     wpfinwgt[UE], na.rm = TRUE),
                  numPooled = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  numEE = sum(wpfinwgt[EE], na.rm = TRUE),
                  numUE = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = "2004") %>%
        bind_rows(incomeChanges)

# Remove from environment
rm(analytic04)

# 2008 Panel --------------------------------------------------------------

if(useSoc2d & useRegResid) {
	analytic08 <- readRDS("./Data/analytic08soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic08 <- readRDS("./Data/analytic08soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic08 <- readRDS("./Data/analytic08Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic08 <- readRDS("./Data/analytic08Raw.RData")   
} else{
	analytic08 <- readRDS("./Data/analytic08.RData")   
}

analytic08$Rec <- with(analytic08, (date>rec_dates[1] & date<rec_dates[2]) | 
                               (date>rec_dates[3] & date<rec_dates[4]))

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic08 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE, date) %>%
        bind_rows(incomePercentiles)

# Calculate mean residual income changes by labor force flow, append
incomeChanges <- analytic08 %>%
        filter(switchedOcc) %>%
        group_by(date) %>%
        summarize(meanIncomeChange = weighted.mean(residWageChange[EE | UE],
                                                   wpfinwgt[EE | UE], na.rm = TRUE),
                  meanIncomeChangeEE = weighted.mean(residWageChange[EE], 
                                                     wpfinwgt[EE], na.rm = TRUE),
                  meanIncomeChangeUE = weighted.mean(residWageChange[UE], 
                                                     wpfinwgt[UE], na.rm = TRUE),
                  numPooled = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  numEE = sum(wpfinwgt[EE], na.rm = TRUE),
                  numUE = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = "2008") %>%
        bind_rows(incomeChanges)

# Remove from environment
rm(analytic08)

# Plots -------------------------------------------------------------------

if(useSoc2d) {
        setwd("./Figures/soc2d")
} else {
        setwd("./Figures/occ")
}

rec_dates   <- as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2009-06-01"))
incomePercentiles$Rec <- ((incomePercentiles$date>rec_dates[1] & incomePercentiles$date<rec_dates[2] ) | 
                     (incomePercentiles$date>rec_dates[3] & incomePercentiles$date<rec_dates[4] ))


incomePercentiles <- incomePercentiles %>%
        filter(!is.na(percentile)) %>%
        group_by(percentile, Rec) %>%
        summarize(prSwitching = weighted.mean(switchedOcc & (UE | EE), wpfinwgt, na.rm = TRUE))

incomeChanges <- incomeChanges %>%
        # why are there NaNs?
        filter(!is.nan(meanIncomeChangeEE)) %>%
        group_by(panel) %>%
        arrange(date) %>%
        do(movingAverage(., "meanIncomeChange", "numPooled", new.name = "maPooled")) %>%
        do(movingAverage(., "meanIncomeChangeEE", "numEE", new.name = "maEE")) %>%
        do(movingAverage(., "meanIncomeChangeUE", "numUE", new.name = "maUE")) %>%
        select(date, starts_with("ma"), panel)

png("prSwitchingByIncome.png", width = 782, height = 569)
ggplot(incomePercentiles, aes(percentile, prSwitching)) + geom_point() + geom_line() +
        ggtitle("P(Switching) by Income Percentile") + 
        ylab("P(Switching)") + xlab("Income percentile within two-digit SOC occupation")
dev.off()

postscript("prSwitchingByIncomeRec.eps", width = 782, height = 569)
ggplot(incomePercentiles, aes(percentile, prSwitching, color = Rec)) + 
	geom_point() + geom_line() +
	ggtitle("P(Switching) by Income Percentile") + 
	ylab("P(Switching)") + xlab("Income percentile within occupation")
dev.off()

incomeChanges <- melt(incomeChanges, id = c("date", "panel"))
levels(incomeChanges$variable) <- c("Pooled", "Employment-to-Employment", "Unemployment-to-Employment")

png("incomeChangesMA.png", width = 782, height = 569)
ggplot(incomeChanges, aes(date, value, color = panel)) + geom_point() + geom_line() +
        facet_grid(. ~ variable) +
        ggtitle("Moving Average 6(1)6 of Residual Log Income Change by LF Flow") +
        xlab("Date") + ylab("Mean Change in Residual Log Income") +
        labs(color = "Panel")
dev.off()

setwd("../../")