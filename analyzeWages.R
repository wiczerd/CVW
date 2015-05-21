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
useSoc2d <- TRUE

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

# 1996 Panel --------------------------------------------------------------

if(useSoc2d) {
        analytic96 <- readRDS("./Data/analytic96soc2d.RData")
} else {
        analytic96 <- readRDS("./Data/analytic96.RData")
}

# Calculate income percentile within two-digit SOC code
incomePercentiles <- analytic96 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE)

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

if(useSoc2d) {
        analytic01 <- readRDS("./Data/analytic01soc2d.RData")
} else {
        analytic01 <- readRDS("./Data/analytic01.RData")
}

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic01 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE) %>%
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
analytic04 <- readRDS("./Data/analytic04.RData")

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic04 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE) %>%
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

if(useSoc2d) {
        analytic08 <- readRDS("./Data/analytic08soc2d.RData")
} else {
        analytic08 <- readRDS("./Data/analytic08.RData")
}

# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic08 %>%
        filter(!is.na(soc2d)) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc, percentile, EE, UE) %>%
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

incomePercentiles <- incomePercentiles %>%
        filter(!is.na(percentile)) %>%
        group_by(percentile) %>%
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