library(dplyr)
library(ggplot2)
library(Hmisc)
library(zoo)
library(xtable)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching Lab/Data/")
setwd('~/workspace/CVW/R')

useSoc2d <- TRUE
useRegResid <- TRUE

if(useSoc2d) {
	setwd("./soc2d")
} else {
	setwd("./occ")
}

if(useRegResid) {
	setwd("./regResid")
} else {
	setwd("./raw")
}

wageChanges <- readRDS("wageChanges.RData")

attach(wageChanges)

# Remove EU switches that don't involve a later UE switch
wageChanges <- wageChanges %>%
	group_by(id) %>%
	mutate(beenEmployed = ifelse(lfStat == 1, TRUE, NA)) %>%
	mutate(beenEmployed = na.locf(beenEmployed, na.rm = FALSE))

wageChanges <- wageChanges %>%
	mutate(wageChange = ifelse(EU & is.na(wageChange_EUE), NA, wageChange)) %>%
	mutate(wageChange = ifelse(UE & is.na(beenEmployed), NA, wageChange))

# Dictionary --------------------------------------------------------------

# wageChange: wage changes for all job switchers EE, EU, and UE. NA if no switch. Changes
# 		from E to U are 0 - wage before switch. Changes from U to E are next reported
# 		wage - 0. Does not measure change between pre-unemployment wage and post-unemployment
# 		wage.
# 
# wageChange_EUE: wage changes for switches that involve a spell of unemployment only. NA for all others.
# 		At the time EU switch occurs, change is next reported wage - previous month's wage, so
# 		periods of 0 wage during unemployment are not counted.
# 
# wageChange_stayer: wage changes for people that don't switch jobs only. NA for all others.
# 
# wageChange_all: wage changes for both job stayers and movers. If moved, uses wageChange.
# 		If not moved, used wageChange_stayer.


# All job changes ---------------------------------------------------------

# Central tendency
fullMean <- wtd.mean(wageChange, wpfinwgt)
fullMed <- wtd.quantile(wageChange, wpfinwgt, probs = 0.5)
p50 <- fullMed

# Dispersion
fullVar <- wtd.var(wageChange, wpfinwgt)
fullIQR <- wtd.quantile(wageChange, wpfinwgt, probs = 0.75) - wtd.quantile(wageChange, wpfinwgt, probs = 0.25)
p90 <- wtd.quantile(wageChange, wpfinwgt, probs = 0.90)
p10 <- wtd.quantile(wageChange, wpfinwgt, probs = 0.10)
full9010 <- p90 - p10

# Skewness
fullKelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)

# EUE changes -------------------------------------------------------------

# Central tendency
EUEMean <- wtd.mean(wageChange_EUE, wpfinwgt)
EUEMed <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = 0.5)
p50 <- EUEMed

# Dispersion
EUEVar <- wtd.var(wageChange_EUE, wpfinwgt)
EUEIQR <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = 0.75) - wtd.quantile(wageChange_EUE, wpfinwgt, probs = 0.25)
p90 <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = 0.90)
p10 <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = 0.10)
EUE9010 <- p90 - p10

# Skewness
EUEKelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)

# EE distribution --------------------------------------------------------

# Central tendency
EEMean <- wtd.mean(wageChange[EE], wpfinwgt[EE])
EEMed <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = 0.5)
p50 <- EEMed

# Dispersion
EEVar <- wtd.var(wageChange[EE], wpfinwgt[EE])
EEIQR <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = 0.75) - wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = 0.25)
p90 <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = 0.90)
p10 <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = 0.10)
EE9010 <- p90 - p10

# Skewness
EEKelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)


# Put all of this into 1 table


ggplot(wageChanges, aes(wageChange)) +
	geom_density()


