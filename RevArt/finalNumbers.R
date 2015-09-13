library(dplyr)
library(ggplot2)
library(Hmisc)
library(zoo)
library(xtable)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching Lab/Data/")
setwd('~/workspace/CVW/R/Data')

useSoc2d <- TRUE
useRegResid <- TRUE

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

wageChanges <- readRDS("wageChanges.RData")

wageChanges$wtchng <- ifelse(wageChanges$EU | wageChanges$UE, wageChanges$wpfinwgt/2,wageChanges$wpfinwgt)

# Make sure EU and UE switches balance
wageChanges <- wageChanges %>%
	group_by(id) %>%
	arrange(id,date) %>%
	mutate(beenEU = ifelse(EU, -1, 0)) %>%
	mutate(willUE = ifelse(UE, 1 , 0)) %>%
	mutate(matchEUUE = cumsum(beenEU + willUE)) %>%
	arrange(id,desc(date)) %>%
	mutate(matchUEEU = cumsum(beenEU + willUE)) %>%
	ungroup

wageChanges <- wageChanges %>%
	mutate(wageChange = ifelse(matchEUUE != 0 & UE , as.numeric(NA), wageChange)) %>%
	mutate(wageChange = ifelse(matchUEEU != 0 & EU , as.numeric(NA), wageChange))

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


attach(wageChanges)
# Central tendency
fullMean <- wtd.mean(wageChange, wtchng)
fullMed <- wtd.quantile(wageChange, wtchng, probs = 0.5)
p50 <- fullMed

# Dispersion
fullVar <- wtd.var(wageChange, wtchng)
fullIQR <- wtd.quantile(wageChange, wtchng, probs = 0.75) - wtd.quantile(wageChange, wtchng, probs = 0.25)
p90 <- wtd.quantile(wageChange, wtchng, probs = 0.90)
p10 <- wtd.quantile(wageChange, wtchng, probs = 0.10)
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


# Put all of the central tendency into 1 table
cent_tend <- c(fullMean)

ggplot(wageChanges, aes(wageChange)) +
	geom_density()


