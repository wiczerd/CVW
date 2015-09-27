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

wageChangesFull <- readRDS("wageChanges.RData")

# Only count EEs and balanced EUs and UEs
wageChangesFull <- wageChangesFull %>%
	group_by(id) %>%
	arrange(date) %>%
	mutate(balancedEU = EU & lead(UE)) %>%
	mutate(balancedUE = UE & lag(EU)) %>%
	filter(EE | balancedEU | balancedUE) %>%
	ungroup

# Change switchedOcc to TRUE for UE if corresponding EU is TRUE
wageChangesFull <- wageChangesFull %>%
	group_by(id) %>%
	arrange(date) %>%
	mutate(switchedOcc = ifelse(UE & lag(switchedOcc), TRUE, switchedOcc)) %>%
	ungroup

# Set HSCol and Young to max over panel to ensure balance
wageChangesFull <- wageChangesFull %>%
	group_by(id) %>%
	mutate(HSCol = max(HSCol)) %>%
	mutate(Young = max(Young))

# create new person weights
wageChangesFull <- wageChangesFull %>%
	group_by(id) %>%
	mutate(newwt = mean(wpfinwgt, na.rm = TRUE)) %>%
	ungroup

# Fallick/Fleischman numbers
EEprob <- 1.3
# halfway between EU and UE prob
UEprob <- 0.85
ratio <- UEprob/EEprob

# re-weight for total distribution
UEweight <- with(wageChangesFull, sum(newwt[UE & !is.na(wageChange)], na.rm = TRUE))
EUweight <- with(wageChangesFull, sum(newwt[EU & !is.na(wageChange)], na.rm = TRUE))
EEweight <- with(wageChangesFull, sum(newwt[EE & !is.na(wageChange)], na.rm = TRUE))

# Force weights to match ratio
multiple <- ratio*EEweight/EUweight
wageChangesFull <- wageChangesFull %>%
	mutate(reweight = ifelse(EU | UE, newwt*multiple*0.5, newwt))

# Scale weights back to original total
weightMultiple <- with(wageChangesFull, sum(reweight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE))
wageChangesFull <- wageChangesFull %>%
	mutate(reweight = reweight/weightMultiple)

# check weights
UEweight <- with(wageChangesFull, sum(reweight[UE & !is.na(wageChange)], na.rm = TRUE))
EUweight <- with(wageChangesFull, sum(reweight[EU & !is.na(wageChange)], na.rm = TRUE))
EEweight <- with(wageChangesFull, sum(reweight[EE & !is.na(wageChange)], na.rm = TRUE))

# new weight ratio == Fallick Fleischmann UE/EE ratio
round((UEweight + EUweight)/EEweight, 10) == round(ratio, 10)
# new weights sum to original weights
with(wageChangesFull, sum(wpfinwgt, na.rm = TRUE) == sum(reweight, na.rm = TRUE))
# EU and UE balance
with(wageChangesFull, sum(reweight[UE & switchedOcc], na.rm = TRUE)) == 
	with(wageChangesFull, sum(reweight[EU & switchedOcc], na.rm = TRUE))

# Dictionary --------------------------------------------------------------

# newwt: new weight with individual's average weight over panel
# reweight: weights that increase representation of EU and UE's to handle truncation
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


# Loop --------------------------------------------------------------------

dataTable <- c()

# Young/Old loop
for(youngIndic in c(TRUE, FALSE)){
	
	wageChangesYoung <- subset(wageChangesFull, Young == youngIndic)
	
	for(switchIndic in c(TRUE, FALSE)) {
		
		wageChanges <- subset(wageChangesYoung, switchedOcc == switchIndic)
		
		# All job changes ---------------------------------------------------------
		
		# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
		attach(wageChanges)
		
		qtls <- c(0.1,0.25,0.5,0.75,0.9)
		# Central tendency
		full.Mean <- wtd.mean(wageChange, reweight)
		full.Qs <- wtd.quantile(wageChange, reweight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wageChange, reweight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( reweight*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(reweight,na.rm=T)
		
		fullRow <- c("full", NA, switchIndic, youngIndic, round(c(full.Mean, full.Med, full.Var, full.IQR, 
			     full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes -------------------------------------------------------------
		
		# Central tendency
		EUE.Mean <- wtd.mean(wageChange_EUE, wpfinwgt)
		EUE.Qs <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wageChange_EUE, wpfinwgt)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(wpfinwgt*(wageChange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(wpfinwgt ,na.rm=T)
		
		EUERow <- c("EUE", NA, switchIndic, youngIndic, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
			    EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution --------------------------------------------------------
		
		# Central tendency
		EE.Mean <- wtd.mean(wageChange[EE], wpfinwgt[EE])
		EE.Qs <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wageChange[EE], wpfinwgt[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(wpfinwgt[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(wpfinwgt[EE] ,na.rm=T)
		
		EERow <- c("EE", NA, switchIndic, youngIndic, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
			   EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wageChanges)
	}
	
}

# HS/College loop
for(HSColIndic in c(0, 1, 2)){
	
	wageChangesHSCol <- subset(wageChangesFull, HSCol == HSColIndic)
	
	for(switchIndic in c(TRUE, FALSE)) {
		
		wageChanges <- subset(wageChangesHSCol, switchedOcc == switchIndic)
		
		# All job changes ---------------------------------------------------------
		
		# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
		attach(wageChanges)
		
		qtls <- c(0.1,0.25,0.5,0.75,0.9)
		# Central tendency
		full.Mean <- wtd.mean(wageChange, reweight)
		full.Qs <- wtd.quantile(wageChange, reweight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wageChange, reweight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( reweight*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(reweight,na.rm=T)
		
		fullRow <- c("full", HSColIndic, switchIndic, NA, round(c(full.Mean, full.Med, full.Var, full.IQR, 
			     full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes -------------------------------------------------------------
		
		# Central tendency
		EUE.Mean <- wtd.mean(wageChange_EUE, wpfinwgt)
		EUE.Qs <- wtd.quantile(wageChange_EUE, wpfinwgt, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wageChange_EUE, wpfinwgt)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(wpfinwgt*(wageChange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(wpfinwgt ,na.rm=T)
		
		EUERow <- c("EUE", HSColIndic, switchIndic, NA, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
			    EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution --------------------------------------------------------
		
		# Central tendency
		EE.Mean <- wtd.mean(wageChange[EE], wpfinwgt[EE])
		EE.Qs <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wageChange[EE], wpfinwgt[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(wpfinwgt[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(wpfinwgt[EE] ,na.rm=T)
		
		EERow <- c("EE", HSColIndic, switchIndic, NA, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
			   EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wageChanges)
	}
	
}

# Results table formatting
dataTable <- data.frame(dataTable)
names(dataTable) <- c("Distribution", "HSCol", "Switched", "Young", "Mean", "Median", "Variance", "IQR", "90-10", "Kelly", "Pearson")

for(varNo in 5:11) {
	dataTable[,varNo] <- as.numeric(as.character(dataTable[,varNo]))
}

YoungOld <- subset(dataTable, !is.na(Young))
YoungOld$HSCol <- NULL

HSCol <- subset(dataTable, !is.na(HSCol))
HSCol$Young <- NULL

YoungOld <- reshape(YoungOld, timevar = "Distrbution", idvar = c("Switched", "Young"), direction = "wide")
HSCol <- reshape(HSCol, timevar = "Distribution", idvar = c("Switched", "HSCol"), direction = "wide")