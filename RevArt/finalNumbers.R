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

wageChangesRead <- readRDS("wageChanges.RData")
UEreadweight <- with(wageChangesRead, sum(newwt[UE & !is.na(wageChange)], na.rm = TRUE))
EUreadweight <- with(wageChangesRead, sum(newwt[EU & !is.na(wageChange)], na.rm = TRUE))
EEreadweight <- with(wageChangesRead, sum(newwt[EE & !is.na(wageChange)], na.rm = TRUE))
UEreadcount <- with(wageChangesRead,sum(UE & !is.na(wageChange)) )
EUreadcount <- with(wageChangesRead,sum(EU & !is.na(wageChange)) )
EEreadcount <- with(wageChangesRead,sum(EE & !is.na(wageChange)) )

wageChangesFull <- wageChangesRead

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


rm(wageChanges)

# Var decomp of wage changes btwn EE, EUE, stayers---------------

analytic9608<-readRDS("analytic9608.RData")

# how many EU's here?
#sum(analytic9608$EU[analytic9608$lfStat==1]*analytic9608$wpfinwgt[analytic9608$lfStat==1],na.rm=T)/sum(analytic9608$wpfinwgt[analytic9608$lfStat==1],na.rm=T)
#0.006766477
analytic9608 <- analytic9608 %>%
	group_by(id) %>%
	mutate(newwt = mean(wpfinwgt, na.rm = TRUE)) %>%
	ungroup
# balance UEs and EUs 
analytic9608 <- analytic9608 %>%
	arrange(id,date) 
analytic9608$lagEUE <- ifelse(analytic9608$EU & is.finite(analytic9608$wageChange_EUE),T,NA)
analytic9608$lagEUE <- ifelse(analytic9608$UE,F,analytic9608$lagEUE)
analytic9608$lagEUE <- na.locf(analytic9608$lagEUE, na.rm=F)
analytic9608$lagEUE <- lag(analytic9608$lagEUE)
sum(analytic9608$lagEUE & analytic9608$UE,na.rm=T)


# Fallick/Fleischman numbers
EEprob <- 1.3
# halfway between EU and UE prob
UEprob <- 0.85
ratio <- UEprob/EEprob


analytic9608$EE_fac <- ifelse(analytic9608$EE,1,0)
analytic9608$EUUE_fac <- ifelse(analytic9608$EU | analytic9608$UE | analytic9608$lfStat == 2,1,0)
analytic9608$EUE_fac <- ifelse(analytic9608$EU,1,0)
analytic9608$swEE_fac <- ifelse(analytic9608$EE & analytic9608$switchedOcc,1,0)
analytic9608$swEUUE_fac <- ifelse((analytic9608$EU|analytic9608$UE|analytic9608$lfStat == 2) & analytic9608$switchedOcc,1,0)
analytic9608$swEUE_fac <- ifelse(analytic9608$EU & analytic9608$switchedOcc,1,0)
analytic9608$stay <- ifelse(analytic9608$EE | analytic9608$EU | analytic9608$UE, 0 ,1)

wageChangeMean = wtd.mean(analytic9608$wageChange_all, analytic9608$newwt, na.rm=T)
wageChangeSS   = with(analytic9608, sum( ( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_stay = with(analytic9608, sum(  stay*( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_EE = with(analytic9608, sum(  EE_fac*( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_EUUE = with(analytic9608, sum(  EUUE_fac*( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))

tot = with(analytic9608, sum(  newwt ,na.rm=T))
pct_stay = with(analytic9608, sum(  stay*newwt ,na.rm=T))
pct_EE = with(analytic9608, sum(  EE_fac* newwt ,na.rm=T))
pct_EUE = with(analytic9608, sum(  EUUE_fac*newwt ,na.rm=T))
pct_swEE = with(analytic9608, sum(  swEE_fac* newwt ,na.rm=T))
pct_swEUE = with(analytic9608, sum(  swEUE_fac* newwt ,na.rm=T))


wageChangeSS_swEE = with(analytic9608, sum(  swEE_fac*( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_swEUUE =with(analytic9608,  sum(  swEUUE_fac*( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_nswEE = with(analytic9608, sum( (1.- swEE_fac)*( EE_fac)*
						  	( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))
wageChangeSS_nswEUUE =with(analytic9608,  sum( (1.- swEUUE_fac)* EUUE_fac*
								( wageChange_all- wageChangeMean)^2* newwt ,na.rm=T))

# not including the U stints
wageChangeEUEMean = wtd.mean(analytic9608$wageChange_EUE, analytic9608$newwt,na.rm=T)
wageChangeEUESS   = with(analytic9608, sum( (wageChange_EUE- wageChangeEUEMean)^2*newwt ,na.rm=T))
wageChangeEUESS_stay = with(analytic9608, sum(  stay*(wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_EE = with(analytic9608, sum(  EE_fac*(wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_EUE = with(analytic9608, sum( EUE_fac*(wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_swEE = with(analytic9608, sum( swEE_fac*(wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_swEUE = with(analytic9608, sum( swEUE_fac*(wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_nswEE = with(analytic9608, sum( (1.-swEE_fac)*swEE_fac*( wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))
wageChangeEUESS_nswEUE = with(analytic9608, sum( (1.-swEUE_fac)*swEUE_fac*( wageChange_EUE- wageChangeEUEMean)^2* newwt ,na.rm=T))

wCh_vardec <- rbind(c(wageChangeEUESS_stay/wageChangeEUESS,wageChangeEUESS_EE/wageChangeEUESS,wageChangeEUESS_EUE /wageChangeEUESS, NA),
					c(wageChangeSS_stay/wageChangeSS,wageChangeSS_EE/wageChangeSS,NA, wageChangeSS_EUUE /wageChangeSS),
					c(pct_stay/tot,pct_EE/tot,pct_EUE/tot, pct_EUE/tot))

wChsw_vardec <-rbind(c(wageChangeEUESS_swEE/wageChangeEUESS_EE,wageChangeEUESS_swEUE/wageChangeEUESS_EUE,,wageChangeEUESS_nswEUE/wageChangeEUESS_EUE),
					 c(pct_swEE/pct_EE,1.-pct_swEE/pct_EE,pct_swEUE/pct_EUE, 1-pct_swEUE/pct_EUE))
