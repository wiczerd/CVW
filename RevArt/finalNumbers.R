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
UEreadweight <- with(wageChangesRead, sum(wpfinwgt[UE & !is.na(wageChange)], na.rm = TRUE))
EUreadweight <- with(wageChangesRead, sum(wpfinwgt[EU & !is.na(wageChange)], na.rm = TRUE))
EEreadweight <- with(wageChangesRead, sum(wpfinwgt[EE & !is.na(wageChange)], na.rm = TRUE))
UEreadcount <- with(wageChangesRead,sum(UE & !is.na(wageChange)) )
EUreadcount <- with(wageChangesRead,sum(EU & !is.na(wageChange)) )
EEreadcount <- with(wageChangesRead,sum(EE & !is.na(wageChange)) )

wageChangesFull <- wageChangesRead

# Only count EEs and balanced EUs and UEs
wageChangesFull <- wageChangesFull %>%
	group_by(id) %>%
	arrange(id, date) %>%
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
	mutate(perwt = mean(wpfinwgt, na.rm = TRUE)) %>%
	ungroup

# Fallick/Fleischman numbers
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
SNpop <- (31.5270935961 + 0.8527093596+1.5384236453+0.8620689655+1.6463054187)/100 #includes transtions to and from N
# probability conditional on bein in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)
ratio <- UEprob/EEprob

# re-weight for total distribution
UEweight <- with(wageChangesFull, sum(perwt[UE & !is.na(wageChange)], na.rm = TRUE))
EUweight <- with(wageChangesFull, sum(perwt[EU & !is.na(wageChange)], na.rm = TRUE))
EEweight <- with(wageChangesFull, sum(perwt[EE & !is.na(wageChange)], na.rm = TRUE))

# Force weights to match ratio
multiple <- ratio*EEweight/EUweight
wageChangesFull <- wageChangesFull %>%
	mutate(balanceWeight = ifelse(EU | UE, perwt*multiple*0.5, perwt))

# Scale weights back to original total
weightMultiple <- with(wageChangesFull, sum(balanceWeight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE))
wageChangesFull <- wageChangesFull %>%
	mutate(balanceWeight = balanceWeight/weightMultiple)

# check weights
UEweight <- with(wageChangesFull, sum(balanceWeight[UE & !is.na(wageChange)], na.rm = TRUE))
EUweight <- with(wageChangesFull, sum(balanceWeight[EU & !is.na(wageChange)], na.rm = TRUE))
EEweight <- with(wageChangesFull, sum(balanceWeight[EE & !is.na(wageChange)], na.rm = TRUE))

# new weight ratio == Fallick Fleischmann UE/EE ratio
round((UEweight + EUweight)/EEweight, 10) == round(ratio, 10)
# new weights sum to original weights
with(wageChangesFull, sum(wpfinwgt, na.rm = TRUE) == sum(balanceWeight, na.rm = TRUE))
# EU and UE balance
with(wageChangesFull, sum(balanceWeight[UE & switchedOcc], na.rm = TRUE)) == 
	with(wageChangesFull, sum(balanceWeight[EU & switchedOcc], na.rm = TRUE))

# Dictionary --------------------------------------------------------------

# perwt: new weight with individual's average weight over panel
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


dataTable <- c()
# Full Sample --------------------
for(switchIndic in c(TRUE, FALSE)) {
	
	wageChanges <- subset(wageChangesFull, switchedOcc == switchIndic)
	
	# All job changes ---------------------------------------------------------
	
	# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
	attach(wageChanges)
	
	qtls <- c(0.1,0.25,0.5,0.75,0.9)
	# Central tendency
	full.Mean <- wtd.mean(wageChange, balanceWeight)
	full.Qs <- wtd.quantile(wageChange, balanceWeight, probs = qtls)
	full.Med <- full.Qs[3]
	p50 <- full.Qs[3]
	
	# Dispersion
	full.Var <- wtd.var(wageChange, balanceWeight)
	full.IQR <- full.Qs[4] - full.Qs[2]
	p90 <- full.Qs[5]
	p10 <- full.Qs[1]
	full.9010 <- p90 - p10
	
	# Skewness
	full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	full.PearsonSkew <- mean( balanceWeight*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceWeight,na.rm=T)
	
	fullRow <- c("full", switchIndic, round(c(full.Mean, full.Med, full.Var, full.IQR, 
						  full.9010, full.Kelly, full.PearsonSkew), 4))
	
	# EUE changes -------------------------------------------------------------
	
	# Central tendency
	EUE.Mean <- wtd.mean(wageChange_EUE, balanceWeight)
	EUE.Qs <- wtd.quantile(wageChange_EUE, balanceWeight, probs = qtls)
	EUE.Med <- EUE.Qs[3]
	p50 <- EUE.Qs[3]
	
	# Dispersion
	EUE.Var <- wtd.var(wageChange_EUE, balanceWeight)
	EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
	p90 <- EUE.Qs[5]
	p10 <- EUE.Qs[1]
	EUE.9010 <- p90 - p10
	
	# Skewness
	EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	EUE.PearsonSkew <- mean(balanceWeight*(wageChange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceWeight ,na.rm=T)
	
	EUERow <- c("EUE", switchIndic, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
						EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
	
	# EE distribution --------------------------------------------------------
	
	# Central tendency
	EE.Mean <- wtd.mean(wageChange[EE], balanceWeight[EE])
	EE.Qs <- wtd.quantile(wageChange[EE], balanceWeight[EE], probs = qtls)
	EE.Med <- EE.Qs[3]
	p50 <- EE.Med
	
	# Dispersion
	EE.Var <- wtd.var(wageChange[EE], balanceWeight[EE])
	EE.IQR <- EE.Qs[4]-EE.Qs[2]
	p90 <- EE.Qs[5]
	p10 <- EE.Qs[1]
	EE.9010 <- p90 - p10
	
	# Skewness
	EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	EE.PearsonSkew <- mean(balanceWeight[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceWeight[EE] ,na.rm=T)
	
	EERow <- c("EE", switchIndic, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
					      EE.9010, EE.Kelly, EE.PearsonSkew), 4))
	
	dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
	
	detach(wageChanges)
}

dataTable <- data.frame(dataTable)
labs <-  c( NA, NA,NA,NA,NA,NA,NA)
names(dataTable) <- c("Distribution", "Switched", "Mean", "Median", "Variance", "IQR", "90-10", "Kelly", "Pearson")
levels(dataTable$Switched) <- c("No Switch", "Switch")
levels(dataTable$Distribution) <- c("EE", "EUE", "Full Labor Force")
xt.dataTable <- xtable(cbind(dataTable[c(-2,-9)]), digits = 1, caption = "Month-to-Month Earnings Changes Among Job Changers", label="tab:wCh_dist")
print(xt.dataTable,file="wCh_dist.tex",table.placement="htb",
	  align = "ll|rr|rrr|r",
	  hline.after=c(-1,-1,0,nrow(xt.dataTable)/2,nrow(xt.dataTable)))

# Young/Old --------------------------------------------------------------------

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
		full.Mean <- wtd.mean(wageChange, balanceWeight)
		full.Qs <- wtd.quantile(wageChange, balanceWeight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wageChange, balanceWeight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( balanceWeight*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceWeight,na.rm=T)
		
		fullRow <- c("full", NA, switchIndic, youngIndic, round(c(full.Mean, full.Med, full.Var, full.IQR, 
			     full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes -------------------------------------------------------------
		
		# Central tendency
		EUE.Mean <- wtd.mean(wageChange_EUE, balanceWeight)
		EUE.Qs <- wtd.quantile(wageChange_EUE, balanceWeight, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wageChange_EUE, balanceWeight)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(balanceWeight*(wageChange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceWeight ,na.rm=T)
		
		EUERow <- c("EUE", NA, switchIndic, youngIndic, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
			    EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution --------------------------------------------------------
		
		# Central tendency
		EE.Mean <- wtd.mean(wageChange[EE], balanceWeight[EE])
		EE.Qs <- wtd.quantile(wageChange[EE], balanceWeight[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wageChange[EE], balanceWeight[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(balanceWeight[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceWeight[EE] ,na.rm=T)
		
		EERow <- c("EE", NA, switchIndic, youngIndic, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
			   EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wageChanges)
	}
	
}



# HS/College loop --------------------------
for(HSColIndic in c(0, 1, 2)){
	
	wageChangesHSCol <- subset(wageChangesFull, HSCol == HSColIndic)
	
	for(switchIndic in c(TRUE, FALSE)) {
		
		wageChanges <- subset(wageChangesHSCol, switchedOcc == switchIndic)
		
		# All job changes ---------------------------------------------------------
		
		# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
		attach(wageChanges)
		
		qtls <- c(0.1,0.25,0.5,0.75,0.9)
		# Central tendency
		full.Mean <- wtd.mean(wageChange, balanceWeight)
		full.Qs <- wtd.quantile(wageChange, balanceWeight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wageChange, balanceWeight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( balanceWeight*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceWeight,na.rm=T)
		
		fullRow <- c("full", HSColIndic, switchIndic, NA, round(c(full.Mean, full.Med, full.Var, full.IQR, 
			     full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes -------------------------------------------------------------
		
		# Central tendency
		EUE.Mean <- wtd.mean(wageChange_EUE, balanceWeight)
		EUE.Qs <- wtd.quantile(wageChange_EUE, balanceWeight, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wageChange_EUE, balanceWeight)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(balanceWeight*(wageChange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceWeight ,na.rm=T)
		
		EUERow <- c("EUE", HSColIndic, switchIndic, NA, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
			    EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution --------------------------------------------------------

		# Central tendency
		EE.Mean <- wtd.mean(wageChange[EE], balanceWeight[EE])
		EE.Qs <- wtd.quantile(wageChange[EE], balanceWeight[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wageChange[EE], balanceWeight[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(balanceWeight[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceWeight[EE] ,na.rm=T)
		
		EERow <- c("EE", HSColIndic, switchIndic, NA, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
			   EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wageChanges)
	}
	
}

# Results table formatting
dataTable <- data.frame(dataTable)
names(dataTable) <- c("Distribution", "HSCol", "Switched", "Young", "Mean", "Median", "Variance", "IQR", "90-10", "Kelly", "Pearson")
levels(dataTable$Switched) <- c("No Switch", "Switch")
levels(dataTable$Young) <- c("Old", "Young")
levels(dataTable$HSCol) <- c("<HS", "HS", "College")

for(varNo in 5:11) {
	dataTable[,varNo] <- as.numeric(as.character(dataTable[,varNo]))
}

YoungOld <- subset(dataTable, !is.na(Young))
YoungOld$HSCol <- NULL

HSCol <- subset(dataTable, !is.na(HSCol))
HSCol$Young <- NULL

library(reshape2)
YoungOld2 <- melt(YoungOld, id = c("Distribution", "Switched", "Young"))
YoungOld2 <- reshape(YoungOld2, timevar = c("Switched"), idvar = c("Distribution", "Young", "variable"), direction = "wide")
YoungOld2 <- reshape(YoungOld2, timevar = c("Young"), idvar = c("Distribution", "variable"), direction = "wide")
YoungOld2 <- reshape(YoungOld2, timevar = c("Distribution"), idvar = c("variable"), direction = "wide")
rownames(YoungOld2) <- NULL

HSCol2 <- melt(HSCol, id = c("Distribution", "Switched", "HSCol"))
HSCol2 <- reshape(HSCol2, timevar = "Switched", idvar = c("Distribution", "HSCol", "variable"), direction = "wide")
HSCol2 <- reshape(HSCol2, timevar = "HSCol", idvar = c("Distribution", "variable"), direction = "wide")
HSCol2 <- reshape(HSCol2, timevar = "Distribution", idvar = c("variable"), direction = "wide")
rownames(HSCol2) <- NULL

new.names = c()
for(name in names(YoungOld2)){
	new.names <- c(new.names, gsub("value.", x = name, replacement = ""))		       
}
names(YoungOld2) <- new.names

new.names = c()
for(name in names(HSCol2)){
	new.names <- c(new.names, gsub("value.", x = name, replacement = ""))		       
}
names(HSCol2) <- new.names

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching Review/Tables")
YoungOld.xt <- xtable(YoungOld2, digits = 3, caption = "Month-to-Month Earnings Changes Among Job Movers by Age")
print(YoungOld.xt,file="YoungOld.tex",hline.after=c(-1,-1,0,nrow(YoungOld2)))

HSCol.xt <- xtable(HSCol2, digits = 3, caption = "Month-to-Month Earnings Changes Among Job Movers by Education")
print(HSCol.xt,file="HSCol.tex",hline.after=c(-1,-1,0,nrow(HSCol2)))
	

# Var decomp of wage changes btwn EE, EUE, stayers ------------------------

analytic9608<-readRDS("analytic9608.RData")

# how many EU's here?
#sum(analytic9608$EU[analytic9608$lfStat==1]*analytic9608$wpfinwgt[analytic9608$lfStat==1],na.rm=T)/sum(analytic9608$wpfinwgt[analytic9608$lfStat==1],na.rm=T)
#0.006766477
analytic9608 <- analytic9608 %>%
	group_by(id) %>%
	mutate(perwt = mean(wpfinwgt, na.rm = TRUE)) %>%
	ungroup
analytic9608$perwt[is.na(analytic9608$perwt)] = 0
analytic9608$EE[is.na(analytic9608$EE)] <-F
analytic9608$UE[is.na(analytic9608$UE)] <-F
analytic9608$EU[is.na(analytic9608$EU)] <-F
# balance UEs and EUs 
analytic9608 <- analytic9608 %>%
	arrange(id,date) 
analytic9608$lagEUE <- ifelse(analytic9608$EU & is.finite(analytic9608$wageChange_EUE),T,NA)
analytic9608$lagEUE <- ifelse(lag(analytic9608$UE) & analytic9608$lfStat==2,F,analytic9608$lagEUE)
analytic9608$lagEUE <- ifelse(lead(analytic9608$id)!=analytic9608$id,F,analytic9608$lagEUE)
analytic9608$lagEUE <- na.locf(analytic9608$lagEUE, na.rm=F)
analytic9608$lagEUE <- lag(analytic9608$lagEUE)
sum(analytic9608$lagEUE & analytic9608$UE,na.rm=T)


# Fallick/Fleischman numbers
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
SNpop <- (31.5270935961 + 0.8527093596+1.5384236453+0.8620689655+1.6463054187)/100 #includes transtions to and from N
# probability conditional on bein in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)

#reweight to match flow rates from CPS
wpfinEE = with(analytic9608, sum(perwt[EE & is.finite(wageChange_all)]))
wpfinEU = with(analytic9608, sum(perwt[EU & is.finite(wageChange_all)]))
wpfinUE = with(analytic9608, sum(perwt[UE & is.finite(wageChange_all)]))
tot = with(analytic9608, sum(perwt[is.finite(wageChange_all)]))
wpfinStay = (1. - wpfinEE/tot-wpfinEU/tot-wpfinUE/tot)

wpfinEE_EUE = with(analytic9608, sum(perwt[EE & is.finite(wageChange_EUE)]))
wpfinEU_EUE = with(analytic9608, sum(perwt[EU & is.finite(wageChange_EUE)]))
wpfinUE_EUE = wpfinEU_EUE
tot_EUE = with(analytic9608, sum(perwt[is.finite(wageChange_EUE)]))
wpfinStay_EUE = (1. - wpfinEE_EUE/tot_EUE-wpfinEU_EUE/tot_EUE-wpfinUE_EUE/tot_EUE)

analytic9608$balwt <- analytic9608$perwt
analytic9608$balwt[analytic9608$EE] <- analytic9608$perwt[analytic9608$EE]*EEprob/(wpfinEE/tot)
analytic9608$balwt[analytic9608$EU] <- analytic9608$perwt[analytic9608$EU]*UEprob/(wpfinEU/tot)
analytic9608$balwt[analytic9608$UE] <- analytic9608$perwt[analytic9608$UE]*UEprob/(wpfinUE/tot)
analytic9608$balwt[!analytic9608$UE & !analytic9608$EU & !analytic9608$EE] <- 
	analytic9608$perwt[!analytic9608$UE & !analytic9608$EU & !analytic9608$EE]/wpfinStay*(1.- EEprob - 2*UEprob)


analytic9608$balEUEwt <- analytic9608$perwt
analytic9608$balEUEwt[analytic9608$EE] <- analytic9608$perwt[analytic9608$EE]*EEprob/(wpfinEE_EUE/tot_EUE)
analytic9608$balEUEwt[analytic9608$EU] <- analytic9608$perwt[analytic9608$EU]*UEprob/(wpfinEU_EUE/tot_EUE)
analytic9608$balEUEwt[analytic9608$UE] <- analytic9608$perwt[analytic9608$UE]*UEprob/(wpfinUE_EUE/tot_EUE)
analytic9608$balEUEwt[!analytic9608$UE & !analytic9608$EU & !analytic9608$EE] <- 
	analytic9608$perwt[!analytic9608$UE & !analytic9608$EU & !analytic9608$EE]/wpfinStay_EUE*(1.- EEprob - 2*UEprob)

analytic9608$EE_fac <- ifelse(analytic9608$EE,1,0)
analytic9608$EUUE_fac <- ifelse(analytic9608$EU | analytic9608$UE | analytic9608$lfStat == 2,1,0)
analytic9608$EUE_fac <- ifelse(analytic9608$EU,1,0)
analytic9608$swEE_fac <- ifelse(analytic9608$EE & analytic9608$switchedOcc,1,0)
analytic9608$swEUUE_fac <- ifelse((analytic9608$EU|analytic9608$UE|analytic9608$lfStat == 2) & analytic9608$switchedOcc,1,0)
analytic9608$swEUE_fac <- ifelse(analytic9608$EU & analytic9608$switchedOcc,1,0)
analytic9608$stay <- ifelse(analytic9608$EE | analytic9608$EU | analytic9608$UE, 0 ,1)

wageChangeMean = wtd.mean(analytic9608$wageChange_all, analytic9608$balwt, na.rm=T)
wageChangeSS   = with(analytic9608, sum( ( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_stay = with(analytic9608, sum(  stay*( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_EE = with(analytic9608, sum(  EE_fac*( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_EUUE = with(analytic9608, sum(  EUUE_fac*( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))



tot = with(analytic9608, sum(  balwt * as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_stay = with(analytic9608, sum(  stay*balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_EE = with(analytic9608, sum(  EE_fac* balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_EUUE = with(analytic9608, sum(  EUUE_fac*balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_EUE = with(analytic9608, sum(  EUE_fac*balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_swEE = with(analytic9608, sum(  swEE_fac* balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_swEUE = with(analytic9608, sum(  swEUE_fac* balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))
pct_swEUUE = with(analytic9608, sum(  swEUUE_fac* balwt* as.integer(is.finite(wageChange_all)) ,na.rm=T))


wageChangeSS_swEE = with(analytic9608, sum(  swEE_fac*( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_swEUUE =with(analytic9608,  sum(  swEUUE_fac*( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_nswEE = with(analytic9608, sum( (1.- swEE_fac)*( EE_fac)*
						  	( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))
wageChangeSS_nswEUUE =with(analytic9608,  sum( (1.- swEUUE_fac)* EUUE_fac*
								( wageChange_all- wageChangeMean)^2* balwt ,na.rm=T))

# not including the U stints
wageChangeEUEMean = wtd.mean(analytic9608$wageChange_EUE, analytic9608$balEUEwt,na.rm=T)
wageChangeEUESS   = with(analytic9608, sum( (wageChange_EUE- wageChangeEUEMean)^2*balEUEwt ,na.rm=T))
wageChangeEUESS_stay = with(analytic9608, sum(  stay*(wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_EE = with(analytic9608, sum(  EE_fac*(wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_EUE = with(analytic9608, sum( EUE_fac*(wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_swEE = with(analytic9608, sum( swEE_fac*(wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_swEUE = with(analytic9608, sum( swEUE_fac*(wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_nswEE = with(analytic9608, sum( (1.-swEE_fac)*swEE_fac*( wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))
wageChangeEUESS_nswEUE = with(analytic9608, sum( (1.-swEUE_fac)*swEUE_fac*( wageChange_EUE- wageChangeEUEMean)^2* balEUEwt ,na.rm=T))

totEUE = with(analytic9608, sum(  balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))
pctEUE_stay = with(analytic9608, sum(  stay*balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))
pctEUE_EE = with(analytic9608, sum(  EE_fac* balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))
pctEUE_EUE = with(analytic9608, sum(  EUE_fac*balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))
pctEUE_swEE = with(analytic9608, sum(  swEE_fac* balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))
pctEUE_swEUE = with(analytic9608, sum(  swEUE_fac* balEUEwt*is.finite(wageChange_EUE) ,na.rm=T))


wCh_vardec <- rbind(c(wageChangeEUESS_stay/wageChangeEUESS,wageChangeEUESS_EE/wageChangeEUESS,wageChangeEUESS_EUE /wageChangeEUESS, NA),
					c(pctEUE_stay/totEUE,pctEUE_EE/totEUE,pctEUE_EUE/totEUE, NA),
					c(wageChangeSS_stay/wageChangeSS,wageChangeSS_EE/wageChangeSS,NA, wageChangeSS_EUUE /wageChangeSS),
					c(pct_stay/tot,pct_EE/tot, NA,pct_EUUE/tot))
wCh_vardec <- wCh_vardec*100
rownames(wCh_vardec) <- c("Contrib Pct, Employed", "Population Pct, Employed", "Contrib Pct, Labor Force", "Population Pct, Labor Force")
colnames(wCh_vardec) <- c("Continuing", "EE", "EUE", "EU,UE")
wCh_vardec.xt <- xtable(wCh_vardec, digits = 1, caption = "Contribution to Variance of Month-to-Month Earnings Changes", label="tab:wCh_vardec")
print(wCh_vardec.xt,file="wCh_vardec.tex",table.placement="htb",
	  align = "l|rrrr",
	  hline.after=c(-1,-1,0,nrow(wCh_vardec)/2,nrow(wCh_vardec)))


wChsw_vardec <-rbind(c(wageChangeEUESS_swEE/wageChangeEUESS_EE,wageChangeEUESS_swEUE/wageChangeEUESS_EUE,NA),
					 c(pctEUE_swEE/pctEUE_EE,pct_swEUE/pct_EUE, NA),
					 c(wageChangeSS_swEE/wageChangeSS_EE, NA ,2*wageChangeSS_swEUUE/ wageChangeSS_EUUE),
					 c(pct_swEE/pct_EE,NA,2*pct_swEUUE/pct_EUUE))
wChsw_vardec <- wChsw_vardec*100

rownames(wChsw_vardec) <- c("Switchers Contrib Pct, Employed", "Switchers Pct, Employed", "Switchers Contrib Pct, Labor Force", "Switchers Pct, Labor Force")
colnames(wChsw_vardec) <- c("EE", "EUE", "EU,UE")
wChsw_vardec.xt <- xtable(wChsw_vardec, digits = 3, caption = "Among Job Changers, Occupation Switchers Contribution to Variance of Month-to-Month Earnings Changes", label="tab:wCh_vardec")
print(wChsw_vardec.xt,file="wChsw_vardec.tex",table.placement="htb",
	  align = "l|rrr", 
	  hline.after=c(-1,-1,0,nrow(wChsw_vardec)/2,nrow(wChsw_vardec)))
