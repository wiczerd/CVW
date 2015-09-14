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
qtls <- c(0.1,0.25,0.5,0.75,0.9)
# Central tendency
full.Mean <- wtd.mean(wageChange, wtchng)
full.Qs <- wtd.quantile(wageChange, wtchng, probs = qtls)
full.Med <- full.Qs[3]
p50 <- full.Qs[3]

# Dispersion
full.Var <- wtd.var(wageChange, wtchng)
full.IQR <- full.Qs[4] - full.Qs[2]
p90 <- full.Qs[5]
p10 <- full.Qs[1]
full.9010 <- p90 - p10

# Skewness
full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
full.PearsonSkew <- mean( wtchng*(wageChange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(wtchng,na.rm=T)

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

# EE distribution --------------------------------------------------------

# Central tendency
EE.Mean <- wtd.mean(wageChange[EE], wpfinwgt[EE])
EE.Qs <- wtd.quantile(wageChange[EE], wpfinwgt[EE], probs = qtls)
EE.Med <- EE.Qs[3]
p50 <- EE.Med

# Dispersion
EE.Var <- wtd.var(wageChange[EE], wpfinwgt[EE])
EE.IQR <- EE.Qs[4]-EE.Q[2]
p90 <- EE.Qs[5]
p10 <- EE.Qs[1]
EE.9010 <- p90 - p10

# Skewness
EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
EE.PearsonSkew <- mean(wpfinwgt[EE]*(wageChange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(wpfinwgt[EE] ,na.rm=T)


# Put all of the central tendency into 1 table
cent <- rbind(c(full.Mean,EE.Mean,EUE.Mean),
				   c(full.Med,EE.Med,EUE.Med))
rownames(cent) <- c("Mean","Median")
colnames(cent) <- c("All","E->E","E->U->E")
cent.xt <- xtable(cent,label="tab:cent",digits=3,caption="Month-to-Month Earnings Changes Among Job Movers, Central Tendency")
print(cent.xt,file="cent.tex",hline.after=c(-1,-1,0,nrow(cent)))

disp <- rbind(c(full.Var,EE.Var,EUE.Var),
			  c(full.IQR,EE.IQR,EUE.IQR),
			  c(full.9010,EE.9010,EE.9010))
rownames(disp) <- c("Variance","IQR","90-10 range")
colnames(disp) <- c("All","E->E","E->U->E")
disp.xt <- xtable(disp,label="tab:disp",digits=3,caption="Average Month-to-Month Earnings Changes Among Job Movers, Dispersion")
print(disp.xt,file="disp.tex",hline.after=c(-1,-1,0,nrow(disp)))

skew <- rbind(c(full.PearsonSkew,EE.PearsonSkew,EUE.PearsonSkew),
			  c(full.Kelly,EE.Kelly,EUE.Kelly))
rownames(skew) <- c("Pearson","Kelly")
colnames(skew) <- c("All","E->E","E->U->E")
skew.xt <- xtable(skew,label="tab:skew",digits=3,caption="Average Month-to-Month Earnings Changes Among Job Movers, Skewness")
print(skew.xt,file="skew.tex",hline.after=c(-1,-1,0,nrow(skew)))


ggplot(wageChanges, aes(wageChange)) +
	geom_density()


