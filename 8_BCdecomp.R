# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)


wagechangesfull <- readRDS("./Data/balancedwagechanges.RData")
wcRec <- subset(wagechangesfull, recIndic==T)
wcExp <- subset(wagechangesfull, recIndic==F)

distptsRec <-wtd.quantile(wcRec$wagechange, wcRec$balanceweight, probs = seq(0.01,0.99,by=0.005))
distptsExp <-wtd.quantile(wcExp$wagechange, wcExp$balanceweight, probs = seq(0.01,0.99,by=0.005))

invdistRec <- approxfun(seq(0.01,0.99,by=0.005),distptsRec)
invdistExp <- approxfun(seq(0.01,0.99,by=0.005),distptsExp)
distRec <- approxfun(distptsRec,seq(0.01,0.99,by=0.005))
distExp <- approxfun(distptsExp,seq(0.01,0.99,by=0.005))


# 6 subgroups, Sw X (EE UE EU)
distptsRecSwEEUEEU[,1] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==T & wcRec$UE==F & wcRec$EU ==F), 
	wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsRecSwEEUEEU[,2] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==F & wcRec$UE==T & wcRec$EU ==F), 
	wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsRecSwEEUEEU[,3] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==F & wcRec$UE==F & wcRec$EU ==T), 
	wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsRecSwEEUEEU[,4] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==T & wcRec$UE==F & wcRec$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsRecSwEEUEEU[,5] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==F & wcRec$UE==T & wcRec$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsRecSwEEUEEU[,6] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==F & wcRec$UE==F & wcRec$EU ==T), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))

distptsExpSwEEUEEU[,1] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==T & wcExp$UE==F & wcExp$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsExpSwEEUEEU[,2] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==F & wcExp$UE==T & wcExp$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsExpSwEEUEEU[,3] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==F & wcExp$UE==F & wcExp$EU ==T), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsExpSwEEUEEU[,4] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==T & wcExp$UE==F & wcExp$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsExpSwEEUEEU[,5] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==F & wcExp$UE==T & wcExp$EU ==F), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))
distptsExpSwEEUEEU[,6] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==F & wcExp$UE==F & wcExp$EU ==T), 
							   wtd.quantile(wagechange, balanceweight, probs = seq(0.01,0.99,by=0.005)))



qtls<- seq(0.1,0.9,by=0.1)

for (qi in qtls){
	
}