# January 28, 2016
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

# run the Daly, Hobijn, Wiles decomposition and various percentiles

cdistTot <- wtd.quantile(wagechangesfull$wagechange_EUE, wagechangesfull$balanceweight, 
						 probs = seq(0.01,0.99,by=0.01))
wagechangesRec <- subset(wagechangesfull, recIndic==T)
wagechangesExp <- subset(wagechangesfull, recIndic==F)

cdistRec <- wtd.quantile(wagechangesRec$wagechange_EUE, wagechangesRec$balanceweight, 
						 probs = seq(0.01,0.99,by=0.01))
cdistExp <- wtd.quantile(wagechangesExp$wagechange_EUE, wagechangesExp$balanceweight, 
						 probs = seq(0.01,0.99,by=0.01))
