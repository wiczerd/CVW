# January 22, 2016
# Balance transitions
# 1) 
library(data.table)
library(zoo)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)


DTall <- readRDS("./Data/DTall_5.RData")

# create data set with defined wage changes only
wagechanges <- DTall[!is.infinite(wagechange) & !is.na(wagechange),]
setkey(wagechanges, id, date)

# sum weights for UE, EU, and EE
UEreadweight <- wagechanges[UE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EUreadweight <- wagechanges[EU & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EEreadweight <- wagechanges[EE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]

# keep only count EEs and balanced EUs and UEs
wagechanges[, balancedEU := EU & shift(UE, 1, type = "lead"), by = id]
wagechanges[, balancedUE := UE & shift(EU, 1, type = "lag"), by = id]
wagechanges <- wagechanges[EE | balancedEU | balancedUE,]

# change switchedocc to TRUE for UE if switchedocc is TRUE for corresponding EU
# Q: Why?
# A: To correct for how we defined the timing of the occupation switch. See 3_createVars.R.
#    Now both the EU and UE will be considered an occupation switch.
wagechanges[UE & shift(switchedOcc, 1, type = "lag"), switchedOcc := TRUE, by = id]

# set HSCol and Young to max over panel to ensure balance
wagechanges[, HSCol := max(HSCol), by = id]
wagechanges[, Young := as.logical(max(Young)), by = id]

# create new person weights
wagechanges[, perwt := mean(wpfinwgt, na.rm = TRUE), by = id]

# Fallick/Fleischman numbers
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
# includes transtions to and from N
SNpop <- (31.5270935961 + 0.8527093596 + 1.5384236453 + 0.8620689655 + 1.6463054187)/100
# probability conditional on being in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)
ratio <- UEprob/EEprob

# re-weight for total distribution
UEweight.SIPP <- wagechanges[UE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EUweight.SIPP <- wagechanges[EU & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EEweight.SIPP <- wagechanges[EE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]

# force weights to match Fallick/Fleischman ratio
multiplier <- ratio*EEweight.SIPP/EUweight.SIPP
wagechanges[, balanceweight := ifelse(EU | UE, perwt*multiplier*0.5, perwt)]

# scale weights back to original total
divisor <- wagechanges[, sum(balanceweight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE)]
wagechanges[, balanceweight := balanceweight/divisor]

# check weights
UEweight <- wagechanges[UE & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]
EUweight <- wagechanges[EU & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]
EEweight <- wagechanges[EE & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]

# new weight ratio == Fallick Fleischmann UE/EE ratio
round((UEweight + EUweight)/EEweight, 10) == round(ratio, 10)

# new weights sum to original weights
wagechanges[, sum(wpfinwgt, na.rm = TRUE)] == wagechanges[, sum(balanceweight, na.rm = TRUE)]

# EU and UE balance
wagechanges[UE & switchedOcc, sum(balanceweight, na.rm = TRUE)] == 
	wagechanges[EU & switchedOcc, sum(balanceweight, na.rm = TRUE)]

# store balanced data
saveRDS(wagechanges, "./Data/balancedwagechanges.RData")
rm(list=ls())