# January 22, 2016
# Balance transitions
# 1) 
library(data.table)
library(zoo)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)


recallRecodeShorTerm <- function(DF){
	DF[job> 0 , jobpos := job]
	# will convert all recall stints as lfstat == NA_integer_
	DF[ , ENEnoseam := lfstat ==2 & ( shift(lfstat,1,type="lead")==1 & shift(lfstat,1,type="lag")==1
									  |  (sum(lfstat ==1, na.rm=T)==2 & (shift(lfstat,1,type="lead")==2 |shift(lfstat)==2)) )
	   , by=list(id,wave)]
	DF[ is.na(ENEnoseam)==T, ENEnoseam:=F ]
	DF[ , recalled := as.numeric( ENEnoseam==T & ( (shift(jobpos)==shift(jobpos,1,type="lead") & maxunempdur==1) |
												   	( maxunempdur==2 & 
												   	  	(shift(jobpos,2,type="lag")==shift(jobpos,1,type="lead")|
												   	  	shift(jobpos,1,type="lag")==shift(jobpos,2,type="lead")) ) )),  
	   by = list( id,wave ) ]
	#impute recalls if there is a seam
	DF[ , ENEwseam := ( maxunempdur<=2 & lfstat ==2 ) &
	   	(wave != shift(wave) | wave != shift(wave,1,type="lead")) ,by=id]
	DF[is.na(ENEwseam), ENEwseam:=F]
	predrecall <- glm( recalled ~ wagechange_EUE + I(wagechange_EUE^2) + recIndic + factor(esr) + union + 
					   	factor(HSCol) + factor(Young)  + I(wagechange_EUE>-.5) + I(wagechange_EUE<.5 & wagechange_EUE>-.5)
					   + switchedOcc + switchedAddress, 
					   na.action = na.exclude, data= subset(DF,ENEnoseam),family=binomial(link="probit"))
	DF[ (ENEwseam | ENEnoseam), 
		fitted_recall := predict(predrecall, type= "response", newdata=subset(DF,(ENEwseam==T | ENEnoseam==T) ))]
	DF[ ENEwseam==T & !is.na(fitted_recall), recalled:= fitted_recall ]
	DF[ , recalled:= ifelse( ENEwseam & (shift(jobpos,1,type="lead") == shift(jobpos,1,type="lag"))
							 , 1,recalled ),  by= id]
	DF[ , recalled:= ifelse( ENEwseam & maxunempdur==2 & 
							 	((shift(jobpos,1,type="lead") == shift(jobpos,2,type="lag")) | (shift(jobpos,2,type="lead") == shift(jobpos,1,type="lag")))
							 , 1,recalled ),  by= id]
	DF[is.na(recalled), recalled := 0.]
	DF[(lfstat>=2 | EU) & !is.na(stintid), recalled := max(recalled,na.rm=T), by=list(id,stintid)]
	DF[ , EU:= ifelse( recalled > 0.75,F,EU), by=id]
	DF[ , UE:= ifelse( recalled > 0.75,F,UE), by=id]
	DF[ , c("jobpos","ENEwseam","ENEnoseam") := NULL]
}

DTall <- readRDS("./Data/DTall_5.RData")

# sum weights for UE, EU, and EE
UEreadweight <- DTall[UE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EUreadweight <- DTall[EU & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EEreadweight <- DTall[EE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]

DTall[lfstat==2, recalled:= 0 ]
recallRecodeShorTerm(DTall)

UEnorecallweight <- DTall[UE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EUnorecallweight <- DTall[EU & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]
EEnorecallweight <- DTall[EE & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]

# create data set with defined wage changes only
wagechanges <- DTall[!is.infinite(wagechange) & !is.na(wagechange),]
setkey(wagechanges, id, date)


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

# Fallick/Fleischman CPS numbers, just for comparison
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
# includes transtions to and from N
SNpop <- (31.5270935961 + 0.8527093596 + 1.5384236453 + 0.8620689655 + 1.6463054187)/100
# probability conditional on being in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)


# re-weight for total distribution
UEweight.balanced <- wagechanges[UE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EUweight.balanced <- wagechanges[EU & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EEweight.balanced <- wagechanges[EE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]

# force weights to match pre-balancing ratio
multiplier <- (EUnorecallweight/EEnorecallweight)*EUweight.balanced/EEweight.balanced
wagechanges[, balanceweight := ifelse(EU | UE, perwt*multiplier*0.5, perwt)]

###############come back here

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