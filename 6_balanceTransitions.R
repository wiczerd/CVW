# January 22, 2016
# Balance transitions
# 1) 
library(data.table)
library(zoo)
library(Hmisc)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

recallRecodeShorTerm <- function(DF){
	DF[job> 0 , jobpos := job]
	DF[ EmpTmrw==T, nextjob := shift(jobpos,1,type="lead"), by=id]
	DF[ is.finite(stintid) & (lfstat>=2 | EU==T) , nextjob := Mode(nextjob), by=list(id,stintid)]
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
	#get the remaining recalled:
	DF[EU , recalled := ifelse(jobpos == nextjob,1,recalled) ]
	DF[(lfstat>=2 | EU) & !is.na(stintid), recalled := max(recalled,na.rm=T), by=list(id,stintid)]
	#DF[ , EU:= ifelse( recalled > 0.75,F,EU), by=id]
	#DF[ , UE:= ifelse( recalled > 0.75,F,UE), by=id]
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

saveRDS(DTall,"./Data/DTall_6.RData")

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




# re-weight for total distribution
UEweight.balanced <- wagechanges[UE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EUweight.balanced <- wagechanges[EU & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EEweight.balanced <- wagechanges[EE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]

# force weights to match pre-balancing ratio
multiplier <- max((EUnorecallweight/EEnorecallweight)*(EEweight.balanced/EUweight.balanced),
				  (UEnorecallweight/EEnorecallweight)*(EEweight.balanced/UEweight.balanced))
wagechanges[, balanceweight := ifelse(EU | UE, perwt*multiplier, perwt)]

# change weight among unemployed to match <15 weeks, 15-26 weeks, 27+ from CPS
CPSunempdur <- readRDS("./InputData/CPSunempDurDist.RData")
wagechanges <- merge(wagechanges, CPSunempdur, by = "date", all.x = TRUE)
wagechanges$year <-as.integer(format(wagechanges$date, "%Y") )
wagechanges[ UE ==T, year := shift(year,1,type="lag")] # be sure UE is in the same year as EU

for( yi in seq(min(wagechanges$year, na.rm=T),max(wagechanges$year, na.rm=T)  )){
	#this is the failure rate
	wagechanges[year == yi,SIPPmax_LT15:=wtd.mean( (maxunempdur<15*12/52) , perwt)]
	wagechanges[year == yi,SIPPmax_15_26:=wtd.mean( (maxunempdur>=15*12/52 & maxunempdur<=26*12/52) , perwt)]
	wagechanges[year == yi,SIPPmax_GT26:=wtd.mean( (maxunempdur>26*12/52) , perwt)]
	#this is the cumulative survival rate
	wagechanges[year == yi,CPSsurv_LT15:=(mean(FRM5_14)+ mean(LT5))/100]
	wagechanges[year == yi,CPSsurv_15_26:=(mean(FRM15_26))/100]
	wagechanges[year == yi,CPSsurv_GT26:=(mean(GT26))/100]
	#1-durwt_LT15*SIPPmax_LT15 = CPSsurv_GT26+ CPSsurv_15_26
	wagechanges[year == yi & (maxunempdur<15*12/52), durwt :=  (1.-CPSsurv_GT26+ CPSsurv_15_26)/SIPPmax_LT15]
	#durwt_FRM1526*SIPPmax_15_26 = CPSsurv_15_26
	wagechanges[year == yi & (maxunempdur>=15*12/52 & maxunempdur<=26*12/52), durwt :=  CPSsurv_15_26/SIPPmax_15_26]
	#durwt_GT26*SIPPmax_GT26 = CPSsurv_GT26
	wagechanges[year == yi & (maxunempdur>26*12/52), durwt :=  CPSsurv_GT26/SIPPmax_GT26]
}
# quantile(wagechanges$durwt, na.rm=T,probs = seq(.1,.9,.1))
# 10%       20%       30%       40%       50%       60%       70%       80%       90% 
# 0.4973496 0.5356923 1.2305479 1.3949508 1.4803896 1.4812019 1.5501993 1.5722958 1.6393264 
wagechanges[is.na(durwt), durwt := 1.]
#do not increase the incidence of unemployment
balwtEU <- sum(wagechanges$balanceweight[wagechanges$EU], na.rm=T)
balwtUE <- sum(wagechanges$balanceweight[wagechanges$UE], na.rm=T)

wagechanges$balanceweight <- wagechanges$balanceweight*wagechanges$durwt
wagechanges$balanceweight[wagechanges$EU] <- wagechanges$balanceweight[wagechanges$EU]/
	sum(wagechanges$balanceweight[wagechanges$EU], na.rm=T)*balwtEU
wagechanges$balanceweight[wagechanges$UE] <- wagechanges$balanceweight[wagechanges$UE]/
	sum(wagechanges$balanceweight[wagechanges$UE], na.rm=T)*balwtUE

wagechanges[,c("year","durwt"):=NULL]


# scale weights back to original total
wtscale <- wagechanges[, sum(balanceweight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE)]
wagechanges[, balanceweight := balanceweight/wtscale]

#- Diagnostics

# check weights
UEweight <- wagechanges[UE & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]
EUweight <- wagechanges[EU & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]
EEweight <- wagechanges[EE & !is.na(wagechange), sum(balanceweight, na.rm = TRUE)]


# new weight ratio diff from Fallick Fleischmann UE/EE ratio
# Fallick/Fleischman CPS numbers, just for comparison
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
# includes transtions to and from N
SNpop <- (31.5270935961 + 0.8527093596 + 1.5384236453 + 0.8620689655 + 1.6463054187)/100
# probability conditional on being in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)
round(UEweight/EEweight, 10) - round(UEprob/EEprob, 10)

# new weights sum to original weights
wagechanges[, sum(wpfinwgt, na.rm = TRUE)] == wagechanges[, sum(balanceweight, na.rm = TRUE)]

# EU and UE balance
wagechanges[UE & switchedOcc, sum(balanceweight, na.rm = TRUE)] == 
	wagechanges[EU & switchedOcc, sum(balanceweight, na.rm = TRUE)]

# store balanced data
saveRDS(wagechanges, "./Data/balancedwagechanges.RData")
#rm(list=ls())