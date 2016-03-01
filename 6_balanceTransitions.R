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

recallRecodeJobID <- function(DF){
	#this function just uses job ID to identify recalls, and may miss especially long-term when job id can reset.  
	DF[job> 0 , jobpos := job]
	DF[ EmpTmrw==T, nextjob := shift(jobpos,1,type="lead"), by=id]
	DF[ is.finite(stintid) & (lfstat>=2 | EU==T) , nextjob := Mode(nextjob), by=list(id,stintid)]
	# will convert all recall stints as lfstat == NA_integer_
	DF[EU==T , recalled := (jobpos == nextjob) ]
	DF[is.finite(stintid) &(EU==T | lfstat>=2 ), recalled:=Mode(recalled) , by=list(id,stintid)]
	DF[EU==T & recalled==T, EU:=F]
	DF[UE==T & recalled==T, UE:=F]
	DF[lfstat>=2 & recalled==T, lfstat:=0]
}

recallRecodeShorTerm <- function(DF){
	#this function closely follows Fujita & Moscarini for identification of short-term recalls
	DF[job> 0 , jobpos := job]
	DF[ EmpTmrw==T, nextjob := shift(jobpos,1,type="lead"), by=id]
	DF[ is.finite(stintid) & (lfstat>=2 | EU==T) , nextjob := Mode(nextjob), by=list(id,stintid)]
	# will convert all recall stints as lfstat == NA_integer_
	DF[ , ENEnoseam := lfstat >=2 & ( shift(lfstat,1,type="lead")==1 & shift(lfstat,1,type="lag")==1
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
	DF[EU==T , recalled := ifelse(jobpos == nextjob,1,recalled) ]
	DF[(lfstat>=2 | EU==T) & !is.na(stintid) & stintid>0, recalled := max(recalled,na.rm=T), by=list(id,stintid)]
	DF[ , EUrecalled:= ( recalled > 0.75 & EU==T), by=id]
	DF[ , UErecalled:= ( recalled > 0.75 & UE==T), by=id]
	
	DF[ recalled > 0.75, EU:= F, by=id]
	DF[ recalled > 0.75, UE:= F, by=id]
	DF[ recalled > 0.75 & lfstat>=2, lfstat := 0L, by=id]
	DF[ , c("jobpos","ENEwseam","ENEnoseam") := NULL]
}

DTall <- readRDS("./Data/DTall_5.RData")

# sum weights for UE, EU, and EE
UEreadweight <- DTall[UE==T & !is.na(wagechange) & is.finite(stintid), sum(wpfinwgt, na.rm = TRUE)]
EUreadweight <- DTall[EU==T & !is.na(wagechange) & is.finite(stintid), sum(wpfinwgt, na.rm = TRUE)]
EEreadweight <- DTall[EE==T & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]

DTall[lfstat==2, recalled:= F ]
recallRecodeJobID(DTall)
#DTall[lfstat==2, recalled:= 0. ]
#recallRecodeShorTerm(DTall)

UEnorecallweight <- DTall[UE==T & !is.na(wagechange) & is.finite(stintid), sum(wpfinwgt, na.rm = TRUE)]
EUnorecallweight <- DTall[EU==T & !is.na(wagechange) & is.finite(stintid), sum(wpfinwgt, na.rm = TRUE)]
EEnorecallweight <- DTall[EE==T & !is.na(wagechange), sum(wpfinwgt, na.rm = TRUE)]

#recall rate:
UEnorecallweight/UEreadweight

#saveRDS(DTall,"./Data/DTall_6.RData")
#some rates: diagnostic
sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE & is.finite(DTall$stintid)==T], na.rm=T)/sum(DTall$wpfinwgt[is.finite(DTall$stintid)==T], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[!is.na(DTall$unempdur)]*DTall$unempdur[!is.na(DTall$unempdur)],na.rm=T)/sum(DTall$wpfinwgt[!is.na(DTall$unempdur)],na.rm=T)
#end diagnostic

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
wagechanges[UE==T & shift(switchedOcc, 1, type = "lag")==T, switchedOcc := T, by = id]
wagechanges[, maxunempdur:= ifelse(UE==T & is.na(maxunempdur), shift(maxunempdur), maxunempdur) , by = id]
wagechanges[, maxunempdur:= ifelse(UE==T & maxunempdur<shift(maxunempdur), shift(maxunempdur), maxunempdur) , by = id]
wagechanges[, maxunempdur:= ifelse(EU==T & maxunempdur<shift(maxunempdur,1,type="lead"), shift(maxunempdur,1,type="lead"), maxunempdur) , by = id]
wagechanges[, stintid := ifelse( EU==T & is.na(stintid), shift(stintid,1,type="lead"), stintid ), by=id]
wagechanges[, stintid := ifelse( UE==T & is.na(stintid), shift(stintid,1,type="lag"), stintid ), by=id]

# set HSCol and Young to max over panel to ensure balance
wagechanges[, HSCol := max(HSCol), by = id]
wagechanges[, Young := as.logical(max(Young)), by = id]

# create new person weights
wagechanges[, perwt := mean(wpfinwgt, na.rm = TRUE), by = id]
DTall[, perwt := mean(wpfinwgt, na.rm = TRUE), by = id]


# re-weight for total distribution
UEweight.balanced <- wagechanges[UE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EUweight.balanced <- wagechanges[EU & !is.na(wagechange), sum(perwt, na.rm = TRUE)]
EEweight.balanced <- wagechanges[EE & !is.na(wagechange), sum(perwt, na.rm = TRUE)]

# re-inflate weights for workers who enter and leave as unemployed & divide by 2 for the whole transition
# this should take the UE, because many EU will leave sample by exit LF
multiplier <- (UEnorecallweight/EEnorecallweight)*(EEweight.balanced/UEweight.balanced)

wagechanges[, balanceweight := ifelse(EU | UE, perwt*multiplier, perwt)]

# change weight among unemployed to match <15 weeks, 15-26 weeks, 27+ from CPS
CPSunempdur <- readRDS("./InputData/CPSunempDurDist.RData")
wagechanges <- merge(wagechanges, CPSunempdur, by = "date", all.x = TRUE)
wagechanges$year <-as.integer(format(wagechanges$date, "%Y") )
wagechanges[ UE ==T, year := shift(year,1,type="lag")] # be sure UE is in the same year as EU

#for( yi in seq(min(wagechanges$year, na.rm=T),max(wagechanges$year, na.rm=T)  )){
	#this is the failure rate
	wagechanges[,SIPPmax_LT15:=wtd.mean( (maxunempdur<15*12/52) , perwt)]
	wagechanges[,SIPPmax_15_26:=wtd.mean( (maxunempdur>=15*12/52 & maxunempdur<=26*12/52) , perwt)]
	wagechanges[,SIPPmax_GT26:=wtd.mean( (maxunempdur>26*12/52) , perwt)]
	#this is the cumulative survival rate
	wagechanges[,CPSsurv_LT15:=(mean(FRM5_14)+ mean(LT5))/100]
	wagechanges[,CPSsurv_15_26:=(mean(FRM15_26))/100]
	wagechanges[,CPSsurv_GT26:=(mean(GT26))/100]
	#1-durwt_LT15*SIPPmax_LT15 = CPSsurv_GT26+ CPSsurv_15_26
	wagechanges[(maxunempdur<15*12/52), durwt :=  (1.-CPSsurv_GT26+ CPSsurv_15_26)/SIPPmax_LT15]
	#durwt_FRM1526*SIPPmax_15_26 = CPSsurv_15_26
	wagechanges[(maxunempdur>=15*12/52 & maxunempdur<=26*12/52), durwt :=  CPSsurv_15_26/SIPPmax_15_26]
	#durwt_GT26*SIPPmax_GT26 = CPSsurv_GT26
	wagechanges[(maxunempdur>26*12/52), durwt :=  CPSsurv_GT26/SIPPmax_GT26]
#}
wagechanges[ (maxunempdur<15*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
wagechanges[ (maxunempdur>=15*12/52 & maxunempdur<=26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
wagechanges[ (maxunempdur>26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]


#wagechanges[is.na(durwt), durwt := 1.]
wagechanges[, durwt := 1.]
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

# create EUE balance weights that give double weight to the EU and zero to the UE
wagechanges[, balanceweightEUE := balanceweight]
wagechanges[EU==T, balanceweightEUE := balanceweight*2]
wagechanges[UE==T, balanceweightEUE := 0.]


wagechangesBalanced<-subset(wagechanges, select=c("id","date","balancedEU","balancedUE","maxunempdur","balanceweight"))

setkey(DTall,id,date)
DTall[is.finite(stintid), balancedEU := max(UE,na.rm=T)==T & EU==T, by = list(id,stintid)]
DTall[is.finite(stintid), balancedUE := max(EU,na.rm=T)==T & UE==T, by = list(id,stintid)]
DTall<- merge(DTall,wagechangesBalanced,by=c("id","date"),all.x=T)

#make the wage changes NA if not-balanced
DTall[UE==T & balancedUE.y==F, c("wagechange_EUE","wagechange_all","wagechange") := NA_real_]
DTall[EU==T & balancedEU.y==F, c("wagechange_EUE","wagechange_all","wagechange") := NA_real_]

DTall[UE==T, UE := balancedUE.y==T]
DTall[EU==T, EU := balancedEU.y==T]
#there are a few from x that don't have stintid, fill these with y.

DTall[ is.finite(stintid), maxunempdur := max(c(maxunempdur.y, maxunempdur.x), na.rm=T), by=list(id,stintid)] #trusting the DTall worked right.
#DTall[ is.finite(stintid) , completestintUE:= as.integer(balancedUE.y==T) ] <- this part is unecessary
#DTall[ is.finite(stintid) , completestintUE:= max(completestintUE) , by = list(id,stintid)]
#DTall[ is.finite(stintid) , completestintEU:= as.integer(balancedEU.y==T) ]
#DTall[ is.finite(stintid) , completestintEU:= max(balancedEU.y), by=list(id,stintid) ]
DTall[ !is.finite(balanceweight), balanceweight:= 0.]
DTall[, balanceweight := max(balanceweight, na.rm=T), by=list(id,stintid)]
DTall[, c("maxunempdur.x","maxunempdur.y") := NULL ]

saveRDS(DTall,"./Data/DTall_6.RData")


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