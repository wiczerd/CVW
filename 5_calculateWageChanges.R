# January 21, 2016
# Calculate wage changes
# 1) create wagechange, wagechange_stayer, wagechange_EUE, and wagechange_all variables
# 2) create occwagechange variable
# 3) save intermediate result, DTall_5.RData
library(data.table)
library(zoo)
library(stats)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

DTall <- readRDS(paste0(datadir,"/DTall_4.RData"))
#drop some un-need variables
DTall[ , nwhite := black==T|hispanic==T]
DTall <- DTall[ , c("coc","syear","smonth","eyear","emonth","esr","epppnum","ajbocc","logearnm","black","hispanic"):=NULL]
DTall <- DTall[ lfstat_wave<3, ]

setkey(DTall, id, date)
#clear out the super-low earnings
DTall[ lfstat==1 & usewage<log(20+(1+20^2)^.5), usewage:=NA]

DTall[ , next.usewage := shift(usewage,type="lead"), by=id]
DTall[ , last.usewage := shift(usewage,type="lag"), by=id]
DTall[ , c("next.earnm","last.earnm") := NULL]

# fill wages upwards to fill in unemployment spells
DTall[EE==T|UE == T, nextwage := next.usewage]
DTall[lfstat>=2 | EU==T, nextwage := Mode(nextwage), by = list(id,ustintid)] #replace within  unemp stints

DTall[, leadwage := shift(occwage, 1, type = "lead"), by = id]
DTall[EE==T|UE == T, nextoccwage := leadwage, by = id]
DTall[lfstat>=2 | EU==T, nextoccwage := Mode(nextoccwage), by = list(id,ustintid)] #replace if it's UE

DTall[, leadwage:=NULL]

DTall[last.lfstat==1, tuw := last.usewage]


# create wagechange variable---------------------------------------------
DTall[EE == T, wagechange := nextwage - tuw] #worksbetter with usewage
DTall[EU == T, wagechange := log(1.0) - tuw]
DTall[UE == T, wagechange := nextwage - log(1.0)]

# create wagechange_stayer variable------------------------------
DTall[ !(EE==T|EU==T|UE==T), wagechange_stayer := next.usewage - last.usewage]
DTall[job == 0 & shift(job,type="lead") == 0,
      wagechange_stayer := 0.0, by = id]

#do not allow stayers to change job listings without having a marked transition
DTall[ !(EE==T|EU==T|UE==T) & jobchng_wave ==T, wagechange_stayer:=NA]
#do not allow stayers to gain/lose more than 200% change in earnings w/o change in status with changes that revert
DTall[ !(EE==T|EU==T|UE==T) & (wagechange_stayer<(-2.) | wagechange_stayer>(2.))&
	   (((usewage - last.usewage)<(-2.)& (next.usewage - usewage)>(2.)) | (usewage - last.usewage)>(2.)& (next.usewage - usewage)<(-2.))
	   , wagechange_stayer_bad :=  T ]
DTall[ is.na(wagechange_stayer_bad), wagechange_stayer_bad :=  F]

DTall[ wagechange_stayer_bad==T, wagechange_stayer := NA]


# create wagechange_EUE variable -----------------------------------
DTall[EU==T | EE==T, wagechange_EUE := nextwage - tuw]
# fill across whole ustintid
DTall[lfstat>=2 | EU==T, wagechange_EUE := Mode(wagechange_EUE), by=list(id,ustintid)]
DTall[!(EU==T | EE==T | UE==T), wagechange_EUE := wagechange_stayer]

# create wagechange_month variable
# if ifelse() condition is NA, end result is NA.
DTall[!(EE | UE | EU), wagechange_month := wagechange_stayer]
DTall[(EE | UE | EU), wagechange_month := wagechange]
DTall[ , c("wagechange_stayer","wagechange","wagechange_stayer_bad"):=NULL]

# create occwagechange variable
#DTall[, occwagechange := as.numeric(ifelse(switchedOcc & !shift(switchedOcc, 1, type = "lag"), 
#				nextoccwage - shift(occwage, 1, type = "lag"), 
#				0.)), by = id]
#DTall[, occwagechange := as.numeric(ifelse(switchedOcc & shift(switchedOcc, 1, type = "lag"),
#				nextoccwage - occwage, 
#				0.)), by = id]

#***********************************************************************
#***********************************************************************
# compute seam-to-seam wage change variable ----------------------------
#***********************************************************************
#***********************************************************************

# residual wages
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]
DTall[ , wavewage := sum(levwage,na.rm=T), by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , wavewage := log(wavewage + (1+wavewage^2)^.5) ]
DTall[ , nawavewage:= all(is.na(usewage) ),by= list(id,wave)]
DTall[ nawavewage==T, wavewage:=NA_real_]
# raw wages
DTall[ , waverawwg := sum(earnm,na.rm=T), by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , waverawwg := log(waverawwg + (1+waverawwg^2)^.5) ]
DTall[ , nawavewage:= all(is.na(waverawwg) ),by= list(id,wave)]
DTall[ nawavewage==T, waverawwg:=NA_real_]

#drop the lowest wages, implies working less than $80 /month:
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
#extrmpctile<-DTall[seam==T & lfstat_wave==1, wtd.quantile(levwage/nmo_lf1,na.rm = T, probs = c(.01,.02,.98,.99),weights = wpfinwgt)]
DTall[ lfstat_wave==1 & wavewage<log(80+(1+80^2)^.5), wavewage:=NA]
DTall[ lfstat_wave==1, earn_imp_wave := sum(earn_imp==1,na.rm=T), by=list(id,wave)]
DTall[ , c("levwage","nawavewage","nmo_lf1"):=NULL]

DTseam <- DTall[ seam==T,]
#need to add change across waves (use wavewage)
DTseam[ , next.wavewage := shift(wavewage,1,type="lead"), by=id]
DTseam[ , last.wavewage := shift(wavewage,1,type="lag"), by=id]

DTseam[ , nw:= next.wavewage]
DTseam[ , tw:= wavewage]

DTseam[ , wagechange_wave := nw - tw]
#for EE that spans waves:
DTseam[ EE_wave==T & EEmon<4 , wagechange_wave:= NA]
DTseam[ EE_wave==T & shift(midEE)==T , wagechange_wave:= (next.wavewage - last.wavewage)]
DTseam[ , next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTseam[EE_wave==T & midEE ==T , wagechange_wave := next.wagechange_wave]

DTseam[ , next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTseam[ , last.wagechange_wave := shift(wagechange_wave, type="lag" ),by=id]

# create wave-level EUE wage change-----------------------------------------
# find wage in period before an EU and one period after UE
DTseam[ , EU_wave_first := EU_wave==T & !(shift(EU_wave)==T), by=id]
DTseam[ , UE_wave_last  := UE_wave==T & !(shift(UE_wave,type="lead")==T), by=id]

DTseam[last.lfstat_wave==1 & EU_wave_first == T, wageAtEU := last.wavewage]
DTseam[, wageAtEU := na.locf(wageAtEU, na.rm = F),by=list(ustintid_wave, id)]
DTseam[next.lfstat_wave==1 & UE_wave_last == T, wageAfterUE :=  next.wavewage]
DTseam[UE_wave_last == T, wagechangeEUE_wave := wageAfterUE - wageAtEU]
DTseam[, wagechangeEUE_wave:= Mode(wagechangeEUE_wave), by=list(ustintid_wave, id)]
DTseam[ EE_wave==T, wagechangeEUE_wave := wagechange_wave]
DTseam[ !(EU_wave|UE_wave|EE_wave), wagechangeEUE_wave := wagechange_wave]

DTseam[ , c("wageAtEU","wageAfterUE","EU_wave_first", "UE_wave_last"):=NULL]


#************************************************************************************
#adding the raw-earnings changes:-----------------------------------------------------
#************************************************************************************

#drop the lowest wages, implies working less than $80 /month:
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
#extrmpctile<-DTall[seam==T & lfstat_wave==1, wtd.quantile(levwage/nmo_lf1,na.rm = T, probs = c(.01,.02,.98,.99),weights = wpfinwgt)]
DTall[ lfstat_wave==1 & waverawwg<log(80+(1+80^2)^.5), waverawwg:=NA]
DTall[ lfstat_wave==1, earn_imp_wave := sum(earn_imp==1,na.rm=T), by=list(id,wave)]
DTall[ , c("nawavewage","nmo_lf1"):=NULL]

DTseam[ , next.waverawwg := shift(waverawwg,1,type="lead"), by=id]
DTseam[ , last.waverawwg := shift(waverawwg,1,type="lag"), by=id]

DTseam[ , nw:= next.waverawwg]
DTseam[ , tw:= waverawwg]

DTseam[ , rawwgchange_wave := nw - tw]
#for EE that spans waves:
DTseam[ EE_wave==T & EEmon<4 , rawwgchange_wave:= NA]
DTseam[ EE_wave==T & shift(midEE)==T , rawwgchange_wave:= (next.waverawwg - last.waverawwg)]
DTseam[ , next.rawwgchange_wave := shift(rawwgchange_wave, type="lead"),by=id]
DTseam[EE_wave==T & midEE ==T , rawwgchange_wave := next.rawwgchange_wave]

DTseam[ , next.rawwgchange_wave := shift(rawwgchange_wave, type="lead"),by=id]
DTseam[ , last.rawwgchange_wave := shift(rawwgchange_wave, type="lag" ),by=id]

# create wave-level EUE wage change-----------------------------------------
# find wage in period before an EU and one period after UE
DTseam[ , EU_wave_first := EU_wave==T & !(shift(EU_wave)==T), by=id]
DTseam[ , UE_wave_last  := UE_wave==T & !(shift(UE_wave,type="lead")==T), by=id]

DTseam[last.lfstat_wave==1 & EU_wave_first == T, wageAtEU := last.waverawwg]
DTseam[, wageAtEU := na.locf(wageAtEU, na.rm = F),by=list(ustintid_wave, id)]
DTseam[next.lfstat_wave==1 & UE_wave_last == T, wageAfterUE :=  next.rawwgchange_wave]
DTseam[UE_wave_last == T, rawwgchangeEUE_wave := wageAfterUE - wageAtEU]
DTseam[, rawwgchangeEUE_wave:= Mode(rawwgchangeEUE_wave), by=list(ustintid_wave, id)]
DTseam[ EE_wave==T, rawwgchangeEUE_wave := rawwgchange_wave]
DTseam[ !(EU_wave|UE_wave|EE_wave), rawwgchangeEUE_wave := rawwgchange_wave]

DTseam[ , c("wageAtEU","wageAfterUE","EU_wave_first", "UE_wave_last"):=NULL]


#cleaning the stayers:-----------------------------------------------------------------
# wagechange wave =NA for large gains or losses that revert

DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_bad2 := (rawwgchange_wave>2|rawwgchange_wave<(-2))  & abs(next.waverawwg - last.waverawwg)<.1 &
	   	(last.lfstat_wave==1) & (next.lfstat_wave==1)] 

DTseam[ is.na(wagechange_wave_bad2)  , wagechange_wave_bad2 :=F] 

#wagechange between 2 0's:
DTseam[lfstat_wave>=2 & next.lfstat_wave>=2 & !(EU_wave==T|UE_wave==T|EE_wave==T) , wagechange_wave_bad2:= T] 

# wagechange wave =NA for job changers without a transition
DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & jobchng_wave==T , wagechange_wave_jcbad := T]
DTseam[ , next.job_wave := shift(job_wave,type="lead"),by=id]
DTseam[ , last.job_wave := shift(job_wave,type="lag"),by=id]
DTseam[ lfstat_wave==1 & next.lfstat==1, no_jobchng_wave:= (job_wave==next.job_wave)&(job_wave==last.job_wave)]
# DTseam[ EE_wave==T & no_jobchng_wave==T, wagechange_wave_jcbad :=T]
DTseam[is.na(wagechange_wave_jcbad )==T , wagechange_wave_jcbad := F]
DTseam[ , c("next.job_wave","last.job_wave"):=NULL]

#wagechanges in the crazy 2004 months:
DTseam[ wave>=7 & panel==2004, wagechange_wave_2004bad :=T]

#lowest/highest wages out:
lowwageqtls= DTseam[ lfstat_wave==1, quantile(waverawwg, na.rm = T, probs=c(.01,.02,.05))]
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_low :=waverawwg<lowwageqtls[2] | next.waverawwg<lowwageqtls[2] ]
highwageqtls= DTseam[ lfstat_wave==1, quantile(waverawwg, na.rm = T, probs=c(.95,.98,.99))]
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_high :=waverawwg>highwageqtls[2] | next.waverawwg>highwageqtls[2] ]

# take out the early attrition
DTseam[ , panelmaxmis:= max(maxmis,na.rm = T), by=panel]
DTseam[ , pctmaxmis:= maxmis/panelmaxmis]

# take out imputed earnings:
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T) & (earn_imp_wave>=1 | shift(earn_imp_wave,type="lead")>=1), wagechange_wave_imp := T, by=id]
DTseam[ is.na(wagechange_wave_imp)==T, wagechange_wave_imp := F]


DTseam<-subset(DTseam, select = c("wagechange_wave","wagechange_wave_jcbad","wagechange_wave_bad2","wagechange_wave_low","wagechange_wave_high","wagechange_wave_imp","pctmaxmis"
								  ,"wagechangeEUE_wave","next.wavewage","last.wagechange_wave","next.wagechange_wave",
								  "rawwgchangeEUE_wave","rawwgchange_wave","id","wave"))
DTall<- merge(DTall,DTseam,by=c("id","wave"),all.x=T)

#looking at occupation-level wage changes
DTall[ , occwage_wave:= sum(1/2*(exp(occwage)-exp(-occwage)),na.rm=T) , by=list(id,wave)]
DTall[ , occwage_wave:= log(occwage_wave + (1+occwage_wave^2)^.5)]
DTall[ , nawavewage:= all(is.na(occwage) ),by= list(id,wave)]
DTall[ nawavewage==T, occwage_wave:=NA_real_]

DTseam <- DTall[ seam==T]
DTseam[ , next.occwage_wave := shift(occwage_wave,type="lead"), by=id]
DTseam[ switchedOcc_wave==T, occwagechange_wave:= next.occwage_wave-occwage_wave]
DTseam[ switchedOcc_wave==F, occwagechange_wave:= 0.]

DTseam[ EU_wave==T , EU_wave_first := is.na(shift(EU_wave)), by=list(id,ustintid_wave)]
DTseam[ UE_wave==T , UE_wave_last := is.na(shift(UE_wave,type="lead")), by=list(id,ustintid_wave)]

DTseam[EU_wave_first == T, occwageAtEU := occwage_wave]
DTseam[, occwageAtEU := na.locf(occwageAtEU, na.rm = F),by=list(ustintid_wave, id)]
DTseam[UE_wave_last == T, occwageAfterUE :=  next.occwage_wave]
DTseam[UE_wave_last == T, occwagechangeEUE_wave := occwageAfterUE - occwageAtEU]
DTseam[, occwagechangeEUE_wave:= Mode(occwagechangeEUE_wave), by=list(ustintid_wave, id)]
DTseam[ EE_wave==T, occwagechangeEUE_wave:= occwagechange_wave]
DTseam[ switchedOcc_wave==F, occwagechangeEUE_wave:= 0.]

DTseam<-subset(DTseam, select = c("occwagechange_wave","occwagechangeEUE_wave","next.occwage_wave","id","wave"))
DTall<- merge(DTall,DTseam,by=c("id","wave"),all.x=T)

saveRDS(DTall, paste0(datadir,"/DTall_5.RData"))

#wc_wave <- DTall[seam==T & lfstat_wave==1 & next.lfstat_wave==1 & wagechange_wave_bad==F & wc_wave<0.2, .(wc_wave = weighted.mean(wagechange_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
#ggplot(wc_wave, aes(date, wc_wave, color = panel, group = panel)) +
#	geom_point() +
#	geom_line() +xlab("") + ylab("mean wage change, stayers wave-frequency")

