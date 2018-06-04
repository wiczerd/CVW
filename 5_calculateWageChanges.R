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
DTall <- DTall[ , c("eyear","next.earnm","logearnm","black","hispanic","filtered.unrate$cycle"):=NULL]
#I don't want to drop these ->  DTall <- DTall[ lfstat_wave<3, ]

setkey(DTall, id, date)
#clear out the super-low earnings
DTall[ lfstat==1 & usewage<log(20+(1+20^2)^.5), usewage:=NA]
DTall[ , calmon := (wave-1)*4+srefmon]
DTall[ calmon+1 == shift(calmon,type="lead"), next.usewage := shift(usewage,type="lead"), by=id]
DTall[ calmon-1 == shift(calmon,type="lag" ), last.usewage := shift(usewage,type="lag"), by=id]

# fill wages upwards to fill in unemployment spells
DTall[EE==T|UE == T, nextwage := next.usewage]
DTall[lfstat>=2 | EU==T, nextwage := Mode(nextwage), by = list(id,ustintid)] #replace within  unemp stints

DTall[ , last.lfstat:= shift(lfstat), by=id]
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
DTall[ , c("wagechange_stayer","wagechange","wagechange_stayer_bad","tuw","last.lfstat"):=NULL]

# create occwagechange variable
#DTall[, leadwage := shift(occwage, 1, type = "lead"), by = id]
#DTall[, last.occwage := shift(occwage, 1, type = "lag"), by = id]
#DTall[EE==T|UE == T, nextoccwage := leadwage, by = id]
#DTall[lfstat>=2 | EU==T, nextoccwage := Mode(nextoccwage), by = list(id,ustintid)] #replace if it's UE
#DTall[, leadwage:=NULL]
# DTall[switchedOcc_wave ==T, occwagechange := nextoccwage - last.occwage]
# DTall[ , c("nextoccwage","last.occwage"):=NULL]

#***********************************************************************
#***********************************************************************
# compute seam-to-seam and year-to-year wage change variable ----------------------------
#***********************************************************************
#***********************************************************************

# residual wages
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]

DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
DTall[ , wavewage := mean(levwage,na.rm=T)*nmo_lf1, by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , wavewage := log(wavewage + (1+wavewage^2)^.5) ]
DTall[ , nawavewage:= all(is.na(usewage) ),by= list(id,wave)]
DTall[ nawavewage==T, wavewage:=NA_real_]
#drop the lowest resid wages, implies working less than $80 /month:
DTall[ lfstat_wave==1 & wavewage<log(80+(1+80^2)^.5), wavewage:=NA]

DTall[ , c("levwage","nawavewage"):=NULL]


# raw wages
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
DTall[ , waverawwg := mean(earnm,na.rm=T)*nmo_lf1, by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , waverawwg := log(waverawwg + (1+waverawwg^2)^.5) ]
DTall[ , nawavewage:= all(is.na(usewage) ),by= list(id,wave)]
DTall[ nawavewage==T, waverawwg:=NA_real_]
#drop the lowest wages, implies working less than $80 /month:
DTall[ lfstat_wave==1 & waverawwg<log(80+(1+80^2)^.5), waverawwg:=NA]
DTall[ lfstat_wave==1, earn_imp_wave := sum(earn_imp==1,na.rm=T), by=list(id,wave)]

DTall[ , c("nawavewage","nmo_lf1"):=NULL]

#************************************************************************************
#wave frequency changes --------------------------------------------
DTseam <- DTall[ seam==T,]
#need to add change across waves (use wavewage)
DTseam[ wave+1==shift(wave,type = "lead"), next.wavewage := shift(wavewage,1,type="lead"), by=id]
DTseam[ wave-1==shift(wave,type = "lag" ), last.wavewage := shift(wavewage,1,type="lag"), by=id]

DTseam[ , nw:= next.wavewage]
DTseam[ , tw:= wavewage]

DTseam[ , wagechange_wave := nw - tw]
#for EE that spans waves:
DTseam[ EE_wave==T & EEmon<seammon , wagechange_wave:= NA]
DTseam[ EE_wave==T & shift(midEE)==T , wagechange_wave:= (next.wavewage - last.wavewage)]
DTseam[ , next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTseam[EE_wave==T & midEE ==T , wagechange_wave := next.wagechange_wave]

DTseam[ , next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTseam[ , last.wagechange_wave := shift(wagechange_wave, type="lag" ),by=id]

# create wave-level EUE wage change-----------------------------------------
# find wage in period before an EU and one period after UE
DTseam[ , EU_wave_first := EU_wave==T & !(shift(EU_wave)==T), by=id]
DTseam[ , UE_wave_last  := UE_wave==T & !(shift(UE_wave,type="lead")==T), by=id]

DTseam[ , last.lfstat_wave:= shift(lfstat_wave), by=id]
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

DTseam[ wave+1==shift(wave,type = "lead"), next.waverawwg := shift(waverawwg,1,type="lead"), by=id]
DTseam[ wave-1==shift(wave,type = "lag" ), last.waverawwg := shift(waverawwg,1,type="lag"), by=id]

DTseam[ , nw:= next.waverawwg]
DTseam[ , tw:= waverawwg]

DTseam[ , rawwgchange_wave := nw - tw]
#for EE that spans waves:
DTseam[ EE_wave==T & EEmon<seammon , rawwgchange_wave:= NA]
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
DTseam[next.lfstat_wave==1 & UE_wave_last == T, wageAfterUE :=  next.waverawwg]
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

DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & jobchng_wave==T , wagechange_wave_jcbad := T]
DTseam[ , next.job_wave := shift(job_wave,type="lead"),by=id]
DTseam[ , last.job_wave := shift(job_wave,type="lag"),by=id]
DTseam[ lfstat_wave==1 & next.lfstat==1, no_jobchng_wave:= (job_wave==next.job_wave)&(job_wave==last.job_wave)]
# DTseam[ EE_wave==T & no_jobchng_wave==T, wagechange_wave_jcbad :=T]
DTseam[is.na(wagechange_wave_jcbad )==T , wagechange_wave_jcbad := F]
DTseam[ , c("next.job_wave","last.job_wave"):=NULL]

#wagechanges in the crazy 2004 months:
# ? DTseam[ wave>=7 & panel==2004, wagechange_wave_2004bad :=T]

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
								  ,"wagechangeEUE_wave","next.wavewage","rawwgchangeEUE_wave","rawwgchange_wave","id","wave"))
DTall<- merge(DTall,DTseam,by=c("id","wave"),all.x=T)

#looking at occupation-level wage changes
DTall[ , occwage_wave:= sum(1/2*(exp(occwage)-exp(-occwage)),na.rm=T) , by=list(id,wave)] #sum in levels
DTall[ , occwage_wave:= log(occwage_wave + (1+occwage_wave^2)^.5)] #convert back to inv hyperbolic sine
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

DTseam<-subset(DTseam, select = c("occwagechange_wave","occwagechangeEUE_wave","id","wave"))
DTall<- merge(DTall,DTseam,by=c("id","wave"),all.x=T)

saveRDS(DTall, paste0(datadir,"/DTall_5.RData"))


#************************************************************************
#Annual wage changes ---------------------------
#************************************************************************

# residual wages
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]

DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,yri)]
DTall[ , annwage := mean(levwage,na.rm = T)*nmo_lf1, by=list(id,yri)]
DTall[ , annwage := log(annwage + (1+annwage^2)^0.5)]
DTall[ , naannwage:= all(is.na(usewage) ),by= list(id,yri)]
DTall[ naannwage==T, annwage:=NA_real_]
DTall[ , c("levwage","naannwage"):=NULL]

# raw wages
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,yri)]
DTall[ , annrawwg := mean(earnm,na.rm = T)*nmo_lf1, by=list(id,yri)]
DTall[ , annrawwg := log(annrawwg + (1+annrawwg^2)^0.5)]
DTall[ , naannwage:= all(is.na(usewage) ),by= list(id,yri)]
DTall[ naannwage==T, annrawwg:=NA_real_]
DTall[ , c("naannwage","nmo_lf1"):=NULL]


#merge together transitions
sipp_ann <- readRDS(paste0(datadir,"/sipp_ann.RData"))
DTall<- merge(DTall,sipp_ann,by=c("id","yri"),all.x=T)

DTall[ , yrseam := (yri !=shift(yri,type="lead") | mis==maxmis), by=id]
DTann <- DTall[ yrseam==T,]
#need to add change across waves (use wavewage)
DTann[ , next.annwage := shift(annwage,1,type="lead"), by=id]
DTann[ , next2.annwage := shift(annwage,2,type="lead"), by=id]
DTann[ , next3.annwage := shift(annwage,3,type="lead"), by=id]
DTann[ , last.annwage := shift(annwage,1,type="lag"), by=id]

DTann[ , wagechange_ann := next.annwage - annwage]
DTann[ , wagechange2_ann := next2.annwage - annwage]
DTann[ , wagechange3_ann := next3.annwage - annwage]

DTann[ , next.wagechange_ann := shift(wagechange_ann, type="lead"),by=id]
DTann[ , last.wagechange_ann := shift(wagechange_ann, type="lag" ),by=id]

DTann[ annwage>0     , rank.annwage := rank(annwage,na.last = "keep",ties.method="average"), by=year]
DTann[ annwage>0     , rank.annwage := (rank.annwage - min(rank.annwage,na.rm = T))/(max(rank.annwage,na.rm = T) - min(rank.annwage,na.rm = T)), by=year]
DTann[ annwage>0     , pctile.annwage := round(rank.annwage*100)]
DTann[ last.annwage>0, last.rank.annwage := rank(last.annwage,na.last = "keep",ties.method="average"), by=year]
DTann[ last.annwage>0, last.rank.annwage := (last.rank.annwage - min(last.rank.annwage,na.rm = T))/(max(last.rank.annwage,na.rm = T) - min(last.rank.annwage,na.rm = T)), by=year]
DTann[ last.annwage>0, last.pctile.annwage := round(last.rank.annwage*100)]

#raw annual wage changes:
DTann[ , next.annrawwg := shift(annrawwg,1,type="lead"),by=id]
DTann[ , next2.annrawwg := shift(annrawwg,2,type="lead"),by=id]
DTann[ , next3.annrawwg := shift(annrawwg,3,type="lead"),by=id]
DTann[ , last.annrawwg := shift(annrawwg,1,type="lag"),by=id]

DTann[ , rawwgchange_ann := next.annrawwg - annrawwg]
DTann[ , rawwgchange2_ann := next2.annrawwg - annrawwg]
DTann[ , rawwgchange3_ann := next3.annrawwg - annrawwg]

DTann[ , next.rawwgchange_ann := shift(rawwgchange_ann)]

DTann[ annrawwg>0     , rank.annrawwg := rank(annrawwg,na.last = "keep",ties.method="average"), by=year]
DTann[ annrawwg>0     , rank.annrawwg := (rank.annrawwg - min(rank.annrawwg,na.rm = T))/(max(rank.annrawwg,na.rm = T) - min(rank.annrawwg,na.rm = T)), by=year]
DTann[ annrawwg>0     , pctile.annrawwg := round(rank.annrawwg*100)]
DTann[ last.annrawwg>0, last.rank.annrawwg := rank(last.annrawwg,na.last = "keep",ties.method="average"), by=year]
DTann[ last.annrawwg>0, last.rank.annrawwg := (last.rank.annrawwg - min(last.rank.annrawwg,na.rm = T))/(max(last.rank.annrawwg,na.rm = T) - min(last.rank.annrawwg,na.rm = T)), by=year]
DTann[ last.annrawwg>0, last.pctile.annrawwg := round(last.rank.annrawwg*100)]

DTann<-subset(DTann, select = c("rawwgchange_ann","wagechange_ann","pctile.annrawwg","pctile.annwage","id","yri"))
DTall<- merge(DTall,DTann,by=c("id","yri"),all.x=T)


saveRDS(DTall, paste0(datadir,"/DTall_5_ann.RData"))

#wc_wave <- DTall[seam==T & lfstat_wave==1 & next.lfstat_wave==1 & wagechange_wave_bad==F & wc_wave<0.2, .(wc_wave = weighted.mean(wagechange_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
#ggplot(wc_wave, aes(date, wc_wave, color = panel, group = panel)) +
#	geom_point() +
#	geom_line() +xlab("") + ylab("mean wage change, stayers wave-frequency")

