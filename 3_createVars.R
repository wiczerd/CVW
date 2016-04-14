# December 21, 2015
# Create labor force flow dummies
# 1) create switchedJob dummy
# 2) create switchedOcc dummy
# 3) create UE, EU, EE dummies
# 4) create unemployment duration variable
# 5) create Young and HSCol variables
# 6) save intermediate result, DTall_3.RData
library(data.table)
library(zoo)


#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

Mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}

DTall <- readRDS("./Data/DTall_2.RData")

setkey(DTall, id, date)

# create switchedJob dummy
#DTall[, switchedJob := (job != shift(job, 1, type = "lead")) &
#    	(shift(job, 1, type = "lag") != shift(job, 1, type = "lead")), by = id]
DTall[, switchedJob := job != shift(job, 1, type = "lead") , by = id]

# create switchedOcc dummy
DTall[, switchedOcc := (occ != shift(occ, 1, type = "lead")) &
     	(shift(occ, 1, type = "lag") != shift(occ, 1, type = "lead")) &
	  	(occ != shift(occ, 2, type = "lead")) &
     	switchedJob, by = id]
# create a switchedAddress dummy
DTall[, switchedAddress := (shhadid != shift(shhadid, 1, type = "lead")) &
	  	(shift(shhadid, 1, type = "lag") != shift(shhadid, 1, type = "lead")) &
	  	(shhadid != shift(shhadid, 2, type = "lead")) &
	  	switchedJob, by = id]
DTall[, switchedInd := (ind != shift(ind, 1, type = "lead")) &
	  	(shift(ind, 1, type = "lag") != shift(ind, 1, type = "lead")) &
	  	(ind != shift(ind, 2, type = "lead")) &
	  	switchedJob, by = id]

# create EE, EU, and UE dummies
DTall[, EE := lfstat == 1 & shift(lfstat, 1, type = "lead") == 1 & switchedJob==T, by = id]
DTall[, EU := lfstat == 1 & shift(lfstat, 1, type = "lead") == 2 & switchedJob==T, by = id]
DTall[, UE := lfstat == 2 & shift(lfstat, 1, type = "lead") == 1 & switchedJob==T, by = id]
# take stint ID into EU and clean it:
DTall[, fstintid:= shift(stintid, 1, type = "lead"), by = id]
DTall[EU==T, stintid := fstintid, by=id]
DTall[, fstintid := NULL]

#check
DTall[UE==T, dupedUE:= duplicated(stintid, na.rm=T), by=id]
DTall[EU==T, dupedEU:= duplicated(stintid, na.rm=T), by=id]
#drop the duplicates:
DTall[ UE==T & dupedUE, UE:=F ]
DTall[ EU==T & dupedEU, EU:=F ]
# it is still possible there are lfstat=2|3 that are associated with a duplicate
DTall[ , c("dupedEU","dupedUE"):=NULL]

# create unemployment duration variable
DTall[, unempdur := seq_len(.N)-1, by = list(id, stintid)]
DTall[is.na(stintid), unempdur := NA]
DTall[!is.na(unempdur) , maxunempdur := max(unempdur, na.rm=T), by= list(id,stintid)]

# drop bad earnings data
# Q: Why is this not in step 2?
# A: It requires UE and UE, which aren't defined yet in step 2.
DTall[, nomearnm := earnm]
DTall[, earnm := earnm/PCEPI*100]
DTall[, badearn := abs(log(shift(earnm, 1, type = "lead")/earnm)) > 2.0 & 
      	abs(log(shift(earnm, 1, type = "lead")/shift(earnm, 1, type = "lag"))) < 0.1, by = id]
DTall[UE | EU | EE, badearn := FALSE]
DTall[(badearn), earnm := NA_real_]
DTall[, c("badearn", "nomearnm") := NULL]

# create Young and HScol variables
DTall[, Young := age < 30]
DTall[, HSCol := (educ >= 4) + (educ >= 2)]


# compute seam-to-seam status change variable ----------------------------
DTall[ , seam:= wave != shift(wave,1,type="lead"), by=id ]

DTall[ is.finite(lfstat), lfstat_wave := max(lfstat,na.rm=T), by=list(id,wave)]
DTall[ is.finite(lfstat), lfstat2_wave := any(lfstat==2, na.rm=T), by=list(id,wave)] #anytime an unemployment shows up, we count it as such
DTall[ is.finite(lfstat), lfstat_wave := ifelse(lfstat2_wave, 2L, lfstat_wave)] #anytime an unemployment shows up, we count it as such

DTall[ , UE_wave := (lfstat_wave==2L & shift(lfstat_wave,1,type="lead")==1L), by=id] #will only count 1 if seam ==1
DTall[ , UE_wave := any(UE_wave, na.rm = T), by=list(id,wave)] 
DTall[ , EU_wave := (lfstat_wave==1L & shift(lfstat_wave,1,type="lead")==2L), by=id] #will only count 1 if seam ==1
DTall[ , EU_wave := any(EU_wave, na.rm = T), by=list(id,wave)] 
DTall[ is.na(UE_wave), UE_wave:=F]
DTall[ is.na(EU_wave), EU_wave:=F]

DTall[ is.finite(EE), EEmax_wave := any(EE==T & seam==F,na.rm=T), by=list(id,wave)] #get whether there was an EE w/in the wave (and not on the edge)
DTall[ seam==T, EEnextmax_wave := shift(EEmax_wave,1,type="lead"), by = id]
DTall[ is.na(EEnextmax_wave), EEnextmax_wave := F]
DTall[ , EEnextmax_wave := any(EEnextmax_wave, na.rm=T), by=list(id,wave)]
DTall[ seam==T, EE_wave := EE]
DTall[ , EE_wave := (EE_wave | EEnextmax_wave )]
DTall[ is.na(EE_wave), EE_wave:=F]
DTall[, EE_wave := any(EE_wave, na.rm=T), by=list(id,wave)]
# HAVE THIS???
DTall[EEmaxwave ==T & EE_wave ==F, EE_wave := NA]
# only EE if nothing else
DTall[UE_wave==T & EE_wave==T, EE_wave :=F]
DTall[EU_wave==T & EE_wave==T, EE_wave :=F]

DTall[ , sOmax_wave := any(switchedOcc), by=list(id,wave)]
DTall[ seam==T, sOnextmax_wave := shift(sOmax_wave,1,type="lead"), by=id]
DTall[ , sOnextmax_wave := any(sOnextmax_wave), by=list(id,wave)]
DTall[ , switchedOcc_wave := (EEnextmax_wave & sOnextmax_wave) | (EE & switchedOcc) | (EU_wave & sOmax_wave) , by=id]
DTall[ , switchedOcc_wave := any(switchedOcc_wave), by=list(id,wave)]

DTall[, c("EEnextmax_wave","EEmax_wave","lfstat2_wave","sOmax_wave","sOnextmax_wave"):= NULL] 

DTall[, recIndic_wave := any(recIndic, na.rm=T), by=list(wave,id)]

#some diagnostics -------------------------------------------
sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE & DTall$stintid>0 & is.finite(DTall$stintid)], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat==2 & is.finite(DTall$stintid)], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$EE & DTall$switchedInd], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedInd], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[!is.na(DTall$unempdur)]*DTall$unempdur[!is.na(DTall$unempdur)],na.rm=T)/sum(DTall$wpfinwgt[!is.na(DTall$unempdur)],na.rm=T)


saveRDS(DTall, "./Data/DTall_3.RData")
rm(list=ls())
