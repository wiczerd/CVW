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

DTall[, linked := is.finite(lfstat) & is.finite( shift(lfstat, 1, type = "lead") ) & lfstat<=3 & shift(lfstat, 1, type = "lead")<=3, by=id]
# create switchedJob dummy
#DTall[, switchedJob := (job != shift(job, 1, type = "lead")) &
#    	(shift(job, 1, type = "lag") != shift(job, 1, type = "lead")), by = id]
DTall[ linked==T , switchedJob := job != shift(job, 1, type = "lead") , by = id]

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
DTall[linked==T, EE := lfstat == 1 & shift(lfstat, 1, type = "lead") == 1 & switchedJob==T, by = id]
DTall[linked==T, EU := lfstat == 1 & shift(lfstat, 1, type = "lead") == 2 & switchedJob==T, by = id]
DTall[linked==T, UE := lfstat == 2 & shift(lfstat, 1, type = "lead") == 1 & switchedJob==T, by = id]
#an alternative measure of EE
DTall[linked==T, EE_alt := shift(estlemp,type="lead")==F, by = id]
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
setkey(DTall, id,wave,date)
DTall[ , seam:= wave != shift(wave,1,type="lead"), by=id ]

DTall[, recIndic_wave := any(recIndic, na.rm=T), by=list(wave,id)]

DTall[ is.finite(lfstat), lfstat_wave := max(lfstat,na.rm=T), by=list(id,wave)]
DTall[ is.finite(lfstat), lfstat2_wave := any(lfstat==2, na.rm=T), by=list(id,wave)] #anytime an unemployment shows up, we count it as such
DTall[ is.finite(lfstat), lfstat_wave := ifelse(lfstat2_wave, 2L, lfstat_wave)] #anytime an unemployment shows up, we count it as such
DTall[ , lfstat2_wave:=NULL]

DTall[seam==T, linked_wave := is.finite(lfstat_wave) & is.finite( shift(lfstat_wave, 1, type = "lead") ) & lfstat_wave<=3 & shift(lfstat_wave, 1, type = "lead")<=3, by=id]
DTall[, linked_wave := any(linked_wave==T), by=list(id,wave)]


DTall[ linked_wave==T, UE_wave := (lfstat_wave==2L & shift(lfstat_wave,1,type="lead")==1L), by=id] #will only count 1 if seam ==1
DTall[ linked_wave==T, UE_wave := any(UE_wave, na.rm = T), by=list(id,wave)] 
DTall[ linked_wave==T, EU_wave := (lfstat_wave==1L & shift(lfstat_wave,1,type="lead")==2L), by=id] #will only count 1 if seam ==1
DTall[ linked_wave==T, EU_wave := any(EU_wave, na.rm = T), by=list(id,wave)] 
#DTall[ is.na(UE_wave), UE_wave:=F]
#DTall[ is.na(EU_wave), EU_wave:=F]

DTall[ linked_wave==T, EE_alt_nwave := any(shift(estlemp & seam==F,type="lead")==F,na.rm=T), by=list(id,wave)] #get whether there was an EE w/in the wave (and not on the edge)
DTall[ linked_wave==T & seam==T, EE_alt_wave := shift(EE_alt_nwave,type="lead"), by = id]
DTall[ linked_wave==T, EE_alt_wave := ifelse(seam==T & EE_alt_wave ==F, EE_alt, EE_alt_wave), by = id]
DTall[ linked_wave==T, EE_alt_wave := any(EE_alt_wave), by=list(id,wave)]

DTall[seam==T & linked_wave==T, EE_wave:= lfstat==1L & shift(lfstat_wave,1,type="lead")==1L & job != shift(job, 1, type = "lead"), by=id]
#DTall[ is.na(EE_wave) , EE_wave:=F]
DTall[ , EE_wave := any(EE_wave==T, na.rm=T), by=list(id,wave)]

DTall[ seam==T & linked_wave==T, EE_lastwave := lfstat_wave==1L & shift(lfstat,1,type="lag")==1L & job != shift(job, 1, type = "lag"), by=id]
DTall[ is.na(EE_lastwave)           , EE_lastwave:=F]
DTall[                              , EE_lastwave := any(EE_lastwave==T, na.rm=T), by=list(id,wave)]
DTall[ seam==T & linked_wave==T, EE_nextwave := shift(EE_wave,1,type="lead"), by=id]
DTall[ seam==T & linked_wave==T, EE_nextwave := any(EE_nextwave,na.rm = T), by=list(id,wave)]
DTall[ , EE_maxwave := any(EE==T, na.rm=T), by=list(id,wave)]

DTall[ seam==T & EE ==T, EE_nomatch := (EE_wave ==F) ]
DTall[ seam==F & EE ==T, EE_lastnomatch := (EE_lastwave ==F) ]
DTall[ , EE_lastnomatch := any(EE_lastnomatch==T , na.rm=T), by=list(id,wave)]
DTall[ seam==T , EE_nomatch := ifelse( shift(EE_lastnomatch ==T,1, type="lead"),T, EE_nomatch )]
DTall[ , EE_nomatch := any(EE_nomatch==T, na.rm=T), by=list(id,wave)]

DTall[ EE_nomatch==T & EE_wave !=T, EE_wave := NA]


# only EE if nothing else
DTall[UE_wave==T & EE_wave==T, EE_wave :=F]
DTall[EU_wave==T & EE_wave==T, EE_wave :=F]

#DTall[ , sOmax_wave := any(switchedOcc), by=list(id,wave)]
#DTall[ seam==T, sOnextmax_wave := shift(sOmax_wave,1,type="lead"), by=id]
#DTall[ , sOnextmax_wave := any(sOnextmax_wave), by=list(id,wave)]
#DTall[ , switchedOcc_wave := (EEnextmax_wave & sOnextmax_wave) | (EE & switchedOcc) | (EU_wave & sOmax_wave) , by=id]
#DTall[ , switchedOcc_wave := any(switchedOcc_wave,na.rm=T), by=list(id,wave)]

#DTall[ , EU_any := any(EU==T,na.rm=T), by=list(id,wave)]
#DTall[ , swOc_any := any(switchedOcc==T,na.rm=T), by=list(id,wave)]
#DTall[ , swOc_any := any(switchedOcc==T & EU==T,na.rm=T), by=list(id,wave)]
#DTall[ , swEUOcc_wave := ifelse(EU_wave==T & EU_any==T, swOc_any,shift(swOc_any,1,type="lead")), by=id]
#DTall[ EU_wave==T, switchedOcc_wave := (swEUOcc_wave|switchedOcc_wave)]

#DTall[ , switchedOcc_wave := ifelse(EU_wave==T,shift(switchedOcc_wave,1,type="lead"),switchedOcc_wave), by=list(id,wave)]
#DTall[, c("sOmax_wave","sOnextmax_wave"):= NULL] 


DTall[ , occ_mis1 := occ]
DTall[ wave==shift(wave), occ_mis1 := 0, by=id]
DTall[ , occ_mis1 := max(occ_mis1), by=list(id,wave)]
DTall[ , occ_mis4 := occ]
DTall[ !(seam==T), occ_mis4 := 0]
DTall[ , occ_mis4 := max(occ_mis4), by=list(id,wave)]
DTall[ seam==T, nextocc_mis4:= shift(occ_mis4,1,type="lead"), by=id]
DTall[ seam==T, switchedOcc_wave:= occ_mis1 != nextocc_mis4]

DTall[, c("nextocc_mis4","occ_mis4","occ_mis1"):= NULL] 

DTall[ lfstat==2, maxunempdur_wave := maxunempdur]
DTall[ lfstat_wave==2, maxunempdur_wave:=max(maxunempdur_wave,na.rm = T), by=list(id,wave)]
DTall[ , maxunempdur_wave := ifelse(EU_wave==T, shift(maxunempdur_wave,1,type="lead"),maxunempdur_wave)]
DTall[ EU_wave==T, maxunempdur_wave:=max(maxunempdur_wave,na.rm = T), by=list(id,wave)]


#some diagnostics -------------------------------------------
sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE & DTall$stintid>0 & is.finite(DTall$stintid)], na.rm=T)/
	sum(DTall$wpfinwgt[DTall$lfstat==2 & is.finite(DTall$stintid)], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$EE & DTall$switchedInd], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedInd], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[!is.na(DTall$unempdur)]*DTall$unempdur[!is.na(DTall$unempdur)],na.rm=T)/sum(DTall$wpfinwgt[!is.na(DTall$unempdur)],na.rm=T)


saveRDS(DTall, "./Data/DTall_3.RData")
rm(list=ls())
