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
dupsUE<- DTall[UE==T, duplicated(stintid, na.rm=T), by=id]
dupsEU<- DTall[EU==T, duplicated(stintid, na.rm=T), by=id]


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

#some diagnostics
sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat ==1], na.rm=T)
sum(DTall$wpfinwgt[DTall$UE & DTall$stintid>0 & is.finite(DTall$stintid)], na.rm=T)/sum(DTall$wpfinwgt[DTall$lfstat==2 & is.finite(DTall$stintid)], na.rm=T)

sum(DTall$wpfinwgt[DTall$EE & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EE], na.rm=T)
sum(DTall$wpfinwgt[DTall$EU & DTall$switchedOcc], na.rm=T)/sum(DTall$wpfinwgt[DTall$EU], na.rm=T)
sum(DTall$wpfinwgt[!is.na(DTall$unempdur)]*DTall$unempdur[!is.na(DTall$unempdur)],na.rm=T)/sum(DTall$wpfinwgt[!is.na(DTall$unempdur)],na.rm=T)

saveRDS(DTall, "./Data/DTall_3.RData")
rm(list=ls())
