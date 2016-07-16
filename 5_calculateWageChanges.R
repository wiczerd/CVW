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

setkey(DTall, id, date)
DTall[ , next.usewage := shift(usewage,type="lead"), by=id]
DTall[ , last.usewage := shift(usewage,type="lag"), by=id]
DTall[ , c("next.earnm","last.earnm") := NULL]

# fill wages upwards to fill in unemployment spells
DTall[EE==T|UE == T, nextwage := next.usewage]
DTall[, leadwage := shift(usewage, 2, type = "lead"), by = id]
DTall[, lead2status := next.lfstat==shift(lfstat,2,type="lead"), by=id]
DTall[(EE==T|UE == T )& is.na(nextwage) & lead2status==T, nextwage := leadwage, by=id] #replace with 2 leads forward
DTall[lfstat==2 | lfstat==3 | EU==T, nextwage := Mode(nextwage), by = list(id,ustintid)] #replace within  unemp stints
#DTall[lfstat==1 & !(EU==T), nextwage := Mode(nextwage), by = list(id,job)] #replace if it's EE

DTall[, leadwage := shift(occwage, 1, type = "lead"), by = id]
DTall[EE==T|UE == T, nextoccwage := leadwage, by = id]
DTall[, leadwage := shift(usewage, 2, type = "lead"), by = id]
DTall[, lead2status := shift(lfstat,1,type="lead")==shift(lfstat,2,type="lead"), by=id]
DTall[(EE==T|UE == T) & is.na(nextoccwage) & lead2status==T, nextoccwage := leadwage, by = id]
DTall[lfstat==2 | lfstat==3 | EU==T, nextoccwage := Mode(nextoccwage), by = list(id,ustintid)] #replace if it's UE
#DTall[lfstat==1 & !(EU==T), nextoccwage := Mode(nextoccwage), by = list(id,job)] #replace if it's EE

DTall[, c("leadwage","lead2status"):=NULL]

DTall[, tuw := ifelse(last.lfstat==1, last.usewage,usewage)]


# create wagechange variable---------------------------------------------
DTall[EE == T, wagechange := nextwage - tuw]
DTall[EU == T, wagechange := log(1.0) - tuw]
DTall[UE == T, wagechange := nextwage - log(1.0)]

# create wagechange_stayer variable------------------------------
DTall[ !(EE==T|EU==T|UE==T), wagechange_stayer := next.usewage - last.usewage]
DTall[job == 0 & shift(job,type="lead") == 0,
      wagechange_stayer := 0.0, by = id]

#do not allow stayers to change job listings without having a marked transition
DTall[ !(EE==T|EU==T|UE==T) & jobchng_wave ==T, wagechange_stayer:=NA]
#do not allow stayers to lose more than 200% change in earnings w/o change in status
DTall[ !(EE==T|EU==T|UE==T) & wagechange_stayer< -2., wagechange_stayer :=  NA ]

# create wagechange_EUE variable -----------------------------------
DTall[EU==T | EE==T, wagechange_EUE := nextwage - tuw]
#DTall[ EU==T|lfstat==2, wagechange_EUE:=ifelse(lfstat==2,
#							shift(wagechange_EUE,1,type="lag"),wagechange_EUE),by=id]
DTall[lfstat==2 | EU==T, wagechange_EUE := Mode(wagechange_EUE), by=list(id,ustintid)]
DTall[!(EU==T | EE==T | UE==T), wagechange_EUE := wagechange_stayer]

# create wagechange_all variable
# if ifelse() condition is NA, end result is NA.
DTall[, wagechange_all := wagechange]
DTall[job == shift(job, 1, type = "lead") & job > 0, wagechange_all := wagechange_stayer, by = id]
DTall[(EE | UE | EU), wagechange_all := wagechange]


# create occwagechange variable
#DTall[, occwagechange := as.numeric(ifelse(switchedOcc & !shift(switchedOcc, 1, type = "lag"), 
#				nextoccwage - shift(occwage, 1, type = "lag"), 
#				0.)), by = id]
#DTall[, occwagechange := as.numeric(ifelse(switchedOcc & shift(switchedOcc, 1, type = "lag"),
#				nextoccwage - occwage, 
#				0.)), by = id]

# compute seam-to-seam wage change variable ----------------------------
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]
DTall[ , wavewage := sum(levwage,na.rm=T), by= list(id,wave)]
DTall[ , wavewage := log(wavewage + (1+wavewage^2)^.5) ]
DTall[is.na(usewage)==T, wavewage:=NA_real_]
DTall[ , levwage:=NULL]
DTall[ seam==T, seamwage := usewage]
DTall[ , seamwage := Mode(seamwage), by=list(id,wave)]
#need to add change across waves (use wavewage)
DTall[ seam==T, wagechange_wave := shift(wavewage,1,type="lead") - wavewage, by=id]
DTall[ , wagechange_wave := Mode(wagechange_wave), by= list(id,wave)]
# wagechange wave =NA for job changers without a transition
DTall[!(EU_wave==T|UE_wave==T|EE_wave==T) & jobchng_wave==T , wagechange_wave := NA]
# wagechange wave =NA for gains that revert
DTall[ seam==T, next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTall[ seam==T, last.wagechange_wave := shift(wagechange_wave, type="lag" ),by=id]
DTall[ , next.wagechange_wave := Mode(next.wagechange_wave), by= list(id,wave)]
DTall[ , last.wagechange_wave := Mode(last.wagechange_wave), by= list(id,wave)]
DTall[!(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_bad:= (wagechange_wave>2 & (last.wagechange_wave< -2. | next.wagechange_wave< -2.)) ] #knocks out 42% of large increases and 1.4% of total change obseravations
# DTall[!(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_bad1:=(wagechange_wave< -2. & (last.wagechange_wave> 2. | next.wagechange_wave> 2.)) ] #knocks out 42% of large increases and 1.4% of total change obseravations
# wagechange wave =NA for loss more than 200% without a transition
DTall[!(EU_wave==T|UE_wave==T|EE_wave==T) & wagechange_wave< -2. , wagechange_wave := NA]
DTall[!(EU_wave==T|UE_wave==T|EE_wave==T) & wagechange_wave_bad == T , wagechange_wave := NA]
DTall[ , wagechange_wave_bad:=NULL]

#comparing seams:
DTall[ seam==T, wagechange_seam := shift(seamwage,1,type="lead") - seamwage, by =id]
DTall[ , wagechange_seam := Mode(wagechange_seam), by= list(id,wave)]

#looking at occupation-level wage changes
DTall[ , occwage_wave:= sum(1/2*(exp(occwage)-exp(-occwage)),na.rm=T) , by=list(id,wave)]
DTall[ , occwage_wave:= log(occwage_wave + (1+occwage_wave^2)^.5)]
DTall[ seam==T & switchedOcc_wave==T, occwagechange_wave:= shift(occwage_wave,type="lead")-occwage_wave,by=id]

saveRDS(DTall, paste0(datadir,"/DTall_5.RData"))
rm(list=ls())
