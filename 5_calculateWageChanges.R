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
setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

DTall <- readRDS("./Data/DTall_4.RData")

setkey(DTall, id, date)

# fill wages upwards to fill in unemployment spells
DTall[, EmpTmrw := EE | UE, by = id]
DTall[, leadwage := shift(usewage, 1, type = "lead"), by = id]
DTall[EmpTmrw == T, nextwage := leadwage]
DTall[, leadwage := shift(usewage, 2, type = "lead"), by = id]
DTall[, lead2status := shift(lfstat,1,type="lead")==shift(lfstat,2,type="lead"), by=id]
DTall[EmpTmrw == T & is.na(nextwage) & lead2status==T, nextwage := leadwage, by=id]
DTall[lfstat==2 | lfstat==3 | EU==T, nextwage := Mode(nextwage), by = list(id,stintid)] #replace within stint for unemp stints
DTall[lfstat==1 & !(EU==T), nextwage := Mode(nextwage), by = list(id,job)] #replace if it's EE

DTall[, leadwage := shift(occwage, 1, type = "lead"), by = id]
DTall[EmpTmrw == T, nextoccwage := leadwage, by = id]
DTall[, leadwage := shift(usewage, 2, type = "lead"), by = id]
DTall[, lead2status := shift(lfstat,1,type="lead")==shift(lfstat,2,type="lead"), by=id]
DTall[EmpTmrw == T & is.na(nextoccwage) & lead2status==T, nextoccwage := leadwage, by = id]
DTall[lfstat==2 | lfstat==3 | EU==T, nextoccwage := Mode(nextoccwage), by = list(id,stintid)] #replace if it's UE
DTall[lfstat==1 & !(EU==T), nextoccwage := Mode(nextoccwage), by = list(id,job)] #replace if it's EE

DTall[, c("leadwage","lead2status"):=NULL]


DTall[, tuw := as.numeric(ifelse(shift(lfstat)==1, shift(usewage, 1, type = "lag"),usewage)), by = id]
#DTall[(!is.finite(tuw) | tuw<=0) & is.finite(usewage) & lfstat == 1, tuw := usewage]
#DTall[, tuw := usewage, by = id]
#DTall[mis>0, tuw := as.numeric(ifelse(!is.finite(tuw) & lfstat == shift(lfstat), shift(usewage, 1, type = "lag"), tuw)), by = id]


# create wagechange variable---------------------------------------------
DTall[EE == T, wagechange := nextwage - tuw]
DTall[EU == T, wagechange := log(1.0) - tuw]
DTall[UE == T, wagechange := nextwage - log(1.0)]

# create wagechange_stayer variable------------------------------
DTall[job == shift(job, 1, type = "lead") & job == shift(job, 1, type = "lag"),
      wagechange_stayer := shift(usewage, 1, type = "lead") - shift(usewage, 1, type = "lag"), by = id]
DTall[job == 0 & shift(job, 1, type = "lead") == 0,
      wagechange_stayer := 0.0, by = id]

# create wagechange_EUE variable -----------------------------------
DTall[EU==T | EE==T, wagechange_EUE := nextwage - tuw]
#DTall[ EU==T|lfstat==2, wagechange_EUE:=ifelse(lfstat==2,
#							shift(wagechange_EUE,1,type="lag"),wagechange_EUE),by=id]
DTall[lfstat==2 | EU==T, wagechange_EUE := Mode(wagechange_EUE), by=list(id,stintid)]
DTall[!(EU==T | EE==T | UE==T), wagechange_EUE := wagechange_stayer]

# create wagechange_all variable
# if ifelse() condition is NA, end result is NA.
DTall[, wagechange_all := wagechange]
DTall[job == shift(job, 1, type = "lead") & job > 0, wagechange_all := wagechange_stayer, by = id]
DTall[(EE | UE | EU), wagechange_all := wagechange]

#do not allow stayers to lose more than 200% change in earnings w/o change in status
DTall[ EE==F&EU==F&UE==F, wagechange_all := ifelse( abs(wagechange_all)>2., NA,wagechange_all )]

# create occwagechange variable
DTall[, occwagechange := as.numeric(ifelse(switchedOcc & !shift(switchedOcc, 1, type = "lag"), 
				nextoccwage - shift(occwage, 1, type = "lag"), 
				0.)), by = id]
DTall[, occwagechange := as.numeric(ifelse(switchedOcc & shift(switchedOcc, 1, type = "lag"),
				nextoccwage - occwage, 
				0.)), by = id]

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
DTall[ seam==T, wagechange_seam := shift(seamwage,1,type="lead") - seamwage, by =id]
DTall[ , wagechange_seam := Mode(wagechange_seam), by= list(id,wave)]


saveRDS(DTall, "./Data/DTall_5.RData")
rm(list=ls())
