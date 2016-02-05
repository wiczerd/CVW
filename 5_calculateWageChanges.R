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


DTall <- readRDS("./Data/DTall_4.RData")

setkey(DTall, id, date)

# fill wages upwards to fill in missing observations
DTall[, EmpTmrw := EE | UE, by = id]
DTall[EmpTmrw == T, nextwage := shift(usewage, 1, type = "lead"), by = id]
DTall[EmpTmrw == T & is.na(nextwage), nextwage := -9e16, by = id]
DTall[, nextwage := na.locf(nextwage, na.rm = FALSE, fromLast = TRUE), by = id]
DTall[EmpTmrw == T, nextoccwage := shift(occwage, 1, type = "lead"), by = id]
DTall[EmpTmrw == T & is.na(nextoccwage), nextoccwage := -9e16, by = id]
DTall[, nextoccwage := na.locf(nextoccwage, na.rm = FALSE, fromLast = TRUE), by = id]
DTall[nextwage <=-8e16, nextwage:=NA_real_]
DTall[nextoccwage <=-8e16, nextoccwage:=NA_real_]

DTall[, tuw := shift(usewage, 1, type = "lag"), by = id]
DTall[!is.finite(tuw) & is.finite(usewage) & lfstat == 1, tuw := usewage]

# create wagechange variable
DTall[EE == T, wagechange := nextwage - tuw]
DTall[EU == T, wagechange := log(1.0) - tuw]
DTall[UE == T, wagechange := nextwage - log(1.0)]

# create wagechange_stayer variable
DTall[job == shift(job, 1, type = "lead") & job == shift(job, 1, type = "lag"),
      wagechange_stayer := shift(usewage, 1, type = "lead") - shift(usewage, 1, type = "lag"), by = id]
DTall[job == 0 & shift(job, 1, type = "lead") == 0,
      wagechange_stayer := 0.0, by = id]

# create wagechange_EUE variable
DTall[EU | EE, wagechange_EUE := nextwage - tuw]
DTall[!(EU | EE | UE), wagechange_EUE := wagechange_stayer]

# create wagechange_all variable
# if ifelse() condition is NA, end result is NA.
DTall[, wagechange_all := wagechange]
DTall[job == shift(job, 1, type = "lead") & job > 0, wagechange_all := wagechange_stayer, by = id]
DTall[(EE | UE | EU), wagechange_all := wagechange]

# create occwagechange variable
DTall[, occwagechange := as.numeric(ifelse(switchedOcc & !shift(switchedOcc, 1, type = "lag"), 
				nextoccwage - shift(occwage, 1, type = "lag"), 
				0.)), by = id]
DTall[, occwagechange := as.numeric(ifelse(switchedOcc & shift(switchedOcc, 1, type = "lag"),
				nextoccwage - occwage, 
				0.)), by = id]

saveRDS(DTall, "./Data/DTall_5.RData")
rm(list=ls())
