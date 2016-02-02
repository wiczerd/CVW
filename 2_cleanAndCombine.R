# December 18, 2015
# Clean data
# 1) restrict sample by age
# 2) remove observations with no employment status information
# 3) recode esr into lfstat, remove esr
# 4) remove suspect earnings jumps
# 5) fix occupation code to be consistent with employment status
# 6) fix job code to be consistent with employment status
# 7) fill in missing occupation codes with next observed occupation
# 8) combine panels into one dataset
# 9) save intermediate result, DTall_2.RData
library(data.table)
library(zoo)

# Q: Why don't we combine panels first, then do all of the cleaning just once?
# A: The dataset might get too big. When the datasets are combined after cleaning,
#    they are much smaller than the original sets. It's worth a try, though.

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)


# 1996 panel --------------------------------------------------------------

DT96 <- readRDS("./Data/DT96.RData")

# restrict sample by age
DT96 <- DT96[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT96 <- DT96[!is.na(esr),]

# recode esr into lfstat
DT96[esr >= 1 & esr <= 5, lfstat := 1]
DT96[esr == 6 | esr == 7, lfstat := 2]
DT96[esr == 8, lfstat := 3]
DT96[, esr := NULL]

# fix occupation code to consistent with employment status, replace soc2d with occ
DT96[, occ := soc2d]
DT96[, soc2d := NULL]
setkey(DT96, id, date)
DT96[lfstat == 2 | lfstat == 3, occ := NA_integer_]
DT96[, occ := na.locf(occ, na.rm = FALSE, fromLast = TRUE), by = id]
DT96 <- DT96[!is.na(occ),]

# fix job code to be consistent with employment status
DT96[is.na(job), job := 0]
DT96[lfstat == 2 | lfstat == 3, job := 0]

#clean out call-backs, where job-id is the same, replace them with NA lfstat
DT96[lfstat==1, jobcheck:= job, by = id]
DT96[ , jobcheck := na.locf( jobcheck, na.rm=F), by = id]
DT96[ lfstat!= 1 & jobcheck== shift(jobcheck,1,type="lead"), lfstat := NA_real_ , by = id]


saveRDS(DT96, "./Data/DT96_2.RData")
rm(DT96)

# 2001 panel --------------------------------------------------------------

DT01 <- readRDS("./Data/DT01.RData")

# restrict sample by age
DT01 <- DT01[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT01 <- DT01[!is.na(esr),]

# recode esr into lfstat
DT01[esr >= 1 & esr <= 5, lfstat := 1]
DT01[esr == 6 | esr == 7, lfstat := 2]
DT01[esr == 8, lfstat := 3]
DT01[, esr := NULL]

# fix occupation code to consistent with employment status, replace soc2d with occ
DT01[, occ := soc2d]
DT01[, soc2d := NULL]
setkey(DT01, id, date)
DT01[lfstat == 2 | lfstat == 3, occ := NA_integer_]
DT01[, occ := na.locf(occ, na.rm = FALSE, fromLast = TRUE), by = id]
DT01 <- DT01[!is.na(occ),]

# fix job code to be consistent with employment status
DT01[is.na(job), job := 0]
DT01[lfstat == 2 | lfstat == 3, job := 0]

#clean out call-backs, where job-id is the same, replace them with NA lfstat
DT01[lfstat==1, jobcheck:= job, by = id]
DT01[ , jobcheck := na.locf( jobcheck, na.rm=F), by = id]
DT01[ lfstat!= 1 & jobcheck== shift(jobcheck,1,type="lead"), lfstat := NA_real_ , by = id]

saveRDS(DT01, "./Data/DT01_2.RData")
rm(DT01)

# 2004 panel --------------------------------------------------------------

DT04 <- readRDS("./Data/DT04.RData")

# restrict sample by age
DT04 <- DT04[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT04 <- DT04[!is.na(esr),]

# recode esr into lfstat
DT04[esr >= 1 & esr <= 5, lfstat := 1]
DT04[esr == 6 | esr == 7, lfstat := 2]
DT04[esr == 8, lfstat := 3]
DT04[, esr := NULL]

# fix occupation code to consistent with employment status, replace soc2d with occ
DT04[, occ := soc2d]
DT04[, soc2d := NULL]
setkey(DT04, id, date)
DT04[lfstat == 2 | lfstat == 3, occ := NA_integer_]
DT04[, occ := na.locf(occ, na.rm = FALSE, fromLast = TRUE), by = id]
DT04 <- DT04[!is.na(occ),]

# fix job code to be consistent with employment status
DT04[is.na(job), job := 0]
DT04[lfstat == 2 | lfstat == 3, job := 0]


#clean out call-backs, where job-id is the same, replace them with NA lfstat
DT04[lfstat==1, jobcheck:= job, by = id]
DT04[ , jobcheck := na.locf( jobcheck, na.rm=F), by = id]
DT04[ lfstat!= 1 & jobcheck== shift(jobcheck,1,type="lead"), lfstat := NA_real_ , by = id]


saveRDS(DT04, "./Data/DT04_2.RData")
rm(DT04)

# 2008 panel --------------------------------------------------------------

DT08 <- readRDS("./Data/DT08.RData")

# restrict sample by age
DT08 <- DT08[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT08 <- DT08[!is.na(esr),]

# recode esr into lfstat
DT08[esr >= 1 & esr <= 5, lfstat := 1]
DT08[esr == 6 | esr == 7, lfstat := 2]
DT08[esr == 8, lfstat := 3]
DT08[, esr := NULL]

# fix occupation code to consistent with employment status, replace soc2d with occ
DT08[, occ := soc2d]
DT08[, soc2d := NULL]
setkey(DT08, id, date)
DT08[lfstat == 2 | lfstat == 3, occ := NA_integer_]
DT08[, occ := na.locf(occ, na.rm = FALSE, fromLast = TRUE), by = id]
DT08 <- DT08[!is.na(occ),]

# fix job code to be consistent with employment status
DT08[is.na(job), job := 0]
DT08[lfstat == 2 | lfstat == 3, job := 0]

#clean out call-backs, where job-id is the same, replace them with NA lfstat
DT08[lfstat==1, jobcheck:= job, by = id]
DT08[ , jobcheck := na.locf( jobcheck, na.rm=F), by = id]
DT08[ lfstat!= 1 & jobcheck== shift(jobcheck,1,type="lead"), lfstat := NA_real_ , by = id]

saveRDS(DT08, "./Data/DT08_2.RData")
rm(DT08)

# Combine panels ----------------------------------------------------------

DT96 <- readRDS("./Data/DT96_2.RData")
DT01 <- readRDS("./Data/DT01_2.RData")
DT04 <- readRDS("./Data/DT04_2.RData")
DT08 <- readRDS("./Data/DT08_2.RData")
DTList <- list(DT96, DT01, DT04, DT08)
DTall <- rbindlist(DTList, use.names = TRUE)
setkey(DTall, id, date)
rm(list = c("DTList", "DT96", "DT01", "DT04", "DT08"))

saveRDS(DTall, "./Data/DTall_2.RData")
rm(list=ls())
