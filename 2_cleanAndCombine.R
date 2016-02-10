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
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}
esrRecode <- function(DF){
	# recode esr into lfstat
	DF[esr == 1 | esr == 2 | esr == 4, lfstat := as.integer(1)] #altenrative esr>=1 & esr<=5
	DF[esr == 3 | esr == 5 | esr == 6 | esr == 7, lfstat := as.integer(2)] #alternative: esr>=6 & esr<=7
	#DF[esr >= 1 & esr <= 5, lfstat := as.integer(1)] #altenrative esr>=1 & esr<=5
	#DF[esr == 6 | esr == 7, lfstat := as.integer(2)] #alternative: esr>=6 & esr<=7
	DF[esr == 8, lfstat := as.integer(3)]
}
estint <- function(DF){
	#sometimes job is not a unique identifier, i.e. if there's been a long nonemployment stint
	DF[lfstat == 1, newemp := mis==1, by=id]
	DF[is.na(newemp), newemp := F ]
	DF[, newemp := (lfstat == 1 & (shift(lfstat,1,type="lag") >=2)) | newemp ==T, by = id]
	DF[lfstat == 1, newemp := ((shift(job,1,type="lag") != job) | newemp==T) , by = id]
	DF[lfstat == 1 & is.na(newemp), newemp:= F ]
	DF[newemp==T, estintid := cumsum(newemp), by = id]
	DF[lfstat == 1, estintid := na.locf(estintid, na.rm = FALSE), by = id]
	DF[, newemp := NULL]
}


# 1996 panel --------------------------------------------------------------

DT96 <- readRDS("./Data/DT96_1.RData")

# restrict sample by age
DT96 <- DT96[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT96 <- DT96[!is.na(esr),]

# recode esr into lfstat
esrRecode(DT96)
# create time w/in id (month in sample)
DT96[, mis := seq_len(.N), by=id]

# replace soc2d with occ and make sure DT is sorted
DT96[, occ := soc2d]
DT96[, soc2d := NULL]
setkey(DT96, id, date)

#sometimes job is not a unique identifier, i.e. if there's been a long nonemployment stint
estint(DT96)

# replace occ with most common observed occ over employment spell
DT96[lfstat == 1, occ := Mode(occ), by = list(id, job,wave)]

# fix occupation code to be consistent with employment status
DT96[lfstat == 2 | lfstat == 3, occ := NA_integer_]

# create unemployment spell id
DT96[, newstint := lfstat == 2 & (shift(lfstat, 1, type = "lag") == 1 | mis == 1), by = id]
DT96[newstint == T, stintid := cumsum(newstint), by = id]
DT96[lfstat == 1, stintid := 0] #now only NA's are 2,3 that are connected to each other
DT96[, stintid := na.locf(stintid, na.rm = FALSE), by = id]
DT96[lfstat == 1, stintid := NA_integer_]

# fill in occupation with next occupation in unemployment stints
DT96[, leadocc := shift(occ, 1, type = "lead"), by = id]
DT96[lfstat==2|lfstat==3, occ := Mode(leadocc), by = list(id, stintid)]
DT96[, lagocc := shift(occ), by = id]
DT96[lfstat==2|lfstat==3, lagocc := Mode(lagocc), by = list(id, stintid)] 
DT96[, c("leadocc", "newstint") := NULL]

DT96 <- DT96[!is.na(occ),]

# fix job code to be consistent with employment status
DT96[is.na(job), job := 0]
DT96[lfstat == 2 | lfstat == 3, job := 0]

saveRDS(DT96, "./Data/DT96_2.RData")
rm(DT96)

# 2001 panel --------------------------------------------------------------

DT01 <- readRDS("./Data/DT01_1.RData")

# restrict sample by age
DT01 <- DT01[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT01 <- DT01[!is.na(esr),]

# recode esr into lfstat
esrRecode(DT01)
# create time w/in id (month in sample)
DT01[, mis := seq_len(.N), by=id]

# replace soc2d with occ and make sure DT is sorted
DT01[, occ := soc2d]
DT01[, soc2d := NULL]
setkey(DT01, id, date)

#sometimes job is not a unique identifier, i.e. if there's been a long nonemployment stint
estint(DT01)

# replace occ with most common observed occ over employment spell
DT01[lfstat == 1, occ := Mode(occ), by = list(id, job,wave)]

# fix occupation code to be consistent with employment status
DT01[lfstat == 2 | lfstat == 3, occ := NA_integer_]

# create unemployment spell id
DT01[, newstint := lfstat == 2 & (shift(lfstat, 1, type = "lag") == 1 | mis == 1), by = id]
DT01[(newstint), stintid := cumsum(newstint), by = id]
DT01[lfstat == 1, stintid := 0] #now only NA's are 2,3 that are connected to each other
DT01[, stintid := na.locf(stintid, na.rm = FALSE), by = id]
DT01[lfstat == 1, stintid := NA_integer_]

# fill in occupation with next occupation in unemployment stints
DT01[, leadocc := shift(occ, 1, type = "lead"), by = id]
DT01[lfstat==2|lfstat==3, occ := Mode(leadocc), by = list(id, stintid)]
DT01[, lagocc := shift(occ), by = id]
DT01[lfstat==2|lfstat==3, lagocc := Mode(lagocc), by = list(id, stintid)] 
DT01[, c("leadocc", "newstint") := NULL]

DT01 <- DT01[!is.na(occ),]

# fix job code to be consistent with employment status
DT01[is.na(job), job := 0]
DT01[lfstat == 2 | lfstat == 3, job := 0]


saveRDS(DT01, "./Data/DT01_2.RData")
rm(DT01)

# 2004 panel --------------------------------------------------------------

DT04 <- readRDS("./Data/DT04_1.RData")

# restrict sample by age
DT04 <- DT04[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT04 <- DT04[!is.na(esr),]

# recode esr into lfstat
esrRecode(DT04)
# create time w/in id (month in sample)
DT04[, mis := seq_len(.N), by=id]

# replace soc2d with occ and make sure DT is sorted
DT04[, occ := soc2d]
DT04[, soc2d := NULL]
setkey(DT04, id, date)

#sometimes job is not a unique identifier, i.e. if there's been a long nonemployment stint
estint(DT04)

# replace occ with most common observed occ over employment spell
DT04[lfstat == 1, occ := Mode(occ), by = list(id, job,wave)]

# fix occupation code to be consistent with employment status
DT04[lfstat == 2 | lfstat == 3, occ := NA_integer_]

# create unemployment spell id
DT04[, newstint := lfstat == 2 & (shift(lfstat, 1, type = "lag") == 1 | mis == 1), by = id]
DT04[(newstint), stintid := cumsum(newstint), by = id]
DT04[lfstat == 1, stintid := 0] #now only NA's are 2,3 that are connected to each other
DT04[, stintid := na.locf(stintid, na.rm = FALSE), by = id]
DT04[lfstat == 1, stintid := NA_integer_]

# fill in occupation with next occupation in unemployment stints
DT04[, leadocc := shift(occ, 1, type = "lead"), by = id]
DT04[lfstat==2|lfstat==3, occ := Mode(leadocc), by = list(id, stintid)]
DT04[, lagocc := shift(occ), by = id]
DT04[lfstat==2|lfstat==3, lagocc := Mode(lagocc), by = list(id, stintid)] 
DT04[, c("leadocc", "newstint") := NULL]

DT04 <- DT04[!is.na(occ),]

# fix job code to be consistent with employment status
DT04[is.na(job), job := 0]
DT04[lfstat == 2 | lfstat == 3, job := 0]

saveRDS(DT04, "./Data/DT04_2.RData")
rm(DT04)

# 2008 panel --------------------------------------------------------------

DT08 <- readRDS("./Data/DT08_1.RData")

# restrict sample by age
DT08 <- DT08[age >= 18 & age <= 65,]

# remove observations with no employment status information
DT08 <- DT08[!is.na(esr),]

# recode esr into lfstat
esrRecode(DT08)
# create time w/in id (month in sample)
DT08[, mis := seq_len(.N), by=id]

# replace soc2d with occ and make sure DT is sorted
DT08[, occ := soc2d]
DT08[, soc2d := NULL]
setkey(DT08, id, date)

#sometimes job is not a unique identifier, i.e. if there's been a long nonemployment stint
estint(DT08)

# replace occ with most common observed occ over employment spell
DT08[lfstat == 1, occ := Mode(occ), by = list(id, job,wave)]

# fix occupation code to be consistent with employment status
DT08[lfstat == 2 | lfstat == 3, occ := NA_integer_]

# create unemployment spell id
DT08[, newstint := lfstat == 2 & (shift(lfstat, 1, type = "lag") == 1 | mis == 1), by = id]
DT08[(newstint), stintid := cumsum(newstint), by = id]
DT08[lfstat == 1, stintid := 0] #now only NA's are 2,3 that are connected to each other
DT08[, stintid := na.locf(stintid, na.rm = FALSE), by = id]
DT08[lfstat == 1, stintid := NA_integer_]

# fill in occupation with next occupation in unemployment stints
DT08[, leadocc := shift(occ, 1, type = "lead"), by = id]
DT08[lfstat==2|lfstat==3, occ := Mode(leadocc), by = list(id, stintid)]
DT08[, lagocc := shift(occ), by = id]
DT08[lfstat==2|lfstat==3, lagocc := Mode(lagocc), by = list(id, stintid)] 
DT08[, c("leadocc", "newstint") := NULL]


DT08 <- DT08[!is.na(occ),]

# fix job code to be consistent with employment status
DT08[is.na(job), job := 0]
DT08[lfstat == 2 | lfstat == 3, job := 0]


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
