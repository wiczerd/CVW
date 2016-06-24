# June 22, 2016
# Daniel Eubanks
# Read in SIPP data and prepare it for analysis
library(foreign)
library(data.table)
library(zoo)
library(ggplot2)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

#rootdir <- "G:/Research_Analyst/Eubanks/Occupation Switching/SIPP to go"
rootdir <- "~/workspace/CVW/R/"
datadir <- paste(rootdir, "/InputDataDE", sep = "")
outputdir <- paste(rootdir, "/Results", sep = "")
setwd(rootdir)

############################## Read data, combine panels ##############################

# create list of variables to extract from the raw data
toKeep <- c("year", 
	    "month",
	    # demographics
	    "age",
	    "educ",
	    "female",
	    "race",
	    # technical and weight variables
	    "id",
	    "wave",
	    "srefmon",
	    "wpfinwgt",
	    # labor force status variables
	    "esr",
	    # job variables
	    "job",
	    "occ",
	    "earnm")

########## read in individual panels, extract variables, and subset sample

setwd(datadir)

# 1996 panel
sipp96 <- read.dta("./sippsets96ABD.dta", convert.factors = FALSE)
sipp96 <- sipp96[toKeep]
sipp96 <- subset(sipp96, !is.na(esr) & (age >= 16))
sipp96 <- data.table(sipp96)
sipp96$panel <- 1996

# 2001 panel
sipp01 <- read.dta("./sippsets01ABD.dta", convert.factors = FALSE)
sipp01 <- sipp01[toKeep]
sipp01 <- subset(sipp01, !is.na(esr) & (age >= 16))
sipp01 <- data.table(sipp01)
sipp01$panel <- 2001

# 2004 panel
sipp04 <- read.dta("./sippsets04ABD.dta", convert.factors = FALSE)
sipp04 <- sipp04[toKeep]
sipp04 <- subset(sipp04, !is.na(esr) & (age >= 16))
sipp04 <- data.table(sipp04)
sipp04$panel <- 2004

# 2008 panel
sipp08 <- read.dta("./sippsets08ABD.dta", convert.factors = FALSE)
sipp08 <- sipp08[toKeep]
sipp08 <- subset(sipp08, !is.na(esr) & (age >= 16))
sipp08 <- data.table(sipp08)
sipp08$panel <- 2008

rm(toKeep)


########## add soc2d codes using crosswalks

setwd(datadir)
occ1990_soc2d <- readRDS("./occ90_soc2d.RData")
coc2000_occ1990 <- readRDS("./coc2000_occ1990.RData")
coc1990_occ1990 <- readRDS("./coc1990_occ1990.RData")

# 1996 panel
sipp96[, coc90 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp96 <- merge(sipp96, coc1990_occ1990, by  = "coc90", all.x = TRUE)
sipp96 <- merge(sipp96, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp96[, c("occ", "coc90") := NULL]
setnames(sipp96, "occ1990", "occ")

# 2001 panel
sipp01[, coc90 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp01 <- merge(sipp01, coc1990_occ1990, by  = "coc90", all.x = TRUE)
sipp01 <- merge(sipp01, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp01[, c("occ", "coc90") := NULL]
setnames(sipp01, "occ1990", "occ")

# 2004 panel
sipp04[, coc2000 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp04 <- merge(sipp04, coc2000_occ1990, by  = "coc2000", all.x = TRUE)
sipp04 <- merge(sipp04, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp04[, c("occ", "coc2000") := NULL]
setnames(sipp04, "occ1990", "occ")

# 2008 panel
sipp08[, coc2000 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp08 <- merge(sipp08, coc2000_occ1990, by  = "coc2000", all.x = TRUE)
sipp08 <- merge(sipp08, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp08[, c("occ", "coc2000") := NULL]
setnames(sipp08, "occ1990", "occ")

rm(list=c("coc1990_occ1990", "coc2000_occ1990", "occ1990_soc2d"))


########## combine panels

sipp <- rbind(sipp96, sipp01, sipp04, sipp08)
rm(list=c("sipp96", "sipp01", "sipp04", "sipp08"))
sipp <- data.table(sipp)
sipp$panel <- as.factor(sipp$panel)

# save intermediate result
setwd(outputdir)
saveRDS(sipp, "./sipp.RData")

################################## Prepare data ###################################

# load intermediate result if starting from here
if(!exists("sipp")) {
	setwd(outputdir)
	sipp <- readRDS("./sipp.RData")
}

########## date, unique id, months in sample

# add date
sipp <- sipp[!is.na(month) & !is.na(year),]
sipp[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]

# create unique ID across all panels
sipp[, idpanel := paste(id, panel, sep = "")]

# sort by idpanel and date
setkey(sipp, idpanel, date)

# create month in sample variable
sipp[, mis := seq_len(.N), by = idpanel]

########## labor force status

# recode esr
sipp[esr <= 3, lfstat := 1]
sipp[esr >= 4 & esr <= 7, lfstat := 2]
sipp[esr == 8, lfstat := 3]

# create lag/lead lfstat
sipp[, last.lfstat := shift(lfstat, 1, type = "lag"), by = idpanel]
sipp[, next.lfstat := shift(lfstat, 1, type = "lead"), by = idpanel]

########## occupation and job

# replace occ with soc2d (occ will now refer to soc2d)
sipp[, occ := soc2d]
sipp[, soc2d := NULL]

# make sure job and occ are NA if not employed
sipp[lfstat != 1, job := NA]
sipp[lfstat != 1, occ := NA]

# replace occ with mode occ over a job spell
sipp[lfstat == 1, occ := Mode(occ), by = list(idpanel, job)]

# create lag/lead job
sipp[, last.job := shift(job, 1, type = "lag"), by = idpanel]
sipp[, next.job := shift(job, 1, type = "lead"), by = idpanel]

########## stint ids

# create an employment stint id
sipp[, newemp := lfstat == 1 & (last.lfstat >= 2 | mis == 1 | job != last.job)]
sipp[lfstat == 1 & is.na(newemp), newemp := FALSE]
sipp[newemp == TRUE, estintid := cumsum(newemp), by = idpanel]
sipp[lfstat == 1, estintid := na.locf(estintid, na.rm = FALSE), by = idpanel]
sipp[, newemp := NULL]

# create an unemployment stint id
sipp[, newunemp := lfstat == 2 & (last.lfstat == 1 | mis == 1)]
sipp[newunemp == TRUE, ustintid := cumsum(newunemp), by = idpanel]
sipp[lfstat != 1, ustintid := na.locf(ustintid, na.rm = FALSE), by = idpanel]
sipp[, newunemp := NULL]

# create unemployment duration variable
sipp[, unempdur := seq_len(.N) - 1, by = list(idpanel, ustintid)]
sipp[is.na(ustintid), unempdur := NA]

# fill in occupation with next occupation during unemployment stints
sipp[, next.occ := shift(occ, 1, type = "lead"), by = idpanel]
sipp[lfstat != 1, occ := Mode(next.occ), by = list(idpanel, ustintid)]
sipp[, last.occ := shift(occ, 1, type = "lag"), by = idpanel]
sipp[lfstat != 1, last.occ := Mode(last.occ), by = list(idpanel, ustintid)]
sipp[, next2.occ := shift(occ, 2, type = "lead"), by = idpanel]

########## transitions

# create switched occupation dummy
#
# universe:
#
# job != next.job: job code changes for those who are continuing employment
#
# lfstat == 1 & next.lfstat == 2: EU transition. We pulled occupation forwards, so 
# this is when the occupation change should happen
#
sipp[job != next.job | (lfstat == 1 & next.lfstat != 2), 
     switchedOcc := (occ != next.occ) & (last.occ != next.occ) & (occ != next2.occ)]

# create EU/UE/EE dummies
sipp[lfstat == 1, EU := (next.lfstat == 2)]
sipp[lfstat == 2, UE := (next.lfstat == 1)]
sipp[lfstat == 1 & next.lfstat == 1, EE := (job != next.job)]
sipp[lfstat == 1 & next.lfstat != 1, EE := FALSE]

# plot transitions time series for sanity check
EU <- sipp[, .(EU = weighted.mean(EU, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EU, aes(date, EU, color = panel, group = panel)) +
	geom_point() +
	geom_line()

UE <- sipp[, .(UE = weighted.mean(UE, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(UE, aes(date, UE, color = panel, group = panel)) +
	geom_point() +
	geom_line()

EE <- sipp[, .(EE = weighted.mean(EE, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EE, aes(date, EE, color = panel, group = panel)) +
	geom_point() +
	geom_line()

switchedOcc <- sipp[, .(switchedOcc = weighted.mean(switchedOcc, wpfinwgt, na.rm = TRUE)),
		    by = list(panel, date)]
ggplot(switchedOcc, aes(date, switchedOcc, color = panel, group = panel)) +
	geom_point() +
	geom_line()

rm(list=c("EU", "UE", "EE", "switchedOcc"))

########## inflation, unemployment, and recession

# create vector of recession dates
recDates <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01"))

# create recession indicator
sipp[, recIndic := (date > recDates[1] & date < recDates[2]) | 
     	(date > recDates[3] & date < recDates[4])]

# get PCE and unemployment rate data
setwd(datadir)
PCE <- readRDS("./PCE.RData")
unrate <- readRDS("./unrate.RData")
setwd(rootdir)

# add PCE and unemployment rate data
sipp <- merge(sipp, PCE, by = "date", all.x = TRUE)
sipp <- merge(sipp, unrate, by = "date", all.x = TRUE)
rm(list = c("PCE", "unrate"))

# make sure order is preserved after merger
setkey(sipp, idpanel, date)

########## earnings

# cerate lag/lead earnings
sipp[, last.earnm := shift(earnm, 1, type = "lag"), by = idpanel]
sipp[, next.earnm := shift(earnm, 1, type = "lead"), by = idpanel]

# remove large jumps in earnings not explained by UE or EU
sipp[, nomearnm := earnm]
sipp[, earnm := earnm/PCEPI*100]
sipp[, badearn := abs(log(next.earnm/earnm)) > 2.0 & 
      	abs(log(next.earnm/last.earnm)) < 0.1]
sipp[UE == TRUE | EU == TRUE, badearn := FALSE]
sipp[badearn == TRUE, earnm := NA]
sipp[, c("badearn", "nomearnm") := NULL]

########## demographic indicators

sipp[, Young := age < 30]
sipp[, HSCol := (educ >= 4) + (educ >= 2)]

########## missing occs

# missing occupation codes cause issues with the occwage regression in step 4.
# remove them here or there
sipp <- sipp[!is.na(occ),]

########## save prepared data
setwd(outputdir)
saveRDS(sipp, "./preparedSipp.RData")
saveRDS(sipp, "./DTall_3.RData")

