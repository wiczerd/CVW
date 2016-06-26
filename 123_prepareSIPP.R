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
		"eyear","emonth","syear","smonth",
		"ersend",
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

################################## Prepare data  ----------------------------

# load intermediate result if starting from here
if(!exists("sipp")) {
	setwd(outputdir)
	sipp <- readRDS("./sipp.RData")
}

########## date, unique id, months in sample ----------------------------

# add date
sipp <- sipp[!is.na(month) & !is.na(year),]
sipp[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]

# create unique ID across all panels
sipp[, id := paste(id, panel, sep = "")]

# sort by id and date
setkey(sipp, id, date)

# create month in sample variable
sipp[, mis := seq_len(.N), by = id]

########## labor force status ----------------------------

# recode esr
sipp[esr <= 3, lfstat := 1]
sipp[esr >= 4 & esr <= 7, lfstat := 2]
sipp[esr == 8, lfstat := 3]

# create lag/lead lfstat
sipp[, last.lfstat := shift(lfstat, 1, type = "lag"), by = id]
sipp[, next.lfstat := shift(lfstat, 1, type = "lead"), by = id]

########## occupation and job ----------------------------

# replace occ with soc2d (occ will now refer to soc2d)
sipp[, occ := soc2d]
sipp[, soc2d := NULL]

# make sure job and occ are NA if not employed
sipp[lfstat != 1, job := NA]
sipp[lfstat != 1, occ := NA]

# replace occ with mode occ over a job spell
sipp[lfstat == 1, occ := Mode(occ), by = list(id, job)]

# create lag/lead job
sipp[, last.job := shift(job, 1, type = "lag"), by = id]
sipp[, next.job := shift(job, 1, type = "lead"), by = id]

########## stint ids ----------------------------

# create an employment stint id
sipp[, newemp := lfstat == 1 & (last.lfstat >= 2 | mis == 1 | job != last.job)]
sipp[lfstat == 1 & is.na(newemp), newemp := FALSE]
sipp[newemp == TRUE, estintid := cumsum(newemp), by = id]
sipp[lfstat == 1, estintid := na.locf(estintid, na.rm = FALSE), by = id]
sipp[, newemp := NULL]

# create an unemployment stint id
sipp[, newunemp := lfstat == 2 & (last.lfstat == 1 | mis == 1)]
sipp[newunemp == TRUE, ustintid := cumsum(newunemp), by = id]
sipp[lfstat != 1, ustintid := na.locf(ustintid, na.rm = FALSE), by = id]
sipp[, newunemp := NULL]

# create unemployment duration variable
sipp[, unempdur := seq_len(.N) - 1, by = list(id, ustintid)]
sipp[is.na(ustintid), unempdur := NA]

# fill in occupation with next occupation during unemployment stints
sipp[, next.occ := shift(occ, 1, type = "lead"), by = id]
sipp[lfstat != 1, occ := Mode(next.occ), by = list(id, ustintid)]
sipp[, last.occ := shift(occ, 1, type = "lag"), by = id]
sipp[lfstat != 1, last.occ := Mode(last.occ), by = list(id, ustintid)]
sipp[, next2.occ := shift(occ, 2, type = "lead"), by = id]

########## transitions ----------------------------


#use start-date & end-date for EE switches
sipp[lfstat == 1 & next.lfstat == 1, Eend:= (year==eyear & month==emonth)]
sipp[is.na(Eend)==T, Eend:= F]

sipp[lfstat == 1 & next.lfstat == 1, Estart:= (year==syear & month==smonth)]
sipp[is.na(Estart)==T, Estart:= F]
sipp[ , next.Estart := shift(Estart, type="lead"), by=id]

# create EU/UE/EE dummies
sipp[lfstat == 1, EU := (next.lfstat == 2)]
sipp[lfstat == 2, UE := (next.lfstat == 1)]
sipp[lfstat == 1 & next.lfstat==1, EE := (Eend==T | next.Estart==T)]

# these are not stayers, but also not switchers either... some miscoding somewhere
sipp[lfstat == 1 & next.lfstat == 1, jobchng := (job != next.job)]
sipp[lfstat == 1 & next.lfstat != 1, jobchng := NA]


# create switched occupation dummy
#
sipp[ , occ_wave := Mode(occ), by=list(id,wave)]
sipp[ srefmon==4, switchedOcc := occ_wave != shift(occ_wave,type="lead")]
sipp[ , switchedOcc := any(switchedOcc,na.rm=T), by=list(wave,id)]
#sipp[, switchedOcc := (occ != next.occ) & (last.occ != next.occ) & (occ != next2.occ)]


# plot transitions time series for sanity check
EU <- sipp[, .(EU = weighted.mean(EU, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EU, aes(date, EU, color = panel, group = panel)) +
	geom_point() +
	geom_line()

UE <- sipp[, .(UE = weighted.mean(UE, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(UE, aes(date, UE, color = panel, group = panel)) +
	geom_point() +
	geom_line()

EE <- sipp[lfstat==1 & next.lfstat==1 & (is.finite(eyear) | is.finite(syear)), .(EE = weighted.mean( EE==T, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EE, aes(date, EE, color = panel, group = panel)) +
	geom_point() +
	geom_line()

switchedOcc <- sipp[EE==T|EU==T|UE==T, .(switchedOcc = weighted.mean(switchedOcc, wpfinwgt, na.rm = TRUE)),
		    by = list(panel, date)]
ggplot(switchedOcc, aes(date, switchedOcc, color = panel, group = panel)) +
	geom_point() +
	geom_line()

Umo <- sipp[, .(Umo = weighted.mean(lfstat==2, wpfinwgt, na.rm = TRUE)),
					by = list(panel, date)]
ggplot(Umo, aes(date, Umo, color = panel, group = panel)) +
	geom_point() +
	geom_line()


#rm(list=c("EU", "UE", "EE", "switchedOcc"))

########## inflation, unemployment, and recession ----------------------------

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
setkey(sipp, id, date)

########## earnings ----------------------------

# cerate lag/lead earnings
sipp[, last.earnm := shift(earnm, 1, type = "lag"), by = id]
sipp[, next.earnm := shift(earnm, 1, type = "lead"), by = id]

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
sipp <- sipp[!is.na(occ)|lfstat>1,]

# save intermediate result
setwd(outputdir)
saveRDS(sipp, "./sipp_2.RData")


########## compute seam-to-seam status change variable ----------------------------

# load intermediate result if starting from here
if(!exists("sipp")) {
	setwd(outputdir)
	sipp <- readRDS("./sipp_2.RData")
}


setkey(sipp, id,date, physical = T)
# sorting check on wave:
sipp[ , wave_sorted:= is.unsorted(wave), by=id]
if( sipp[ , mean(wave_sorted,na.rm=T)]!=0 ){
	warning("not sorted on wave")
}
sipp[ , wave_sorted:=NULL]

sipp[ , maxmis:= max(mis), by=id]
sipp[ , seam:= (wave != shift(wave,1,type="lead") | mis==maxmis) , by=id ]
sipp[ , maxmis:=NULL]
sipp[, recIndic_wave := any(recIndic, na.rm=T), by=list(wave,id)]

sipp[ , lfstat_wave := max(lfstat,na.rm=T), by=list(id,wave)]
#sipp[ , lfstat2_wave := any(lfstat==2, na.rm=T), by=list(id,wave)] #anytime an unemployment shows up, we count it as such
#sipp[ , lfstat_wave := ifelse(lfstat2_wave, 2L, lfstat_wave)] #anytime an unemployment shows up, we count it as such
#sipp[ , lfstat2_wave:=NULL]

sipp[ , Eend_wave:= any(Eend, na.rm=T), by=list(wave,id)]
sipp[ , Estart_wave:= any(next.Estart, na.rm=T), by=list(wave,id)]


sipp[ seam ==T & is.finite(lfstat), lfstat_seam := lfstat]
sipp[ , lfstat_seam := max(lfstat_seam,na.rm=T), by=list(id,wave)]

#create leads and lags using subsetted dataset
sipp_wave <- subset(sipp, seam==T)
setkey(sipp_wave, id,wave)
sipp_wave[ , job_wave := job]
sipp_wave[ lfstat_wave==2 | lfstat_wave ==3, job_wave := NA]
sipp_wave[ , next.lfstat_wave := shift(lfstat_wave,type="lead"), by=id]
sipp_wave[ , last.lfstat_wave := shift(lfstat_wave,type="lag") , by=id]
sipp_wave[ , next.job_wave    := shift(job_wave,type="lead")   , by=id]
sipp_wave[ , last.job_wave    := shift(job_wave,type="lag")    , by=id]
sipp_wave[ , next.occ_wave    := shift(occ,type="lead")        , by=id]
sipp_wave[ , last.occ_wave    := shift(occ,type="lag")         , by=id]
# create EU/UE/EE dummies
sipp_wave[lfstat_wave == 1, EU_wave := (next.lfstat_wave == 2)]
sipp_wave[lfstat_wave == 2, UE_wave := (next.lfstat_wave == 1)]
sipp_wave[lfstat_wave == 1 & next.lfstat_wave == 1, EE_wave := (job_wave != next.job_wave) & (last.job_wave != next.job_wave)]
sipp_wave[lfstat_wave == 1 & next.lfstat_wave != 1, EE_wave := NA]
sipp_wave[(EU_wave ==T | UE_wave==T | EE_wave==T) , ]

sipp_wave[job_wave != next.job_wave | (lfstat_wave == 1 & next.lfstat_wave == 2), 
	 switchedOcc_wave := (occ != next.occ_wave)]

sipp_wave <- subset(sipp_wave, select=c("next.lfstat_wave","last.lfstat_wave","next.job_wave","last.job_wave","job_wave","next.occ_wave","last.occ_wave",
										"EE_wave","EU_wave","UE_wave","switchedOcc_wave","wave","id"))

sipp <- merge(sipp,sipp_wave, by=c("id","wave"), all=T)

# plot transitions time series for sanity check
EU_wave <- sipp[, .(EU_wave = weighted.mean(EU_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EU_wave, aes(date, EU_wave, color = panel, group = panel)) +
	geom_point() +ylim(c(0,.1))+
	geom_line()

UE_wave <- sipp[, .(UE_wave = weighted.mean(UE_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(UE_wave, aes(date, UE_wave, color = panel, group = panel)) +
	geom_point() +
	geom_line()

EE_wave <- sipp[, .(EE_wave = weighted.mean(EE_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(EE_wave, aes(date, EE_wave, color = panel, group = panel)) +
	geom_point() + ylim(c(0,.1))+
	geom_line()

U_wave <- sipp[ , .(U_wave =  weighted.mean(lfstat_wave==2, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(U_wave, aes(date, U_wave, color = panel, group = panel)) +
	geom_point() + ylim(c(0,.10))+
	geom_line()
U_mo <- sipp[ , .(U_mo =  weighted.mean(lfstat==2, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
ggplot(U_mo, aes(date, U_mo, color = panel, group = panel)) +
	geom_point() + ylim(c(0,.15))+
	geom_line()



########## save prepared data
setwd(outputdir)
saveRDS(sipp, "./preparedSipp.RData")
saveRDS(sipp, "./DTall_3.RData")

