# June 22, 2016
# Daniel Eubanks
# Read in SIPP data and prepare it for analysis
library(foreign)
library(data.table)
library(zoo)
library(lubridate)
library(ggplot2)
library(mFilter)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

Max_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		max(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA_real_
		}
	}
}
Min_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		min(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA_real_
		}
	}
}
Any_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		any(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA
		}
	}
}

#rootdir <- "G:/Research_Analyst/Eubanks/Occupation Switching/SIPP to go"
rootdir <- "~/workspace/CVW/R/"
datadir <- paste(rootdir, "/InputDataDE", sep = "")
outputdir <- paste(rootdir, "/Results", sep = "")
figuredir <- paste(rootdir, "/Figures", sep = "")
intermed_plots = F
final_plots = F
max_wavefreq = 1 # controls whether take max over months in wave or wave-frequency observation
setwd(rootdir)

############################## Read data, combine panels ##############################

# create list of variables to extract from the raw data
toKeep <- c("year", 
	    "month",
	    # demographics
	    "age",
	    "educ",
	    "female",
	    "race","state",
	    # technical and weight variables
	    "id",
	    "wave",
	    "srefmon",
	    "wpfinwgt",
		"epppnum",
	    # labor force status variables
	    "esr","rwkesr2",
	    # job variables
	    "job",
		"eyear","emonth","syear","smonth",
		"ersend","estlemp",
	    "occ","ind","ajbocc",
		# income variables
	    "earnm","earn_imp","ui_a","ui_r","thearn","thtotinc"
		)

########## read in individual panels, extract variables, and subset sample

setwd(datadir)

# 1996 panel
sipp96 <- read.dta("./sippsets96ABDFG.dta", convert.factors = FALSE)
sipp96 <- sipp96[toKeep]
sipp96 <- subset(sipp96, !is.na(esr) & (age >= 16) & (age <= 65))
sipp96 <- data.table(sipp96)
sipp96$panel <- 1996

# 2001 panel
sipp01 <- read.dta("./sippsets01ABDFG.dta", convert.factors = FALSE)
sipp01 <- sipp01[toKeep]
sipp01 <- subset(sipp01, !is.na(esr) & (age >= 16) & (age <= 65))
sipp01 <- data.table(sipp01)
sipp01$panel <- 2001

# 2004 panel
sipp04 <- read.dta("./sippsets04ABDFG.dta", convert.factors = FALSE)
sipp04 <- sipp04[toKeep]
sipp04 <- subset(sipp04, !is.na(esr) & (age >= 16) & (age <= 65) )
sipp04 <- data.table(sipp04)
sipp04$panel <- 2004

# 2008 panel
sipp08 <- read.dta("./sippsets08ABDFG.dta", convert.factors = FALSE)
sipp08 <- sipp08[toKeep]
sipp08 <- subset(sipp08, !is.na(esr) & (age >= 16) & (age <= 65))
sipp08 <- data.table(sipp08)
sipp08$panel <- 2008

rm(toKeep)


########## add soc2d codes using crosswalks

setwd(datadir)
occ1990_soc2d <- readRDS("./occ90_soc2d.RData")
coc2000_occ1990 <- readRDS("./coc2000_occ1990.RData")
coc1990_occ1990 <- readRDS("./coc1990_occ1990.RData")

# 1996 panel
sipp96[, coc:= occ]
sipp96[, coc90 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp96 <- merge(sipp96, coc1990_occ1990, by  = "coc90", all.x = TRUE)
sipp96 <- merge(sipp96, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp96[, c("occ", "coc90") := NULL]
setnames(sipp96, "occ1990", "occ")

# 2001 panel
sipp01[, coc:= occ]
sipp01[, coc90 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp01 <- merge(sipp01, coc1990_occ1990, by  = "coc90", all.x = TRUE)
sipp01 <- merge(sipp01, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp01[, c("occ", "coc90") := NULL]
setnames(sipp01, "occ1990", "occ")

# 2004 panel
sipp04[, coc:= occ]
sipp04[, coc2000 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
sipp04 <- merge(sipp04, coc2000_occ1990, by  = "coc2000", all.x = TRUE)
sipp04 <- merge(sipp04, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp04[, c("occ", "coc2000") := NULL]
setnames(sipp04, "occ1990", "occ")

# 2008 panel
sipp08[, coc:=occ]
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
sipp[, id := as.integer(factor(id)) ]

# sort by id and date
setkey(sipp, id, date)

# create month in sample variable
sipp[, mis := seq_len(.N), by = id]

########## labor force status ----------------------------

#ESR 1 -- With job entire month, worked all weeks.
#ESR 2 -- With job entire month, missed one or more weeks but not because of a layoff.
#ESR 3 -- With job entire month, missed 1 or more weeks because of layoff.
#ESR 4 -- With job part of month, but not because of layoff or looking for work.
#ESR 5 -- With job part of month, some time spent on layoff or looking for work.
#ESR 6 -- No job in month, spent entire month on layoff or looking for work.
#ESR 7 -- No job in month, spent part of month on layoff or looking for work.
#ESR 8 -- No job in month, no time spent on layoff or looking for work

# recode esr
sipp[esr <= 5, lfstat := 1]
sipp[esr >= 6 & esr <= 7, lfstat := 2]
#alt:
#sipp[esr <= 3, lfstat := 1]
#sipp[esr >= 4 & esr <= 7, lfstat := 2]
sipp[esr == 8, lfstat := 3]
# CPS-style
#sipp[rwkesr2 >= 1 & rwkesr2 <= 3, lfstat := as.integer(1)] 
#sipp[rwkesr2 == 4, lfstat := as.integer(2)] 
#sipp[rwkesr2 == 5, lfstat := as.integer(3)]
#sipp[, rwkesr2:= NULL]


#clean the 2->3 and 3->2 into 2
#sipp[ lfstat==3 & (shift(lfstat, 1, type = "lag")==2 | shift(lfstat, 1, type = "lead")==2), lfstat:=2, by=id]

# create lag/lead lfstat
sipp[, last.lfstat := shift(lfstat, 1, type = "lag"), by = id]
sipp[, next.lfstat := shift(lfstat, 1, type = "lead"), by = id]

########## occupation and job ----------------------------

# replace occ with soc2d (occ will now refer to soc2d)
sipp[, occ90 := occ] 
sipp[, occ   := soc2d]
sipp[ajbocc>0, occ := NA ]
sipp[ajbocc>0, coc := NA ]

# very frequently occ is missing during an employment spell. Fill it forward then back where there're no conflicts:
sipp[ lfstat==1 & is.finite(job), occ:= na.locf(occ,na.rm = F), by=list(id,job)]
sipp[ lfstat==1 & is.finite(job), occ:= na.locf(occ,na.rm = F,fromLast = T), by=list(id,job)]
sipp[ lfstat==1 & is.finite(job), ind:= na.locf(ind,na.rm = F), by=list(id,job)]
sipp[ lfstat==1 & is.finite(job), ind:= na.locf(ind,na.rm = F,fromLast = T), by=list(id,job)]

# make sure job and occ are NA if not employed
sipp[lfstat != 1, job := NA]
sipp[lfstat != 1, occ := NA]
sipp[lfstat != 1, ind := NA]
#save the variable before it's been filled over lfstat>1 gaps
sipp[ , occ_nofill := occ]
sipp[ , ind_nofill := occ]

########## stint ids ----------------------------

# create an employment stint id
sipp[, last.job := shift(job, 1, type = "lag"), by = id]
sipp[, newemp := lfstat == 1 & (last.lfstat >= 2 | mis == 1 | job != last.job)]
sipp[lfstat == 1 & is.na(newemp), newemp := FALSE]
sipp[newemp == TRUE, estintid := cumsum(newemp), by = id]
sipp[lfstat == 1, estintid := na.locf(estintid, na.rm = FALSE), by = id]
sipp[, newemp := NULL]

# create an unemployment stint id
sipp[, newunemp := lfstat >= 2 & (last.lfstat == 1 | mis == 1)]
sipp[newunemp == TRUE, ustintid := seq_len(.N), by = id]
sipp[lfstat == 1, ustintid := 9999]
sipp[, ustintid := na.locf(ustintid, na.rm = FALSE), by = id]
sipp[lfstat == 1 | ustintid  == 9999, ustintid := NA]
sipp[, newunemp := NULL]

# create unemployment duration variable
sipp[lfstat==2, unempdur := seq_len(.N) - 1, by = list(id, ustintid)]  #count U entries in each ustintid
sipp[is.na(ustintid), unempdur := NA]
# create max unemployment duration
sipp[lfstat==2, max.unempdur := max(unempdur,na.rm = T), by = list(id, ustintid)]

# This removed occupation codes that changed during a job spell:
#sipp[ , varoccjob := var(occ,na.rm = T), by=list(id,estintid)] #in about 6\% of cases, occupation changes during job
#sipp[ , varjobidjob := var(job,na.rm = T), by=list(id,estintid)]
#sipp[ lfstat==1, occ_mode := Mode(occ), by=list(id,estintid)]
#sipp[ varoccjob>0 & varjobidjob==0 & lfstat==1, occ:=occ_mode, by=list(id,estintid)]
# some employment not list occupation
#sipp[ lfstat==1 & is.na(occ), occ:= occ_mode]
#sipp[ , c("varoccjob","varjobidjob","occ_mode"):=NULL]

# fill in occupation with next occupation during unemployment stints
sipp[, next.occ := shift(occ, 1, type = "lead"), by = id]
sipp[lfstat != 1, occ := Mode(next.occ), by = list(id, ustintid)]
sipp[, next.ind := shift(ind, 1, type = "lead"), by = id]
sipp[lfstat != 1, ind := Mode(next.ind), by = list(id, ustintid)]
sipp[, last.occ := shift(occ, 1, type = "lag"), by = id]
sipp[lfstat != 1, last.occ := Mode(last.occ), by = list(id, ustintid)]

sipp[, next.occ := shift(occ, 1, type = "lead"), by = id]
sipp[, next.ind := shift(ind, 1, type = "lead"), by = id]
sipp[, switchedOcc := occ !=next.occ]
sipp[, switchedInd := ind !=next.ind]


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

sipp[is.finite(lfstat) & is.finite(next.lfstat) & is.na(EU), EU :=F ]
sipp[is.finite(lfstat) & is.finite(next.lfstat) & is.na(UE), UE :=F ]
sipp[is.finite(lfstat) & is.finite(next.lfstat) & is.na(EE), EE :=F ]

# create job change id.  Job id generally measured only at the wave frequency, mostly this shows up in srefmon=4
# sipp[srefmon==4, next4.job := (job != next.job)]
# sipp[lfstat == 1 & next.lfstat != 1, jobchng := NA]
sipp[ , next.job := shift(job,type="lead"),by=id]
sipp[lfstat==1 & next.lfstat==1 , jobchng := (job!=next.job)]

#add EU to ustintid
sipp[ , next.ustintid := shift(ustintid,type="lead"), by=id]
sipp[ EU==T, ustintid:= next.ustintid]
sipp[ , next.ustintid:=NULL]

#match UE and EU:
sipp[UE==T, mis_UE:=mis]
sipp[EU==T, mis_EU:=mis]
sipp[ ustintid>0, mis_UE := Max_narm(mis_UE), by =list(id,ustintid)]
sipp[ ustintid>0, mis_EU := Min_narm(mis_EU), by =list(id,ustintid)]
sipp[ ustintid>0, matched_EUUE := (mis_EU<=mis_UE) ]
sipp[ ustintid>0 & is.na(matched_EUUE)==T, matched_EUUE:=F ]
sipp[, c("mis_UE","mis_EU"):=NULL]


# plot transitions time series for sanity check
if(intermed_plots==T){
	EU <- sipp[lfstat==1, .(EU = weighted.mean(EU, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(EU, aes(date, EU, color = panel, group = panel)) +
		geom_point() + theme_bw()+
		geom_line() + xlab("") + ylab("EU monthly rate")
	ggsave(filename=paste0(outputdir,"/EUmo.eps"),height= 5,width=10)
	
	UE <- sipp[lfstat==2, .(UE = weighted.mean(UE, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(UE, aes(date, UE, color = panel, group = panel)) +
		geom_point() +
		geom_line()+ xlab("") + ylab("UE monthly rate")+theme_bw()
	ggsave(filename=paste0(outputdir,"/UEmo.eps"),height= 5,width=10)
	
	EE <- sipp[lfstat==1 & next.lfstat==1 & (is.finite(eyear) | is.finite(syear)), .(EE = weighted.mean( EE==T, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(EE, aes(date, EE, color = panel, group = panel)) +
		geom_point() +
		geom_line()+xlab("")+ ylab("EE monthly rate")+theme_bw()
	ggsave(filename=paste0(outputdir,"/EEmo.eps"),height= 5,width=10)
	
	Umo <- sipp[lfstat<3, .(Umo = weighted.mean(lfstat==2, wpfinwgt, na.rm = TRUE)),
						by = list(panel, date)]
	ggplot(Umo, aes(date, Umo, color = panel, group = panel)) +
		geom_point() +
		geom_line()
	
	EUsw <- sipp[lfstat==1 & EU==T, .(EU = weighted.mean(switchedOcc, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(EU, aes(date, EU, color = panel, group = panel)) +
		geom_point() + theme_bw()+
		geom_line() + xlab("") + ylab("EU monthly rate")
	ggsave(filename=paste0(outputdir,"/EUmo.eps"),height= 5,width=10)
	
	rm(list=c("EU", "UE", "EE"))
}

########## inflation, unemployment, and recession ----------------------------

# create vector of recession dates
recDates <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01"))
# create vector of recession dates : above 6%
recDates2 <- as.Date(c("2003-04-01", "2003-10-01","2008-08-01", "2014-09-01"))

# get PCE and unemployment rate data
setwd(datadir)
PCE <- readRDS("./PCE.RData")
unrate <- readRDS("./unrate.RData")
unrate <- readRDS("./CPSunempRt.RData")
setwd(rootdir)

filtered.unrate <- hpfilter(unrate$unrt,type = "lambda",freq=14400)
unrate<-cbind(unrate,filtered.unrate$cycle)
# create recession indicator
sipp[, recIndic := (date > recDates[1] & date < recDates[2]) | 
	 	(date > recDates[3] & date < recDates[4])]
#recDates2 <- as.Date(unrate$date[unrate$unrateSA>6.5])
recDates2 <- unrate$date[ unrate$`filtered.unrate$cycle`>=quantile(unrate$`filtered.unrate$cycle`,probs = 4/5) ]
sipp[, recIndic2 := date %in% recDates2]


# add PCE and unemployment rate data
sipp <- merge(sipp, PCE, by = "date", all.x = TRUE)
sipp <- merge(sipp, unrate, by = "date", all.x = TRUE)
rm(list = c("PCE", "unrate"))

# make sure order is preserved after merger
setkey(sipp, id, date)

########## earnings ----------------------------

# create lag/lead earnings
sipp[, last.earnm := shift(earnm, 1, type = "lag"), by = id]
sipp[, next.earnm := shift(earnm, 1, type = "lead"), by = id]

# remove large jumps in earnings not explained by UE or EU
sipp[, nomearnm := earnm]
sipp[, earnm := earnm/PCEPI*100]
sipp[, badearn := abs(log(next.earnm/earnm)) > 2.0 & 
      	abs(log(next.earnm/last.earnm)) < 0.1]
sipp[UE == T | EU == T | EE==T, badearn := F]
sipp[badearn == T, earnm := NA]
sipp[, c("badearn", "nomearnm") := NULL]

########## demographic indicators/ selection

sipp <- subset(sipp, age >=16 & age<=65)

sipp[, Young := age <= 30]
sipp[age < 31, ageGrp := 1]
sipp[age >= 31 & age < 56, ageGrp := 2]
sipp[age >= 56, ageGrp := 3]
sipp[, HSCol := (educ >= 4) + (educ >= 2)]



########## attrition weights
sipp[ , maxmis:= max(mis), by=id]
# sipp[, panelmaxmis := max(maxmis), by=panel]
# sipp[, attrition := maxmis < 0.9*panelmaxmis]
# # use LPM to model probability of dropping out
# attritionModel <- lm(attrition ~ age + female + earnm + factor(educ) + factor(soc2d),
# 					 data=sipp, weights=1/maxmis, na.action=na.exclude)
#sipp[, probAttrition := predict(attritionModel)]
#sipp[, attritionWeights := wpfinwgt/(1-probAttrition)]
#sipp[is.na(attritionWeights), attritionWeights := wpfinwgt]


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
sipp[ , seam:= (wave != shift(wave,1,type="lead") | mis==maxmis) , by=id ]
sipp[ , seammon := max(srefmon,na.rm = T), by=list(id,wave)]

sipp[, recIndic_wave := any(recIndic, na.rm=T), by=list(wave,id)]
sipp[, recIndic2_wave := any(recIndic2, na.rm=T), by=list(wave,id)]

sipp[ , lfstat_wave := as.integer(max(lfstat,na.rm=T)), by=list(id,wave)]
sipp[ lfstat_wave>1, lfstat_wave := 3L-any(lfstat==2), by=list(id,wave)]

sipp[lfstat_wave==2 , max.unempdur_wave := max(max.unempdur, na.rm = T) , by=list(id,wave)]

sipp[ , Eend_wave:= any(Eend, na.rm=T), by=list(wave,id)]
sipp[ , Estart_wave:= any(Estart, na.rm=T), by=list(wave,id)]

sipp[ EE==T, EEmon := srefmon]
sipp[ is.na(EEmon), EEmon := 0]
sipp[ , EEmon := max(EEmon), by=list(id,wave)]
sipp[ EU==T, EUmon := srefmon]
sipp[ is.na(EUmon), EUmon := 0]
sipp[ , EUmon := max(EUmon), by=list(id,wave)]
sipp[ UE==T, UEmon := srefmon]
sipp[ is.na(UEmon), UEmon := 0]
sipp[ , UEmon := max(UEmon), by=list(id,wave)]
#now check EU, UE, EE as max over months in the wave
sipp[ , UE_max := any(UE,na.rm=T), by=list(id,wave)]
sipp[ , EU_max := any(EU,na.rm=T), by=list(id,wave)]
sipp[ , EE_max := any(EE,na.rm=T), by=list(id,wave)]
sipp[ , switchedOcc_max := Any_narm(switchedOcc), by=list(id,wave)]
sipp[ , matched_EUUE_max := any(matched_EUUE,na.rm = T), by=list(id,wave)]

sipp[ EU==T, occ_wave := occ]
sipp[ EU==T, ind_wave := ind]
sipp[ UE==T, occ_wave := occ]
sipp[ UE==T, ind_wave := ind]
sipp[  (UE_max==T|EU_max==T), occ_wave := Mode(occ_wave), by=list(id,wave)]
sipp[  (UE_max==T|EU_max==T), ind_wave := Mode(ind_wave), by=list(id,wave)]
sipp[ !(UE_max==T|EU_max==T) & seam==T, occ_wave := occ]
sipp[ !(UE_max==T|EU_max==T) & seam==T, ind_wave := ind]


sipp[ , jobchng_max := any(jobchng,na.rm=T), by=list(id,wave)]

#**************************************************************
#sipp_wave begins here----------------------

sipp_wave <- subset(sipp, seam==T)
setkeyv(sipp_wave, c("id","wave"))
sipp_wave[ , job_wave := job]

#create leads and lags using subsetted dataset
sipp_wave[                    , next.wave := shift(wave,type = "lead"),by=id]
sipp_wave[                    , last.wave := shift(wave,type = "lag" ),by=id]
sipp_wave[ next.wave == wave+1, next.lfstat_wave := shift(lfstat_wave,type="lead"), by=id]
sipp_wave[ last.wave == wave-1, last.lfstat_wave := shift(lfstat_wave,type="lag") , by=id]
sipp_wave[ next.wave == wave+1, next.job_wave    := shift(job_wave,type="lead")   , by=id]
sipp_wave[ last.wave == wave-1, last.job_wave    := shift(job_wave,type="lag")    , by=id]
sipp_wave[ next.wave == wave+1, next.occ_wave    := shift(occ_wave,type="lead")   , by=id]
sipp_wave[ last.wave == wave-1, last.occ_wave    := shift(occ_wave,type="lag")    , by=id]
sipp_wave[ next.wave == wave+1, next.ind_wave    := shift(ind_wave,type="lead")   , by=id]
sipp_wave[ last.wave == wave-1, last.ind_wave    := shift(ind_wave,type="lag")    , by=id]

##create wave-level ustintid_wave
#sipp[ , ustintid_wave := Mode(ustintid), by=list(wave,id)]
#test: ~0.8% of waves have mutliple spells sipp[ , varustintwave := var(ustintid,na.rm = T), by=list(wave,id)]
# or compute akin to how we did earlier
sipp_wave[ , wis := seq_len(.N), by = id]
sipp_wave[ , maxwis := Max_narm(wis), by = id]
sipp_wave[, newunemp_wave := lfstat_wave >1 & (last.lfstat_wave == 1 | wis == 1)]
sipp_wave[newunemp_wave == T, ustintid_wave := cumsum(newunemp_wave), by = id]
sipp_wave[lfstat_wave   == 1, ustintid_wave := 9999]
sipp_wave[                  , ustintid_wave := na.locf(ustintid_wave, na.rm = F), by = id]
sipp_wave[lfstat_wave == 1 | ustintid_wave  == 9999, ustintid_wave := NA]
sipp_wave[, newunemp_wave := NULL]
sipp_wave[ , next.ustintid_wave:= shift(ustintid_wave,type="lead")]

# create EU/UE/EE dummies
if(max_wavefreq==2){
	sipp_wave[lfstat_wave == 1, EU_wave := (next.lfstat_wave == 2)]
	sipp_wave[lfstat_wave == 2, UE_wave := (next.lfstat_wave == 1)]
	sipp_wave[lfstat_wave == 1 & next.lfstat_wave == 1 , EE_wave := (Eend_wave == T | Estart_wave==T)]
	sipp_wave[lfstat_wave == 1 & next.lfstat_wave == 1 & EEmon==4, EE_wave := (Eend_wave == T | next.Estart_wave==T)]
}else{ #max_wavefreq=1
	sipp_wave[ , UE_wave := UE_max]
	sipp_wave[ , EU_wave := EU_max]
	sipp_wave[ , EE_wave := EE_max]
}

#set up jobchng_wave and EU_wave with ustintid
sipp_wave[ EE_wave==T & (lfstat_wave>1 | next.lfstat_wave>1), EE_wave:=F]
if(max_wavefreq==1){
	sipp_wave[last.lfstat_wave == 1 & next.lfstat_wave == 1 , jobchng_wave := jobchng_max ] #& (last.job_wave != next.job_wave)
	sipp_wave[ , next.jobchng_wave := shift(jobchng_wave,type="lead"), by=id]
	#sipp_wave[is.na(last.lfstat_wave) == T & is.na(lfstat_wave) == F &
	#		  	next.lfstat_wave == 1 , jobchng_wave := (last.job_wave != next.job_wave) ] #in case it's the first observation
	#This cleaning will happen in 5_ with the rest of the cleaning sipp_wave[EE_wave==T & !(jobchng_wave==T) , EE_wave := NA] 
	#add EU_wave, EUmon==4 to next ustintid_wave
	sipp_wave[EU_wave==T & EUmon==seammon, ustintid_wave:=next.ustintid_wave  ]
}else{#wavefreq==2
	sipp_wave[lfstat_wave == 1 & next.lfstat_wave == 1 , jobchng_wave := (job_wave != next.job_wave) ] #& (last.job_wave != next.job_wave)
	#sipp_wave[EE_wave==T & !(jobchng_wave==T) , EE_wave := NA] #knocks out ~5% of the changes
}
#add max.unempdur_wave to EU_wave
sipp_wave[ (EU_wave==T|lfstat_wave==2) , max.unempdur_wave:= Max_narm(max.unempdur_wave ), by=list(id,ustintid_wave)]

#fill in occupation over u spells and compute switching
#sipp_wave[ !(EU_wave ==T |UE_wave==T), occtest:=na.locf(occtest,na.rm=F,fromLast=T),by=id]
sipp_wave[(EU_wave==T|lfstat_wave==1|UE_wave==T) & is.na(occ_wave) ,occ_wave:=999] #don't fill in missing EU's 
sipp_wave[(EU_wave==T|lfstat_wave==1|UE_wave==T) & is.na(ind_wave) ,ind_wave:=999] #don't fill in missing EU's 
sipp_wave[ ,occ_wave:=na.locf(occ_wave,na.rm = F,fromLast = T), by=id]
sipp_wave[ ,ind_wave:=na.locf(ind_wave,na.rm = F,fromLast = T), by=id]
sipp_wave[(EU_wave==T|lfstat_wave==1|UE_wave==T) & occ_wave>=990 ,occ_wave:=NA] #don't fill in missing EU's 
sipp_wave[(EU_wave==T|lfstat_wave==1|UE_wave==T) & ind_wave>=990 ,ind_wave:=NA] #don't fill in missing EU's 

sipp_wave[ wave+1 == shift(wave,type="lead"), next.occ_wave := shift(occ_wave,type="lead") , by=id] #recompute now that I've filled back U-spells
sipp_wave[ wave-1 == shift(wave,type="lag" ), last.occ_wave := shift(occ_wave,type="lag" ) , by=id] #recompute now that I've filled back U-spells
sipp_wave[ EU_wave==T & EUmon<seammon , switchedOcc_wave := last.occ_wave != next.occ_wave]
sipp_wave[ EU_wave==T & EUmon==seammon, switchedOcc_wave := occ_wave != next.occ_wave]
sipp_wave[ EE_wave==T & EEmon<seammon , switchedOcc_wave := last.occ_wave != next.occ_wave]
sipp_wave[ EE_wave==T & EEmon==seammon, switchedOcc_wave := occ_wave != next.occ_wave]
sipp_wave[ lfstat_wave==1 & next.lfstat_wave==1 & !(EE_wave==T|EU_wave==T|UE_wave==T), switchedOcc_wave := occ_wave != next.occ_wave]

sipp_wave[ wave+1 == shift(wave,type="lead"), next.ind_wave := shift(ind_wave,type="lead") , by=id]
sipp_wave[ wave-1 == shift(wave,type="lag" ), last.ind_wave := shift(ind_wave,type="lag" ) , by=id] #recompute now that I've filled back U-spells
sipp_wave[ EU_wave==T & EUmon<seammon , switchedInd_wave := last.ind_wave != next.ind_wave]
sipp_wave[ EU_wave==T & EUmon==seammon, switchedInd_wave := ind_wave != next.ind_wave]
sipp_wave[ EE_wave==T & EEmon<seammon , switchedInd_wave := last.ind_wave != next.ind_wave]
sipp_wave[ EE_wave==T & EEmon==seammon, switchedInd_wave := ind_wave != next.ind_wave]
sipp_wave[ lfstat_wave==1 & next.lfstat_wave==1 & !(EE_wave==T|EU_wave==T|UE_wave==T), switchedInd_wave := ind_wave != next.ind_wave]

sipp_wave[ , next.switchedOcc_wave:= shift(switchedOcc_wave,type="lead"),by=id]
sipp_wave[ , next.switchedInd_wave:= shift(switchedInd_wave,type="lead"),by=id]

#compute matched EUUE by ustintid ---- this will replace the one taken from monthly, level
sipp_wave[ UE_wave==T, wis_UE_wave:= wis]
sipp_wave[ EU_wave==T, wis_EU_wave:= wis]
sipp_wave[ ustintid_wave>0, wis_UE_wave := Max_narm(wis_UE_wave), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0, wis_EU_wave := Min_narm(wis_EU_wave), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0, any_EU_wave := any(EU_wave), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0, any_UE_wave := any(UE_wave), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0, matched_EUUE_wave := wis_EU_wave<=wis_UE_wave]
sipp_wave[ is.na(matched_EUUE_wave), matched_EUUE_wave:= F]


#correct for w/in wave transitions
sipp_wave[ , next.EEmon   := shift( EEmon  , type="lead"), by=id]# adjust because EE_wave will be counted before
sipp_wave[ , next.EE_wave := shift( EE_wave, type="lead"), by=id]
sipp_wave[ , next.recIndic_wave := shift( recIndic_wave, type="lead"), by=id]
sipp_wave[ , next.recIndic2_wave := shift( recIndic2_wave, type="lead"), by=id]
sipp_wave[ next.EEmon>0 & next.EEmon<4 & next.EE_wave==T & lfstat_wave==1 , midEE  :=T  ]
sipp_wave[ is.na(midEE)==T, midEE:=F]
sipp_wave[ midEE  ==T , EEmon := next.EEmon]
sipp_wave[ midEE  ==T , EE_wave:=T  ]
sipp_wave[ midEE  ==T , switchedOcc_wave:=next.switchedOcc_wave  ]
sipp_wave[ midEE  ==T , switchedInd_wave:=next.switchedInd_wave  ]
sipp_wave[ midEE  ==T , jobchng_wave:=next.jobchng_wave  ]
sipp_wave[ midEE  ==T , recIndic_wave := next.recIndic_wave  ]
sipp_wave[ midEE  ==T , recIndic2_wave:= next.recIndic2_wave  ]


sipp_wave[ , next.EUmon   := shift( EUmon  , type="lead" ), by= id]
sipp_wave[ , next.EU_wave := shift( EU_wave, type="lead" ), by= id]
sipp_wave[ , next.max.unempdur_wave := shift( max.unempdur_wave, type="lead" ), by= id]
sipp_wave[ next.EUmon>0 & next.EUmon<seammon & next.EU_wave ==T & lfstat_wave==1 , midEU  :=T  ]
sipp_wave[ is.na(midEU)==T, midEU:=F]
sipp_wave[ midEU ==T , EUmon := next.EUmon]
sipp_wave[ midEU ==T , EU_wave:=T  ]
sipp_wave[ midEU ==T , max.unempdur_wave:=next.max.unempdur_wave  ]
sipp_wave[ midEU  ==T , recIndic_wave := next.recIndic_wave  ]
sipp_wave[ midEU  ==T , recIndic2_wave:= next.recIndic2_wave  ]

sipp_wave[ , next.UEmon   := shift( UEmon  , type="lead" ), by= id]
sipp_wave[ , next.UE_wave := shift( UE_wave, type="lead" ), by= id]
sipp_wave[ next.UEmon>0 & next.UEmon<seammon & next.UE_wave ==T & lfstat_wave>=2 , midUE  :=T  ]
sipp_wave[ is.na(midUE)==T, midUE:=F]
sipp_wave[ midUE == T , UEmon := next.UEmon]
sipp_wave[ midUE  ==T , recIndic_wave := next.recIndic_wave  ]
sipp_wave[ midUE  ==T , recIndic2_wave:= next.recIndic2_wave  ]
sipp_wave[ midUE == T , UE_wave:=T  ]

# clean-up the EU, UE in 1 wave
sipp_wave[ EU_wave==T & UE_wave==T & midUE==T , UE_wave := F]
sipp_wave[ EU_wave==T & UE_wave==T & midEU==T , EU_wave := F]
sipp_wave[ UE_wave==T & EU_wave==T & UEmon>EUmon & shift(EU_wave,type="lag")==T, EU_wave:=F]
sipp_wave[ UE_wave==T & EU_wave==T & UEmon>EUmon & shift(UE_wave,type="lead")==T, UE_wave:=F]

sipp_wave[ UE_wave==T | EU_wave==T, UEfollows:=shift(UE_wave,type="lead")==T , by=id]
sipp_wave[ UE_wave==T | EU_wave==T, EUpreceds:=shift(EU_wave,type="lag")==T , by=id]
sipp_wave[ UE_wave==T & EU_wave==T & EUpreceds==T & UEfollows==T & UEmon<EUmon, mkUEEUf :=T ] #combine unemp spells with a UE EU in the middle of a wave
sipp_wave[ mkUEEUf ==T, EU_wave :=F ]
sipp_wave[ mkUEEUf ==T, UE_wave :=F ]
sipp_wave[ UE_wave==T & EU_wave==T & UEfollows==T, UE_wave :=F ] #combine unemp spells with a UE EU in the middle of a wave
sipp_wave[ UE_wave==T & EU_wave==T & EUpreceds==T, EU_wave :=F ] #combine unemp spells with a UE EU in the middle of a wave

#associate ustintid_wave for the middle ones
sipp_wave[UE_wave==T & midUE ==T , ustintid_wave   :=next.ustintid_wave  ]
sipp_wave[EU_wave==T & midEU ==T , ustintid_wave   :=next.ustintid_wave  ]

#EU,UE trumps EE
sipp_wave[ (EU_wave==T | UE_wave==T) & EE_wave==T, EE_wave:=F]

#add switchedOcc_wave to whole ustintid and create occL, occD
sipp_wave[ ustintid_wave>0 & EU_wave!=T & midEU!=T, switchedOcc_wave:=NA]
sipp_wave[ ustintid_wave>0 , switchedOcc_wave := Any_narm(switchedOcc_wave), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0 & EU_wave!=T & midEU!=T, switchedInd_wave:=NA]
sipp_wave[ ustintid_wave>0 , switchedInd_wave := Any_narm(switchedInd_wave), by=list(id,ustintid_wave)]
sipp_wave[ EU_wave==T & midEU!=T & EUmon< seammon, occL:= last.occ_wave]
sipp_wave[ EU_wave==T & midEU!=T & EUmon==seammon, occL:= occ_wave]
sipp_wave[ EU_wave==T & midEU!=T, occD := next.occ_wave]
sipp_wave[ ustintid_wave>0, occL := Max_narm(occL), by=list(id,ustintid_wave)]
sipp_wave[ ustintid_wave>0, occD := Max_narm(occD), by=list(id,ustintid_wave)]
sipp_wave[ EE_wave ==T & EEmon<seammon, occL:= last.occ_wave]
sipp_wave[ EE_wave ==T & EEmon==seammon, occL:= occ_wave]
sipp_wave[ EE_wave ==T , occD:= next.occ_wave]
sipp_wave[ !(EE_wave ==T|EU_wave==T|UE_wave==T) , occD:= next.occ_wave]
sipp_wave[ !(EE_wave ==T|EU_wave==T|UE_wave==T) , occL:= occ_wave]

#create recIndic_stint
sipp_wave[is.na(ustintid_wave)|ustintid_wave==0 , recIndic_stint := recIndic_wave]
sipp_wave[is.finite(ustintid_wave), recIndic_stint := any(recIndic_wave,na.rm=T), by=list(id,ustintid_wave)]
sipp_wave[is.na(ustintid_wave)|ustintid_wave==0 , recIndic2_stint := recIndic2_wave]
sipp_wave[is.finite(ustintid_wave), recIndic2_stint := any(recIndic2_wave,na.rm=T), by=list(id,ustintid_wave)]
#count by UE or EU
sipp_wave[UE_wave==T & midUE==F, recIndic_UE := recIndic_wave]
sipp_wave[is.finite(ustintid_wave), recIndic_UE := any(recIndic_UE,na.rm=T), by=list(id,ustintid_wave)]
sipp_wave[is.na(ustintid_wave)|ustintid_wave==0 , recIndic_UE := recIndic_wave]
sipp_wave[EU_wave==T & midEU==F, recIndic_EU := recIndic_wave]
sipp_wave[is.finite(ustintid_wave), recIndic_EU := any(recIndic_EU,na.rm=T), by=list(id,ustintid_wave)]
sipp_wave[is.na(ustintid_wave)|ustintid_wave==0 , recIndic_EU := recIndic_wave]



#save intermediate result:
saveRDS(sipp_wave, file=paste0(outputdir,"/sipp_wave.RData"))

sipp_wave <- subset(sipp_wave, select=c("next.lfstat_wave","next.job_wave","job_wave","next.occ_wave","occ_wave","occL","occD","ind_wave",
										"jobchng_wave","EE_wave","EU_wave","UE_wave","matched_EUUE_wave","midEE","midEU","midUE","EEmon","UEmon","EUmon","ustintid_wave",
										"recIndic_stint","recIndic2_stint","recIndic_EU","recIndic_UE","max.unempdur_wave","switchedOcc_wave","switchedInd_wave","wave","id"))

sipp[ , c("EEmon","EUmon","UEmon","max.unempdur_wave"):=NULL]


sipp <- merge(sipp,sipp_wave, by=c("id","wave"), all=T)

#*********************************************************
#compute annual versions of things    ----------------------

sipp[ , c("EE_max","EU_max","UE_max"):=NULL]
sipp[ , date0 := min(date), by=id]
sipp[ date>= date0         & date<(date0+365  ), yri := 1]
sipp[ date>= (date0+365)   & date<(date0+365*2), yri := 2]
sipp[ date>= (date0+365*2) & date<(date0+365*3), yri := 3]
sipp[ date>= (date0+365*3) & date<(date0+365*4), yri := 4]
sipp[ date>= (date0+365*4) & date<(date0+365*5), yri := 5]
sipp[ date>= (date0+365*5) & date<(date0+365*6), yri := 6]

sipp[ , yrseam := (yri !=shift(yri,type="lead")), by=id]
sipp[ , month0:= month(date0)]
sipp[ , yrrefmon := month(date) - month0+1]
sipp[ yrrefmon<=0, yrrefmon := 12+yrrefmon]

#annual transitions
sipp[, recIndic_ann := any(recIndic, na.rm=T), by=list(yri,id)]
sipp[, recIndic2_ann := any(recIndic2, na.rm=T), by=list(yri,id)]

sipp[ , lfstat_ann := as.integer(max(lfstat,na.rm=T)), by=list(id,yri)]
sipp[ lfstat_ann>1, lfstat_ann := 3L-any(lfstat==2), by=list(id,yri)]

sipp[ EE==T, EEyrmon := yrrefmon]
sipp[ is.na(EEyrmon), EEyrmon := 0]
sipp[ , EEyrmon := max(EEyrmon), by=list(id,yri)]
sipp[ EU==T, EUyrmon := yrrefmon]
sipp[ is.na(EUyrmon), EUyrmon := 0]
sipp[ , EUyrmon := max(EUyrmon), by=list(id,yri)]
sipp[ UE==T, UEyrmon := yrrefmon]
sipp[ is.na(UEyrmon), UEyrmon := 0]
sipp[ , UEyrmon := max(UEyrmon), by=list(id,yri)]
#now check EU, UE, EE as max over months in the year
sipp[ , UE_max := any(UE,na.rm=T), by=list(id,yri)]
sipp[ , EU_max := any(EU,na.rm=T), by=list(id,yri)]
sipp[ , EE_max := any(EE,na.rm=T), by=list(id,yri)]

sipp[ , jobchng_ann := any(jobchng,na.rm=T), by=list(id,yri)]

sipp[ , switchedOcc_ann := any(switchedOcc_wave,na.rm = T)]
sipp[ , switchedInd_ann := any(switchedInd_wave,na.rm = T)]

#**************************************************************
#sipp_ann begins here----------------------

sipp_ann <- subset(sipp, yrseam==T)
setkey(sipp_ann, id,yri)

sipp_ann[ , EU_ann := EU_max]
sipp_ann[ , UE_ann := UE_max]
sipp_ann[ , EE_ann := EE_max]

sipp_ann[ , next.jobchng_ann:= shift(jobchng_ann,type="lead"),by=id]

#correct for w/in wave transitions
sipp_ann[ , next.EEyrmon   := shift( EEyrmon  , type="lead"), by=id]# adjust because EE_ann will be counted before
sipp_ann[ , next.EE_ann := shift( EE_ann, type="lead"), by=id]
sipp_ann[ , next.recIndic_ann := shift( recIndic_ann, type="lead"), by=id]
sipp_ann[ , next.recIndic2_ann := shift( recIndic2_ann, type="lead"), by=id]
sipp_ann[ next.EEyrmon>0 & next.EEyrmon<12 & next.EE_ann==T & lfstat_ann==1 , midEE_ann :=T  ]
sipp_ann[ is.na(midEE_ann)==T, midEE_ann:=F]
sipp_ann[ midEE_ann  ==T , EEyrmon := next.EEyrmon]
sipp_ann[ midEE_ann  ==T , EE_ann:=T  ]
sipp_ann[ midEE_ann  ==T , jobchng_ann:=next.jobchng_ann  ]
#sipp_ann[ midEE_ann  ==T , recIndic_ann := next.recIndic_ann  ]
#sipp_ann[ midEE_ann  ==T , recIndic2_ann:= next.recIndic2_ann  ]


sipp_ann[ , next.EUyrmon   := shift( EUyrmon  , type="lead" ), by= id]
sipp_ann[ , next.EU_ann := shift( EU_ann, type="lead" ), by= id]
sipp_ann[ next.EUyrmon>0 & next.EUyrmon<12 & next.EU_ann ==T & lfstat_ann==1 , midEU_ann :=T  ]
sipp_ann[ is.na(midEU_ann)==T, midEU_ann:=F]
sipp_ann[ midEU_ann ==T , EUyrmon := next.EUyrmon]
sipp_ann[ midEU_ann ==T , EU_ann:=T  ]
#sipp_ann[ midEU_ann  ==T , recIndic_ann := next.recIndic_ann  ]
#sipp_ann[ midEU_ann  ==T , recIndic2_ann:= next.recIndic2_ann  ]

sipp_ann[ , next.UEyrmon   := shift( UEyrmon  , type="lead" ), by= id]
sipp_ann[ , next.UE_ann := shift( UE_ann, type="lead" ), by= id]
sipp_ann[ next.UEyrmon>0 & next.UEyrmon<12 & next.UE_ann ==T & lfstat_ann>=2 , midUE_ann :=T  ]
sipp_ann[ is.na(midUE_ann)==T, midUE_ann:=F]
sipp_ann[ midUE_ann == T , UEyrmon := next.UEyrmon]
#sipp_ann[ midUE_ann  ==T , recIndic_ann := next.recIndic_ann  ]
#sipp_ann[ midUE_ann  ==T , recIndic2_ann:= next.recIndic2_ann  ]
sipp_ann[ midUE_ann == T , UE_ann:=T  ]

#EU,UE trumps EE
sipp_ann[ (EU_ann==T | UE_ann==T) & EE_ann==T, EE_ann:=F]

#merge them back
sipp_ann <- subset(sipp_ann, select=c(	"jobchng_ann","EE_ann","EU_ann","UE_ann","midEE_ann","midEU_ann","midUE_ann","EEyrmon","UEyrmon","EUyrmon","yri","id"))

sipp[ , c("EEyrmon","EUyrmon","UEyrmon","jobchng_ann"):=NULL]

sipp <- merge(sipp,sipp_ann, by=c("id","yri"), all=T)

sipp[ EU ==T, EUmis:= mis ]
sipp[ is.na(EUmis), EUmis:= 0 ]
sipp[ ustintid>0, EUmis:= max(EUmis,na.rm = T), by=list(id,ustintid)]
sipp[ , EUmis:= max(EUmis,na.rm=T), by=list(id,wave)]
sipp[ (EU_wave==T|UE_wave==T)&!(midEU|midUE) & EUmis<=maxmis-12 & matched_EUUE_wave==T, validEUUE:=T]
sipp[ is.na(validEUUE), validEUUE:=F]
sipp[ , EUmis:=NULL]

#sipp[ , ersend_wave:= Mode(ersend), by=list(id,wave)]
sipp[ ustintid_wave>0 , ersend_wave:=Mode(ersend), by=list(id,ustintid_wave)]
sipp[ !is.finite(ustintid_wave), ersend_wave := Mode(ersend), by=list(id,wave)]
sipp[ ustintid_wave>0 , var_ersend:=var(ersend,na.rm = T), by=list(id,ustintid_wave)]
sipp[ var_ersend>0 , ersend_wave:=NA]
sipp[ , var_ersend:=NULL]

sipp[ is.finite(ersend_wave) , displaced:= (ersend_wave>=9 & ersend_wave<=10)|(ersend_wave==1)|(ersend_wave==13)]

sipp[is.finite(ersend_wave), displaced_layoff := ersend_wave==1]
sipp[is.finite(ersend_wave), displaced_empclosed := ersend_wave==9|ersend_wave==10]
sipp[is.finite(ersend_wave), displaced_slackbiz := ersend_wave==13]
########## save prepared data--------------
#a bit of cleanup
sipp[ , c("ajbocc","date0","jobchng_max","EE_max","EU_max","UE_max","PCEPI","last.earnm","last.job","last.job_wave","last.occ","last.lfstat",
		  "smonth","syear","epppnum","next.Estart","ui_r","coc","soc2d"):=NULL]

setwd(outputdir)
#saveRDS(sipp, "./preparedSipp.RData")
saveRDS(sipp, "./DTall_3.RData")


if(final_plots==T){
	# plot transitions time series for sanity check------------------
	EU_wave <- sipp[lfstat_wave==1 &  is.finite(next.lfstat_wave), .(EU_wave = weighted.mean(EU_wave , wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(EU_wave, aes(date, EU_wave, color = panel, group = panel)) +
		geom_point()+ylim(c(0,.035))+
		geom_line() +xlab("") + ylab("EU wave-frequency")
	
	EU_wave <- sipp[lfstat_wave==1 & is.finite(next.lfstat_wave), .(EU_wave = weighted.mean(EU_wave, wpfinwgt, na.rm = TRUE)), by = list(date)]
	setkey(EU_wave,date)
	EU_wave[ , EU_wave_apx:= na.approx(EU_wave,na.rm=F)]
	ggplot(EU_wave, aes(date, EU_wave)) + theme_bw()+
		geom_point()+ylim(c(0.,0.035))+
		geom_line( aes(date,EU_wave_apx)) +xlab("") + ylab("EU wave-frequency")
	ggsave(filename = paste0(figuredir,"/EU_wave.eps"),height= 5,width=10)
	ggsave(filename = paste0(figuredir,"/EU_wave.png"),height= 5,width=10)
	
	
	UE_wave <- sipp[lfstat_wave==2 & is.finite(next.lfstat_wave), .(UE_wave = weighted.mean(UE_wave & !midUE, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(UE_wave, aes(date, UE_wave, color = panel, group = panel)) +
		geom_point() +
		geom_line() +xlab("") + ylab("UE wave-frequency")
	
	UE_wave <- sipp[lfstat_wave==2 & is.finite(next.lfstat_wave), .(UE_wave = weighted.mean(UE_wave, wpfinwgt, na.rm = TRUE)), by = list(date)]
	setkey(UE_wave,date)
	UE_wave[ , UE_wave_apx:= na.approx(UE_wave,na.rm=F)]
	ggplot(UE_wave, aes(date, UE_wave)) + theme_bw()+
		geom_point()+ylim(c(0.,0.5))+
		geom_line( aes(date,UE_wave_apx)  ) +xlab("") + ylab("UE wave-frequency")
	ggsave(filename = paste0(figuredir,"/UE_wave.eps"),height= 5,width=10)
	ggsave(filename = paste0(figuredir,"/UE_wave.png"),height= 5,width=10)
	
	
	EE_wave <- sipp[lfstat_wave==1 & next.lfstat_wave==1 , .(EE_wave = weighted.mean(EE_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(EE_wave, aes(date, EE_wave, color = panel, group = panel)) +
		geom_point() +
		geom_line() +xlab("") + ylab("EE wave-frequency")
	
	EE_wave <- sipp[lfstat_wave==1 & next.lfstat_wave==1 & !(panel==2004 & year<2005), .(EE_wave = weighted.mean(EE_wave, wpfinwgt, na.rm = TRUE)), by = list(date)]
	setkey(EE_wave,date)
	EE_wave[ , EE_wave_apx:= na.approx(EE_wave,na.rm=F)]
	ggplot(EE_wave, aes(date, EE_wave)) + theme_bw()+
		geom_point()+ylim(c(0.,0.035))+
		geom_line( aes(date,EE_wave_apx)  ) +xlab("") + ylab("EE wave-frequency")
	ggsave(filename = paste0(figuredir,"/EE_wave.eps"),height= 5,width=10)
	ggsave(filename = paste0(figuredir,"/EE_wave.png"),height= 5,width=10)
	
	
	
	JC_wave <- sipp[, .(JC_wave = weighted.mean(jobchng_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(JC_wave, aes(date, JC_wave, color = panel, group = panel)) + ylim(0,0.1)+
		geom_point() +
		geom_line()
	
	U_wave <- sipp[ , .(U_wave =  weighted.mean(lfstat_wave==2, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(U_wave, aes(date, U_wave, color = panel, group = panel)) +
		geom_point() + ylim(c(0,.10))+
		geom_line()
	ui_wave <- sipp[ lfstat==2, .(ui_wave = weighted.mean((ui_a>0), wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(ui_wave, aes(date, ui_wave, color = panel, group = panel)) +
		geom_point() + ylim(c(0,.40))+
		geom_line()
	
	
	swOc_wave <- sipp[ (EE_wave==T|EU_wave==T|UE_wave==T) & is.finite(switchedOcc_wave), .(swOc = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(swOc_wave, aes(date, swOc, color = panel, group = panel)) +
		geom_point() + 
		geom_smooth()
	
	swOcEE_wave <- sipp[EE_wave==T & is.finite(switchedOcc_wave) , .(swOcEE = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(swOcEE_wave, aes(date, swOcEE, color = panel)) +
		geom_point() + 
		geom_smooth()
	
	swOcEUUE_wave <- sipp[(EU_wave==T|UE_wave==T) &  is.finite(switchedOcc_wave) , .(swOc = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
	ggplot(swOcEUUE_wave, aes(date, swOc, color = panel)) +
		geom_point() + 
		geom_smooth() + ggtitle("Occupational Switching | EU,UE counted at EU")
	
	 
	swOc_wave <- sipp[EE_wave==T& is.finite(occ_wave) & is.finite(next.occ_wave) & !(panel=="2004" & (year<2005 | year>=2007)) , .(swOc_wave = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), by = date]
	setkey(swOc_wave,date)
	swOc_wave[ , swOc_wave_apx:= na.approx(swOc_wave,na.rm=F)]
	ggplot(swOc_wave, aes(date, swOc_wave)) +
		geom_point()  + theme_bw()+ylim(c(0.45,0.6))+
		geom_line(aes(date,swOc_wave_apx)) + xlab("") + ylab("Occupational Switching Frequency")
	ggsave(filename = paste0(figuredir,"/sWocEE_wave.eps"),height= 5,width=10)
	ggsave(filename = paste0(figuredir,"/sWocEE_wave.png"),height= 5,width=10)
	
	misSwOcEU_wave <- sipp[ EU_wave==T & matched_EUUE_wave==T, .(misSwOcEU_wave = mean(is.na(switchedOcc_wave))), by=list(panel,date)]
	ggplot( misSwOcEU_wave, aes(date, misSwOcEU_wave, color=panel, group=panel) ) +
		geom_point() +
		geom_line()
	misSwOcUE_wave <- sipp[ UE_wave==T , .(misSwOcUE_wave = mean(is.na(switchedOcc_wave))), by=list(panel,date)]
	ggplot( misSwOcUE_wave, aes(date, misSwOcUE_wave, color=panel, group=panel) ) +
		geom_point() +
		geom_line()
	
	misSwOcUE_wave <- sipp[matched_EUUE_wave==T & UE_wave==T , .(misSwOcUE_wave = mean(is.na(switchedOcc_wave))), by=list(panel,date)]
	ggplot( misSwOcUE_wave, aes(date, misSwOcUE_wave, color=panel, group=panel) ) +
		geom_point() +
		geom_line()
	
	
	sipp_sw_state_EE <- sipp[ EE_wave==T, .(swEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw_state_EU <- sipp[ EU_wave==T, .(swEU=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw_state_EE,sipp_sw_state_EU,by= "state")
	sipp_sw_state_EUEE <- sipp[ EE_wave==T|EU_wave==T, .(swEUEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw,sipp_sw_state_EUEE,by= "state")
	sipp_sw19962013 <- merge(sipp_sw,state_codes,by.x= "state", by.y="code")
	
	#1998-2000
	sipp_sw_state_EE <- sipp[ EE_wave==T & year>=1998 & year<=2000, .(swEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw_state_EU <- sipp[ EU_wave==T & year>=1998 & year<=2000, .(swEU=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw_state_EE,sipp_sw_state_EU,by= "state")
	sipp_sw_state_EUEE <- sipp[ EE_wave==T|EU_wave==T & year>=1998 & year<=2000, .(swEUEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw,sipp_sw_state_EUEE,by= "state")
	sipp_sw2000 <- merge(sipp_sw,state_codes,by.x= "state", by.y="code")
	#
	sipp_sw_state_EE <- sipp[ EE_wave==T & year>=2011 & year<=2013, .(swEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw_state_EU <- sipp[ EU_wave==T & year>=2011 & year<=2013, .(swEU=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw_state_EE,sipp_sw_state_EU,by= "state")
	sipp_sw_state_EUEE <- sipp[ EE_wave==T|EU_wave==T & year>=2011 & year<=2013, .(swEUEE=wtd.mean(switchedOcc_wave,na.rm = T,weights=wpfinwgt)), by=state ]
	sipp_sw <- merge(sipp_sw,sipp_sw_state_EUEE,by= "state")
	sipp_sw2013 <- merge(sipp_sw,state_codes,by.x= "state", by.y="code")
}