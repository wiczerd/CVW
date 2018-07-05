
# Read in SIPP data and prepare it for analysis
library(foreign)
library(readstata13)
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
datadir <- paste(rootdir, "InputDataDE", sep = "")
outputdir <- paste(rootdir, "Results", sep = "")
figuredir <- paste(rootdir, "Figures", sep = "")
intermed_plots = T
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
			"occ","ind", #"ajbocc",
			# income variables
			"earnm","earn_imp","ui_a","ui_r","thearn","thtotinc"
)

########## read in individual panels, extract variables, and subset sample

setwd(datadir)

# 1990 panel
sipp90 <- read.dta13("./sippsets90ABDFG.dta", convert.factors = FALSE)
sipp90 <- sipp90[toKeep]
sipp90 <- subset(sipp90, !is.na(esr) & (age >= 16) & (age <= 65))
sipp90 <- data.table(sipp90)
sipp90$panel <- 1990

# 1991 panel
sipp91 <- read.dta13("./sippsets91ABDFG.dta", convert.factors = FALSE)
sipp91 <- sipp91[toKeep]
sipp91 <- subset(sipp91, !is.na(esr) & (age >= 16) & (age <= 65))
sipp91 <- data.table(sipp91)
sipp91$panel <- 1991

# 1992 panel
sipp92 <- read.dta13("./sippsets92ABDFG.dta", convert.factors = FALSE)
sipp92 <- sipp92[toKeep]
sipp92 <- subset(sipp92, !is.na(esr) & (age >= 16) & (age <= 65))
sipp92 <- data.table(sipp92)
sipp92$panel <- 1992

# 1993 panel
sipp93 <- read.dta13("./sippsets93ABDFG.dta", convert.factors = FALSE)
sipp93 <- sipp93[toKeep]
sipp93 <- subset(sipp93, !is.na(esr) & (age >= 16) & (age <= 65))
sipp93 <- data.table(sipp93)
sipp93$panel <- 1993


# 1996 panel
sipp96 <- read.dta13("./sippsets96ABDFG.dta", convert.factors = FALSE)#read.dta("./sippsets96ABDFG.dta", convert.factors = FALSE)
sipp96 <- sipp96[toKeep]
sipp96 <- subset(sipp96, !is.na(esr) & (age >= 16) & (age <= 65))
sipp96 <- data.table(sipp96)
sipp96$panel <- 1996

# 2001 panel
sipp01 <- read.dta13("./sippsets01ABDFG.dta", convert.factors = FALSE)
sipp01 <- sipp01[toKeep]
sipp01 <- subset(sipp01, !is.na(esr) & (age >= 16) & (age <= 65))
sipp01 <- data.table(sipp01)
sipp01$panel <- 2001

# 2004 panel
sipp04 <- read.dta13("./sippsets04ABDFG.dta", convert.factors = FALSE)
sipp04 <- sipp04[toKeep]
sipp04 <- subset(sipp04, !is.na(esr) & (age >= 16) & (age <= 65) )
sipp04 <- data.table(sipp04)
sipp04$panel <- 2004

# 2008 panel
sipp08 <- read.dta13("./sippsets08ABDFG.dta", convert.factors = FALSE)
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
coc1980_occ1990 <- readRDS("./coc1980_occ1990.RData")

# 1990 panel
sipp90[, coc:= occ]
sipp90[, coc80 := occ]
sipp90 <- merge(sipp90, coc1980_occ1990, by  = "coc80", all.x = TRUE)
sipp90[ occ1990>=999 , occ1990:=NA ]
sipp90 <- merge(sipp90, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp90[, c("occ", "coc80") := NULL]
setnames(sipp90, "occ1990", "occ")

# 1991 panel
sipp91[, coc:= occ]
sipp91[, coc80 := occ]
sipp91 <- merge(sipp91, coc1980_occ1990, by  = "coc80", all.x = TRUE)
sipp91[ occ1990>=999 , occ1990:=NA ]
sipp91 <- merge(sipp91, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp91[, c("occ", "coc80") := NULL]
setnames(sipp91, "occ1990", "occ")

# 1992 panel
sipp92[, coc:= occ]
sipp92[, coc80 := occ]
sipp92 <- merge(sipp92, coc1980_occ1990, by  = "coc80", all.x = TRUE)
sipp92[ occ1990>=999 , occ1990:=NA ]
sipp92 <- merge(sipp92, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp92[, c("occ", "coc80") := NULL]
setnames(sipp92, "occ1990", "occ")

# 1993 panel
sipp93[, coc:= occ]
sipp93[, coc80 := occ]
sipp93 <- merge(sipp93, coc1980_occ1990, by  = "coc80", all.x = TRUE)
sipp93 <- merge(sipp93, occ1990_soc2d, by  = "occ1990", all.x = TRUE)
sipp93[, c("occ", "coc80") := NULL]
setnames(sipp93, "occ1990", "occ")

# 1996 panel
sipp96[, coc:= occ]
sipp96[, coc90 := occ]
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

sipp <- rbind(sipp90,sipp91,sipp92,sipp93,sipp96, sipp01, sipp04, sipp08)
rm(list=c("sipp90","sipp91","sipp92","sipp93","sipp96", "sipp01", "sipp04", "sipp08"))
sipp <- data.table(sipp)
sipp$panel <- as.factor(sipp$panel)

### convert some types to save space
sipp[ , occ:=as.integer(occ)]
sipp[ , soc2d:=as.integer(soc2d)]
sipp[ , educ:=as.integer(educ)]
sipp[ , state:=as.integer(state)]
sipp[ , female:= as.logical(female)]



# save intermediate result
#setwd(outputdir)
saveRDS(sipp, paste0(outputdir,"/sipp.RData"))
setwd(rootdir)
