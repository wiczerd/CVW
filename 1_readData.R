# December 17, 2015
# Read in data and add basic variables
# 1) read data into R
# 2) select variables of interest
# 3) create date variable
# 4) create recession indicator
# 5) add PCE data
# 6) add SOC2d codes
library(foreign)
library(data.table)

# eventually move data files to master directory
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

# specify which variables to keep from CEPR sets A, B, and D
keepVars <- c("age", 
			  "educ", 
			  "female", 
			  "id",
			  "race", 
			  "month",
			  "srefmon",
			  "wpfinwgt", 
			  "year", 
			  "earnm",
			  "occ", 
			  "job", 
			  "esr","ersend","estlemp","rwkesr2",
			  "union",
			  "shhadid",
			  #"occ14", 
			  "ind", "ind23",
			  #"state", 
			  "wave")

# create vector of recession dates
recDates <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01"))

# get PCE data
PCE <- readRDS("InputData/PCE.RData")

# get crosswalk data
occ90_soc2d <- readRDS(paste0(xwalkdir,"/occ90_soc2d.RData"))
coc2000_occ1990 <- readRDS(paste0(xwalkdir,"/coc2000_occ1990.RData"))
coc1990_occ1990 <- readRDS(paste0(xwalkdir,"/coc1990_occ1990.RData"))


CIC2002_CIC2000<- readRDS(paste0(xwalkdir,"/CIC2002_2_CIC2000.RData"))
CIC2000_2_CIC1990<- readRDS(paste0(xwalkdir,"/CIC2000_2_CIC1990.RData"))


# 1996 panel --------------------------------------------------------------

sipp96 <- read.dta("./Data/sippsets96ABD.dta", convert.factors = FALSE)
sipp96 <- sipp96[keepVars]
DT96 <- data.table(sipp96)
rm(sipp96)

# create date variable, remove year and month to save space
DT96[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]
DT96[, c("year", "month") := NULL]

# create recession indicator
DT96[, recIndic := (date > recDates[1] & date < recDates[2]) | 
     		   (date > recDates[3] & date < recDates[4])]

# add PCE data
DT96 <- merge(DT96, PCE, by = "date", all.x = TRUE)

# add soc2d codes
DT96[, coc90 := occ]
DT96 <- merge(DT96, coc1990_occ1990, by  = "coc90", all.x = TRUE)
DT96 <- merge(DT96, occ90_soc2d, by  = "occ1990", all.x = TRUE)
DT96[, c("occ","coc90") := NULL]
setnames(DT96, "occ1990", "occ")

#recode industry as ind23
DT96[, ind:= ind23]
DT96[, ind23:=NULL]

saveRDS(DT96, file("./Data/DT96_1.RData"))
rm(DT96)

# 2001 panel --------------------------------------------------------------

sipp01 <- read.dta("./Data/sippsets01ABD.dta", convert.factors = FALSE)
sipp01 <- sipp01[keepVars]
DT01 <- data.table(sipp01)
rm(sipp01)

# create date variable, remove year and month to save space
DT01[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]
DT01[, c("year", "month") := NULL]

# create recession indicator
DT01[, recIndic := (date > recDates[1] & date < recDates[2]) | 
     		   (date > recDates[3] & date < recDates[4])]

# add PCE data
DT01 <- merge(DT01, PCE, by = "date", all.x = TRUE)

# add unemployment data
# FILL IN LATER

# add soc2d codes
DT01[, coc90 := occ]
DT01 <- merge(DT01, coc1990_occ1990, by  = "coc90", all.x = TRUE)
DT01 <- merge(DT01, occ90_soc2d, by  = "occ1990", all.x = TRUE)
DT01[, c("occ","coc90") := NULL]
setnames(DT01, "occ1990", "occ")
#recode industry as ind23
DT01[, ind:= ind23]
DT01[, ind23:=NULL]

saveRDS(DT01, file("./Data/DT01_1.RData"))
rm(DT01)

# 2004 panel --------------------------------------------------------------

sipp04 <- read.dta("./Data/sippsets04ABD.dta", convert.factors = FALSE)
sipp04 <- sipp04[keepVars]
DT04 <- data.table(sipp04)
rm(sipp04)

# create date variable, remove year and month to save space
DT04[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]
DT04[, c("year", "month") := NULL]

# create recession indicator
DT04[, recIndic := (date > recDates[1] & date < recDates[2]) | 
     		   (date > recDates[3] & date < recDates[4])]

# add PCE data
DT04 <- merge(DT04, PCE, by = "date", all.x = TRUE)

# add unemployment data
# FILL IN LATER

# add soc2d codes
DT04[, coc2000 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
DT04 <- merge(DT04, coc2000_occ1990, by  = "coc2000", all.x = TRUE)
DT04 <- merge(DT04, occ90_soc2d, by  = "occ1990", all.x = TRUE)
DT04[, c("occ", "coc2000") := NULL]
setnames(DT04, "occ1990", "occ")

#add conversion to ind23
#DT04 <- merge(DT04,CIC2002_CIC2000, by = "ind", all.x=T)
#DT04[ , ind := CIC2000]
DT04[ , ind := as.integer(ind/10)]
DT04[ , ind23 := NULL]
DT04 <- merge(DT04,CIC2000_2_CIC1990, by = "ind", all.x=T)
DT04[ , ind := ind23]
DT04[, ind23:=NULL]

saveRDS(DT04, file("./Data/DT04_1.RData"))
rm(DT04)

# 2008 panel --------------------------------------------------------------

sipp08 <- read.dta("./Data/sippsets08ABD.dta", convert.factors = FALSE)
sipp08 <- sipp08[keepVars]
DT08 <- data.table(sipp08)
rm(sipp08)

# create date variable, remove year and month to save space
DT08[, date := as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")]
DT08[, c("year", "month") := NULL]

# create recession indicator
DT08[, recIndic := (date > recDates[1] & date < recDates[2]) | 
     		   (date > recDates[3] & date < recDates[4])]

# add PCE data
DT08 <- merge(DT08, PCE, by = "date", all.x = TRUE)

# add unemployment data
# FILL IN LATER

# add soc2d codes
DT08[, coc2000 := as.integer(ifelse(occ >= 1000,  occ/10, occ))]
DT08 <- merge(DT08, coc2000_occ1990, by  = "coc2000", all.x = TRUE)
DT08 <- merge(DT08, occ90_soc2d, by  = "occ1990", all.x = TRUE)
DT08[, c("occ", "coc2000") := NULL]
setnames(DT08, "occ1990", "occ")

#add conversion to ind23
#DT08 <- merge(DT08,CIC2002_CIC2000, by = "ind", all.x=T)
#DT08[ , ind := CIC2000]
DT08[ , ind := as.integer(ind/10)]
DT08[ , ind23 := NULL]
DT08 <- merge(DT08,CIC2000_2_CIC1990, by = "ind", all.x=T)
DT08[ , ind := ind23]
DT08[, ind23:=NULL]

saveRDS(DT08, file("./Data/DT08_1.RData"))
rm(list=ls())
