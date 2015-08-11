# March 26, 2015
# Prepare SIPP data for analysis
# Fill in missing occupation codes, generate dummy variables for 
# switching occupations and labor force flow, add SOC2d codes.
# Save processed files to ./Data/ directory.
# Precondition: readDTAs.R has been run.
library(dplyr)
library(zoo)
library(foreign)
library(stats)
library(reshape2)

setwd("~/workspace/CVW/R")
#setwd("G:/Research_Analyst/Eubanks/Occupation switching")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T

# Read crosswalk files
coc2000_to_occ1990 <- read.dta("./Crosswalks/coc2000_2_occ1990.dta")
occ1990_to_SOC2d <- read.dta("./Crosswalks/occ90_2_soc2d.dta", convert.underscore = TRUE) %>%
        select(-.merge.occs)

setwd("./Data")

# Functions ---------------------------------------------------------------

# Generate recession indicator
genRec <- function(df) {
	rec_dates <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01"))
	result <- df %>%
		mutate(recIndic = (date > rec_dates[1] & date < rec_dates[2]) |
			   	(date > rec_dates[3] & date < rec_dates[4]))
	result <- group_by(result, id, wave) %>%
		mutate(waveRec = as.logical(max(recIndic))) %>%
		ungroup
	return(result)
}

# Generate LF status variable from esr
genLFStat <- function(df) {#
	result <- df %>%
		# 1: employed
		mutate(lfStat = ifelse(esr >= 1 & esr <= 5, 1, 0)) %>%
		# 2: unemployed
		mutate(lfStat = ifelse(esr == 6 | esr == 7, 2, lfStat)) %>%
		# 3: NILF
		mutate(lfStat = ifelse(esr == 8, 3, lfStat))
		# unemp indicator
		# mutate(unemployed = (lfStat == 2))
		return(result)
}

# Correct occupation code
fixOccCode <- function(df) {
	if(useSoc2d) {
		# drop occ and replace with soc2d
		df <- df %>%
			# this is to fill soc2d if the occupation merge didn't work
			mutate(soc2d = ifelse( is.na(soc2d) & !is.na(occ) & (occ>0 & occ<=900), 99, soc2d ) ) %>%
			select(-occ) %>%
			mutate(occ = soc2d)
	}
	result <- df %>%
		group_by(id) %>%
		arrange(date) %>%
		#                 # replace occ with NA if unemployed or NILF
		#                 mutate(occ = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, occ))) %>%
		#                 # carry forward last observation of occ to fill NAs
		#                 mutate(occ = na.locf(occ, na.rm = FALSE, fromLast = TRUE)) %>%
		# replace NA job codes with 0
		mutate(job = as.integer(ifelse(is.na(job), 0, job))) %>%
		# replace job code with 0 if unemployed or NILF 
		mutate(job = as.integer(ifelse(lfStat == 2 | lfStat == 3, 0, job))) %>%
		# replace ind23 with NA if unemployed or NILF
		mutate(ind23 = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, ind23))) %>%
		# carry backward last observation of ind23 to fill NAs
		mutate(ind23 = na.locf(ind23, na.rm = FALSE, fromLast = TRUE)) %>%
		# replace occupation codes with NA if unemp or NILF	
		mutate(occ = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, occ))) %>%
		# carry backward last observation of ind23 to fill NAs
		mutate(occ = na.locf(occ, na.rm = FALSE, fromLast = TRUE))
	return(result)
}


# Generate occupation switching and LF flow dummies
genFlowDummies <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(date) %>%
		mutate(switchedJob = (job != lead(job)) & !(is.na(job) | is.na(lead(job))) ) %>%
		mutate(EE = lfStat == 1 & lead(lfStat) == 1 & switchedJob) %>%
		# now EU
		mutate(UE = lfStat == 1 & lead(lfStat) == 2 & switchedJob) %>%
		mutate(switchedOcc = (occ != lead(occ)) & switchedJob) %>%
		mutate(switchedInd = (ind23 != lead(ind23)) & switchedJob)         
	return(result)
}

# Generate unemployment duration
# If respondent enters panel unemployed, duration will be NA for that spell
# can only call after genLFStat
genUnempDuration <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(date) %>%
		# generate dummy for unemployed
		mutate(unemployed = lfStat >= 2)
	result <- result %>%
		# generate unique id for each period respondent enters unemployment
		mutate(spellID = as.integer(ifelse(unemployed & !lag(unemployed), 1:n(), NA))) %>%
		# carryforward unique id
		mutate(spellID = na.locf(spellID, na.rm = FALSE)) %>%
		group_by(id, spellID) %>%
		# in each spell, calculate cumulative sum of unemployed
		mutate(unempDur = as.integer(ifelse(unemployed & !is.na(spellID), cumsum(unemployed), 0))) %>%
		mutate(unempDur = as.integer(ifelse(unemployed & !is.na(spellID), max(unempDur), 0))) %>%
		# THIS SEEMS NOT TO WORK: fill it so that this is in the period prior to a spell
		mutate(unempDur = as.integer(ifelse(lfStat == 1 & lead(unemployed), lead(unempDur), unempDur))) %>%
		mutate(unempDur = as.integer(ifelse(lfStat == 1 & !lead(unemployed), 0., unempDur))) %>%
		ungroup %>%
		select(-spellID, -unemployed) %>%
		ungroup
	return(result)
}

# Remove data from non-interview months to keep only the "seams"
# 4th observation in each wave-id is srefmon == 4, interview month
collapseSeams <- function(df) {
	result <- df %>%
		group_by(id, wave) %>%
		summarize(age = nth(age, 4),
				  educ = nth(educ, 4),
				  female = nth(female, 4),
				  race = nth(race, 4),
				  wpfinwgt = nth(wpfinwgt, 4),
				  earnm = nth(earnm, 4),
				  occ = nth(occ, 4),
				  soc2d = nth(soc2d, 4),
				  job = nth(job, 4),
				  esr = nth(esr, 4),
				  occ14 = nth(occ14, 4),
				  ind23 = nth(ind23, 4),
				  srefmon = nth(srefmon, 4),
				  date = nth(date, 4),
				  qtrdate = nth(qtrdate, 4),
				  #lfStat = nth(lfStat, 4),
				  recIndic = as.logical(max(recIndic)),
				  lfStat = max(lfStat)
		)
	return(result)     
}

# Remove NILFs and NAs
sampleSelect <- function(df) {
	result <- df %>%
		filter(age >= 18 & age <= 65) %>%
		filter(!is.na(job)) %>%
		filter(!is.na(lead(job))) %>%
#		filter(!is.na(occ)) %>%
#		filter(!is.na(lead(occ)))
#		filter(!is.na(esr)) %>%
#		filter(!is.na(occ)) %>%
		filter(!is.na(lfStat)) %>%
		filter(lfStat != 3) %>%
		filter(!(is.na(earnm) & lfStat==1)) %>%
		filter(!(earnm <0 & lfStat==1))
	return(result)
}

# 1996 Panel --------------------------------------------------------------

sipp96 <- readRDS("sipp96.RData")

# add soc2d codes
processed96 <- sipp96 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000, 
                                          occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occDiv" = "occ1990")) %>%
        select(-occDiv, -occ2000)

# generate variables for analysis
processed96 <- processed96 %>%
        genRec(.) %>%
        genLFStat(.) %>%                                        
        fixOccCode(.) %>%                                      
        genFlowDummies(.) %>%                           
        genUnempDuration(.) %>%
		sampleSelect(.)
processed96 <- processed96 %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(unempDur = as.numeric( ifelse(UE, lead(unempDur), unempDur) ) )%>%
	mutate(unempDur = as.numeric( ifelse(EE, 0., unempDur) )) %>%
	ungroup

if(useSoc2d) {
        setwd("./soc2d")
} else {
        setwd("./occ")
}
saveRDS(processed96, "processed96.RData")

rm(list = c("sipp96", "processed96"))

setwd("..")

# 2001 Panel --------------------------------------------------------------

sipp01 <- readRDS("sipp01.RData")

# add soc2d codes
processed01 <- sipp01 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occDiv" = "occ1990")) %>%
        select(-occDiv, -occ2000)

# generate variables for analysis
processed01 <- processed01 %>%        
	genRec(.) %>%
	genLFStat(.) %>%
	fixOccCode(.) %>%
	genFlowDummies(.) %>%
	genUnempDuration(.) %>%
	sampleSelect(.)
processed01 <- processed01 %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(unempDur = as.numeric( ifelse(UE, lead(unempDur), unempDur) ) )%>%
	mutate(unempDur = as.numeric( ifelse(EE, 0., unempDur) )) %>%
	ungroup

if(useSoc2d) {
        setwd("./soc2d")
} else {
        setwd("./occ")
}
saveRDS(processed01, "processed01.RData")

rm(list = c("sipp01", "processed01"))

setwd("..")

# 2004 Panel --------------------------------------------------------------

sipp04 <- readRDS("sipp04.RData")

# add soc2d codes
processed04 <- sipp04 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(coc2000_to_occ1990, 
                  by = c("occDiv" = "coc2000")) %>%
        left_join(occ1990_to_SOC2d, by = "occ1990") %>%
        select(-occDiv, -occ2000, -occ1990)
        
# generate variables for analysis
processed04<- processed04 %>%
	genRec(.) %>%
	genLFStat(.) %>%
	fixOccCode(.) %>%
	genFlowDummies(.) %>%
	genUnempDuration(.) %>%
	sampleSelect(.)
processed04 <- processed04 %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(unempDur = as.numeric( ifelse(UE, lead(unempDur), unempDur) ) )%>%
	mutate(unempDur = as.numeric( ifelse(EE, 0., unempDur) )) %>%
	ungroup

if(useSoc2d) {
        setwd("./soc2d")
} else {
        setwd("./occ")
}
saveRDS(processed04, "processed04.RData")

rm(list = c("sipp04", "processed04"))

setwd("..")

# 2008 Panel --------------------------------------------------------------

sipp08 <- readRDS("sipp08.RData")

# add soc2d codes
processed08 <- sipp08 %>%
	mutate(occDiv = as.integer(ifelse(occ >= 1000, 
									  occ/10, occ))) %>%
	left_join(coc2000_to_occ1990, 
			  by = c("occDiv" = "coc2000")) %>%
	left_join(occ1990_to_SOC2d, by = "occ1990") %>%
	select(-occDiv, -occ2000, -occ1990)

# generate variables for analysis
processed08 <- processed08 %>%
	genRec(.) %>%
	genLFStat(.) %>%
	fixOccCode(.) %>%
	genFlowDummies(.) %>%
	genUnempDuration(.) %>%
	sampleSelect(.)
processed08 <- processed08 %>%
	group_by(id) %>%
	arrange(id, date) %>%
	mutate(unempDur = as.numeric( ifelse(UE, lead(unempDur), unempDur) ) )%>%
	mutate(unempDur = as.numeric( ifelse(EE, 0., unempDur) )) %>%
	ungroup

if(useSoc2d) {
	setwd("./soc2d")
} else {
	setwd("./occ")
}
saveRDS(processed08, "processed08.RData")

rm(list = c("sipp08", "processed08"))