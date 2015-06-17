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

# Functions ---------------------------------------------------------------

# Generate LF status variable from esr
genLFStat <- function(df) {#
        result <- df %>%
                # 1: employed
                #mutate(lfStat = ifelse(esr == 1 | esr == 2 | esr == 4 | esr == 3 | esr == 5, 1, 0)) %>%
        		mutate(lfStat = ifelse(esr >= 1 & esr <= 5, 1, 0)) %>%
                # 2: unemployed
                mutate(lfStat = ifelse(esr == 6 | esr == 7, 2, lfStat)) %>%
                # 3: NILF
                mutate(lfStat = ifelse(esr == 8, 3, lfStat))
        return(result)
}

# Correct occupation code
fixOccCode <- function(df) {
	if(useSoc2d) {
		# drop occ and replace with soc2d
		df <- df %>%
			select(-occ) %>%
			rename(occ = soc2d)
	}
	result <- df %>%
		group_by(id) %>%
		arrange(date) %>%
		# replace occ with NA if unemployed or NILF
		mutate(occ = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, occ))) %>%
		# carry forward last observation of occ to fill NAs
		mutate(occ = na.locf(occ, na.rm = FALSE)) %>%
		# replace NA job codes with 0
		mutate(job = as.integer(ifelse(is.na(job), 0, job))) %>%
		# replace job code with 0 if unemployed or NILF 
		mutate(job = as.integer(ifelse(lfStat == 2 | lfStat == 3, 0, job)))
	return(result)
}

# Generate occupation switching and LF flow dummies
genFlowDummies <- function(df) {
        result <- df %>%
                group_by(id) %>%
                arrange(date) %>%
                mutate(switchedJob = job != lead(job)) %>%
                mutate(switchedOcc = (occ != lead(occ)) & switchedJob &
                	   	!is.na(occ) & !is.na(lead(occ)) ) %>%
                mutate(EE = lfStat == 1 & lead(lfStat) == 1 & switchedJob &
                               !is.na(occ) & !is.na(lead(occ)) ) %>%
                mutate(UE = lfStat == 2 & lead(lfStat) == 1 & switchedJob &
                               !is.na(occ) & !is.na(lead(occ)))
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
                mutate(unemployed = lfStat == 2)
        result <- result %>%
                # generate unique id for each period respondent enters unemployment
                mutate(spellID = as.integer(ifelse(unemployed & !lag(unemployed), 1:n(), NA))) %>%
                # carryforward unique id
                mutate(spellID = na.locf(spellID, na.rm = FALSE)) %>%
                group_by(id, spellID) %>%
                # in each spell, calculate cumulative sum of unemployed
                mutate(unempDur = as.integer(ifelse(unemployed & !is.na(spellID), cumsum(unemployed), NA))) %>%
                select(-spellID)
        return(result)
}

# 1996 Panel --------------------------------------------------------------

sipp96 <- readRDS("./Data/sipp96.RData")

# add soc2d codes
processed96 <- sipp96 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000, 
                                          occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occDiv" = "occ1990")) %>%
        select(-occDiv, -occ2000)

# generate variables for analysis
processed96 <- processed96 %>%
        genLFStat(.) %>%                                        
        fixOccCode(.) %>%                                      
        genFlowDummies(.) %>%                           
        genUnempDuration(.)

if(useSoc2d) {
        saveRDS(processed96, "./Data/processed96soc2d.RData")
} else {
        saveRDS(processed96, "./Data/processed96.RData")
}
rm(list = c("sipp96", "processed96"))

# 2001 Panel --------------------------------------------------------------

sipp01 <- readRDS("./Data/sipp01.RData")

# add soc2d codes
processed01 <- sipp01 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occDiv" = "occ1990")) %>%
        select(-occDiv, -occ2000)

# generate variables for analysis
processed01 <- processed01 %>%        
        genLFStat(.) %>%
        fixOccCode(.) %>%
        genFlowDummies(.) %>%
        genUnempDuration(.)

if(useSoc2d) {
        saveRDS(processed01, "./Data/processed01soc2d.RData")
} else {
        saveRDS(processed01, "./Data/processed01.RData")
}
rm(list = c("sipp01", "processed01"))

# 2004 Panel --------------------------------------------------------------

sipp04 <- readRDS("./Data/sipp04.RData")

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
        genLFStat(.) %>%
        fixOccCode(.) %>%
        genFlowDummies(.) %>%
        genUnempDuration(.)

if(useSoc2d) {
        saveRDS(processed04, "./Data/processed04soc2d.RData")
} else {
        saveRDS(processed04, "./Data/processed04.RData")
}
rm(list = c("sipp04", "processed04"))

# 2008 Panel --------------------------------------------------------------

sipp08 <- readRDS("./Data/sipp08.RData")

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
        genLFStat(.) %>%
        fixOccCode(.) %>%
        genFlowDummies(.) %>%
        genUnempDuration(.)

if(useSoc2d) {
        saveRDS(processed08, "./Data/processed08soc2d.RData")
} else {
        saveRDS(processed08, "./Data/processed08.RData")
}
rm(list = c("sipp08", "processed08"))