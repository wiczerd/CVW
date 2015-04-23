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

# Read crosswalk files
coc2000_to_occ1990 <- read.dta("./Crosswalks/coc2000_2_occ1990.dta")
occ1990_to_SOC2d <- read.dta("./Crosswalks/occ90_2_soc2d.dta", convert.underscore = TRUE) %>%
        select(-.merge.occs)

# Generate LF status variable from esr
genLFStat <- function(df) {#
        # 1: employed
        # 2: unemployed
        # 3: NILF
        mutate(df, 
               lfStat = ifelse(esr == 1 | esr == 2 | esr == 4, 1, 0),
               lfStat = ifelse(esr == 3 | esr == 5 | esr == 6 | esr == 7, 2, lfStat),
               lfStat = ifelse(esr == 8, 3, lfStat))
}

# Correct occupation code
fixOccCode <- function(df) {
        df %>%
                group_by(id) %>%
                arrange(date) %>%
                mutate(occ = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, occ))) %>%
        		mutate(occ = na.locf(occ, na.rm = FALSE)) %>%
        		mutate(job = as.integer(ifelse(is.na(job), 0, job))) %>%
        		mutate(job = as.integer(ifelse(lfStat == 2 | lfStat == 3, 0, job))) 
}

# Generate occupation switching and LF flow dummies
# df must be grouped and arranged
# genFlowDummies <- function(df) {
#         mutate(df,
#                switchedJob = job != lead(job),
#                switchedOcc = (occ != lead(occ)),
#                EE = lfStat == 1 & lead(lfStat) == 1 & switchedJob &
#                        !is.na(occ) & !is.na(lead(occ)),
#                UE = lfStat == 2 & lead(lfStat) == 1 & switchedJob &
#                        !is.na(occ) & !is.na(lead(occ)))
#                #EEraw = lfStat == 1 & lead(lfStat) == 1 & switchedJob,
#                #UEraw = lfStat == 2 & lead(lfStat) == 1 & switchedJob)
# }
genFlowDummies <- function(df) {
	df$switchedJob <- df$job != lead(df$job)
	df$switchedOcc <- (df$occ != lead(df$occ))
	df$EE <- df$lfStat == 1 & lead(df$lfStat) == 1 & df$switchedJob &	
		!is.na(df$occ) & !is.na(lead(df$occ))
	df$UE <- df$lfStat == 2 & lead(df$lfStat) == 1 & df$switchedJob &
		!is.na(df$occ) & !is.na(lead(df$occ))
	return(df)     

}



# 1996
sipp96 <- readRDS("./Data/sipp96.RData")
processed96 <- sipp96 %>%
        genLFStat(.) %>%                                        # generate LF status variable
        fixOccCode(.) %>%                                       # fix occupation codes
        genFlowDummies(.) %>%                                   # generate flow dummies
        mutate(occ = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occ" = "occ1990")) %>%                # add SOC codes
        group_by(id) %>%                                        # group by id
        arrange(date) %>%                                       # sort by date witin id
        mutate(switched2d = (soc2d != lead(soc2d))) %>%         # generate adjusted switch dummies
        select(-occ2000)
saveRDS(processed96, "./Data/processed96.RData")
rm(list = c("sipp96", "processed96"))

# 2001
sipp01 <- readRDS("./Data/sipp01.RData")
processed01 <- sipp01 %>%
        genLFStat(.) %>%                                        # generate LF status variable
        fixOccCode(.) %>%                                       # fix occupation code
        genFlowDummies(.) %>%                                   # generate flow dummies
        mutate(occ = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, 
                  by = c("occ" = "occ1990")) %>%                # add SOC codes
        group_by(id) %>%                                        # group by id
        arrange(date) %>%                                       # sort by date witin id
        mutate(switched2d = (soc2d != lead(soc2d))) %>%         # generate adjusted switch dummies
        select(-occ2000)
saveRDS(processed01, "./Data/processed01.RData")
rm(list = c("sipp01", "processed01"))

# 2004
sipp04 <- readRDS("./Data/sipp04.RData")
processed04 <- sipp04 %>%
        genLFStat(.) %>%                                        # generate LF status variable
        fixOccCode(.) %>%                                       # fix occupation code
        genFlowDummies(.) %>%                                   # generate flow dummies
        mutate(occ = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(coc2000_to_occ1990, 
                  by = c("occ" = "coc2000")) %>%                # convert codes to 1990
        left_join(occ1990_to_SOC2d, by = "occ1990") %>%         # add SOC codes
        group_by(id) %>%                                        # group by id
        arrange(date) %>%                                       # sort by date witin id
        mutate(switched2d = (soc2d != lead(soc2d))) %>%         # generate adjusted switch dummies
        select(-occ1990, -occ2000)
saveRDS(processed04, "./Data/processed04.RData")
rm(list = c("sipp04", "processed04"))

# 2008
sipp08 <- readRDS("./Data/sipp08.RData")
processed08 <- sipp08 %>%
        genLFStat(.) %>%                                        # generate LF status variable
        fixOccCode(.) %>%                                       # fix occupation code
        genFlowDummies(.) %>%                                   # generate flow dummies
        mutate(occ = as.integer(ifelse(occ >= 1000, 
                                       occ/10, occ))) %>%
        left_join(coc2000_to_occ1990, 
                  by = c("occ" = "coc2000")) %>%                # convert codes to 1990
        left_join(occ1990_to_SOC2d, by = "occ1990") %>%         # add SOC codes
        group_by(id) %>%                                        # group by id
        arrange(date) %>%                                       # sort by date witin id
        mutate(switched2d = (soc2d != lead(soc2d))) %>%         # generate adjusted switch dummies
        select(-occ1990, -occ2000)
saveRDS(processed08, "./Data/processed08.RData")
rm(list = c("sipp08", "processed08"))