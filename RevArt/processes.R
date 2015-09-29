# July 22, 2015
# Process data by calling functions in functions.R
library(dplyr)
library(foreign)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching Lab/")
setwd('~/workspace/CVW/R')
source("./RevArt/functions.R")

useSoc2d <- TRUE

# Wrapper function: runs all processing functions in correct order
processWrapper <- function(df) {
	result <- df %>%
		sampleSelect %>%
		genRec %>%
		genLFStat %>%
		cleanEarn %>%
		fixOccCode %>%
		genFlowDummies %>%
		nextLastOcc
	return(result)
}

# import PCE deflator
PCE <- read.csv("./Data/PCE.csv")
PCE$date <- as.Date(PCE$date)


# Read crosswalk files
coc2000_to_occ1990 <- read.dta("./Crosswalks/coc2000_2_occ1990.dta")
occ1990_to_SOC2d <- read.dta("./Crosswalks/occ90_2_soc2d.dta", convert.underscore = TRUE) %>%
        select(-.merge.occs)

setwd("./Data")

# 1996 Panel --------------------------------------------------------------

sipp96 <- readRDS("sipp96.RData")

# add soc2d codes & PCE
processed96 <- sipp96 %>%
        mutate(occDiv = as.integer(ifelse(occ >= 1000,  occ/10, occ))) %>%
        left_join(occ1990_to_SOC2d, by = c("occDiv" = "occ1990")) %>%
        select(-occDiv, -occ2000)

processed96 <- left_join(processed96, PCE, by = "date")

processed96 <- processWrapper(processed96)

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

processed01 <- left_join(processed01, PCE, by = "date")

processed01 <- processWrapper(processed01)

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

processed04 <- left_join(processed04, PCE, by = "date")

processed04<- processWrapper(processed04)

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

processed08 <- left_join(processed08, PCE, by = "date")

processed08 <- processWrapper(processed08)

if(useSoc2d) {
        setwd("./soc2d")
} else {
        setwd("./occ")
}
saveRDS(processed08, "processed08.RData")

rm(list = c("sipp08", "processed08"))

