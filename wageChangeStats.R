# April 17, 2015
# Calculate wage change statistics:
# 1) mean, standard deviation, median wage change for non-occupation-switch job changes
# 2) mean, standard deviation, median wage change for EE, UE, pooled occupation switches
# 3) quantile regressions
# 4) fraction of workers with positive and negative wage changes
# 5) correlation between unemployment rate and fraction of positive changes
library(Hmisc)
library(dplyr)
library(ggplot2)
library(xlsx)
library(quantreg)
library(reshape2)

setwd("~/workspace/CVW/R")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T
useRegResid <- T

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
                   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
        mutate(month = format(date, "%m"),
               year = format(date, "%Y"),
               date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        select(-year, -month)

toKeep <- c("wpfinwgt", "switchedOcc","soc2d", "EE", "UE",  
            "residWageChange", "residWageChange_wU", "residWageChange_stayer","lfStat", "date",
            "residWageChange_q", "residWageChange_q_wU","occWage","occWageChange","resid",
			"lastResidWage","lastOccWage")

detach("package:xlsx")
detach("package:xlsxjars")
detach("package:rJava")

# Load data --------------------------------------------------------------

if(useSoc2d & useRegResid) {
	analytic96 <- readRDS( "./Data/analytic96soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic96 <- readRDS("./Data/analytic96soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic96 <- readRDS("./Data/analytic96Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic96 <- readRDS("./Data/analytic96Raw.RData")   
} else{
	analytic96 <- readRDS("./Data/analytic96.RData")   
}


wageChanges <- analytic96 %>%
        select(one_of(toKeep))
rm(analytic96)


if(useSoc2d & useRegResid) {
	analytic01 <- readRDS( "./Data/analytic01soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic01 <- readRDS("./Data/analytic01soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic01 <- readRDS("./Data/analytic01Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic01 <- readRDS("./Data/analytic01Raw.RData")   
} else{
	analytic01 <- readRDS("./Data/analytic01.RData")   
}

wageChanges <- analytic01 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic01)


if(useSoc2d & useRegResid) {
	analytic04 <- readRDS( "./Data/analytic04soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic04 <- readRDS("./Data/analytic04soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic04 <- readRDS("./Data/analytic04Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic04 <- readRDS("./Data/analytic04Raw.RData")   
} else{
	analytic04 <- readRDS("./Data/analytic04.RData")   
}

wageChanges <- analytic04 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic04)


if(useSoc2d & useRegResid) {
	analytic08 <- readRDS("./Data/analytic08soc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	analytic08 <- readRDS("./Data/analytic08soc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	analytic08 <- readRDS("./Data/analytic08Resid.RData")   
} else if(!useSoc2d & !useRegResid){
	analytic08 <- readRDS("./Data/analytic08Raw.RData")   
} else{
	analytic08 <- readRDS("./Data/analytic08.RData")   
}

wageChanges <- analytic08 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic08)

#merge in unemployment
wageChanges <- left_join(wageChanges,haver, by="date")

# throw out infinity and missing values
wageChanges <- wageChanges %>%
	# drop the ones with no wage change (i.e missing values)
	filter(!is.infinite(residWageChange) & !is.na(residWageChange_wU) & !is.nan(residWageChange_wU) )


# store full set
if(useSoc2d & useRegResid) {
	saveRDS(wageChanges,"./Data/wageChangesSoc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	saveRDS(wageChanges,"./Data/wageChangesSoc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	saveRDS(wageChanges,"./Data/wageChangesResid.RData")   
} else if(!useSoc2d & !useRegResid){
	saveRDS(wageChanges,"./Data/wageChangesRaw.RData")   
} else{
	saveRDS(wageChanges,"./Data/wageChanges.RData")   
}

# Summary statistics --------------------------------------------------------------

if(useSoc2d & useRegResid) {
	wageChanges<-readRDS("./Data/wageChangesSoc2dResid.RData")
} else if(useSoc2d & !useRegResid){
	wageChanges<-readRDS("./Data/wageChangesSoc2dRaw.RData")   
} else if(!useSoc2d & useRegResid){
	wageChanges<-readRDS("./Data/wageChangesResid.RData")   
} else if(!useSoc2d & !useRegResid){
	wageChanges<-readRDS("./Data/wageChangesRaw.RData")   
} else{
	wageChanges<-readRDS("./Data/wageChanges.RData")   
}


wageChanges <- wageChanges %>%
	mutate(posChange = (residWageChange > 0)) %>%
	mutate(negChange = (residWageChange < 0))

wageChangesQrtile <-with(wageChanges, wtd.quantile(residWageChange[EE | UE], wpfinwgt[EE | UE], na.rm = TRUE,probs=c(.25,.5,.75,.9) ) )
with(wageChanges, wtd.mean(switchedOcc[(EE | UE) & residWageChange>wageChangesQrtile[3]], 
						   wpfinwgt[(EE | UE) & residWageChange>wageChangesQrtile[3]], na.rm = TRUE))
with(wageChanges, wtd.mean(switchedOcc[(EE | UE) & residWageChange>wageChangesQrtile[4]], 
						   wpfinwgt[(EE | UE) & residWageChange>wageChangesQrtile[4]], na.rm = TRUE))
with(wageChanges, wtd.mean(switchedOcc[(EE | UE) & residWageChange<wageChangesQrtile[1]], 
						   wpfinwgt[(EE | UE) & residWageChange<wageChangesQrtile[1]], na.rm = TRUE))



# Mean wage changes
with(wageChanges, wtd.mean(residWageChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE))
# no change
with(wageChanges, wtd.mean(residWageChange[!switchedOcc & (EE | UE)], 
						   wpfinwgt[!switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[!switchedOcc & EE], 
						   wpfinwgt[!switchedOcc & EE], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[!switchedOcc & UE], 
						   wpfinwgt[!switchedOcc & UE], na.rm = TRUE))
# tot
with(wageChanges, wtd.mean(residWageChange[ (EE | UE)], 
						   wpfinwgt[ (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[ EE], 
						   wpfinwgt[ EE], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[ UE], 
						   wpfinwgt[ UE], na.rm = TRUE))


# Standard deviation of wage changes
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], na.rm = TRUE)))

# Median of wage changes
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], probs = c(.25,.5,0.75 ) ))
#no switch
with(wageChanges, wtd.quantile(residWageChange[!switchedOcc & (EE | UE)], 
							   wpfinwgt[!switchedOcc & (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(residWageChange[!switchedOcc & EE], 
							   wpfinwgt[!switchedOcc & EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(residWageChange[!switchedOcc & UE], 
							   wpfinwgt[!switchedOcc & UE], probs = c(.25,.5,0.75 ) ))
#tot
with(wageChanges, wtd.quantile(residWageChange[ (EE | UE)], 
							   wpfinwgt[ (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(residWageChange[ EE], 
							   wpfinwgt[ EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(residWageChange[ UE], 
							   wpfinwgt[ UE], probs = c(.25,.5,0.75 ) ))

# Fraction of workers with positive and negative wage changes
# explicitly calculate negative
with(wageChanges, wtd.mean(posChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(posChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE)) 
with(wageChanges, wtd.mean(posChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE))
#no swtich
with(wageChanges, wtd.mean(posChange[!switchedOcc & (EE | UE)], 
						   wpfinwgt[!switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(posChange[!switchedOcc & EE], 
						   wpfinwgt[!switchedOcc & EE], na.rm = TRUE)) 
with(wageChanges, wtd.mean(posChange[!switchedOcc & UE], 
						   wpfinwgt[!switchedOcc & UE], na.rm = TRUE))
#tot
with(wageChanges, wtd.mean(posChange[(EE | UE)], 
						   wpfinwgt[(EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(posChange[ EE], 
						   wpfinwgt[ EE], na.rm = TRUE))
with(wageChanges, wtd.mean(posChange[ UE], 
						   wpfinwgt[ UE], na.rm = TRUE))


# Correlation
dirWageChanges <- wageChanges %>%
        group_by(date) %>%
        summarize(pctPos = wtd.mean(posChange[switchedOcc & (EE | UE)], 
                                    wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE),
                  pctPosEE = wtd.mean(posChange[switchedOcc & EE], 
                                     wpfinwgt[switchedOcc & EE], na.rm = TRUE),
                  pctPosUE = wtd.mean(posChange[switchedOcc & UE], 
                                      wpfinwgt[switchedOcc & UE], na.rm = TRUE),
                  unrateNSA = first(unrateNSA))

with(dirWageChanges, cor(pctPos, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosEE, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosUE, unrateNSA, use = "complete.obs"))
dirWageChanges <- wageChanges %>%
	group_by(date) %>%
	summarize(pctPos = wtd.mean(posChange[!switchedOcc & (EE | UE)], 
								wpfinwgt[!switchedOcc & (EE | UE)], na.rm = TRUE),
			  pctPosEE = wtd.mean(posChange[!switchedOcc & EE], 
			  					wpfinwgt[!switchedOcc & EE], na.rm = TRUE),
			  pctPosUE = wtd.mean(posChange[!switchedOcc & UE], 
			  					wpfinwgt[!switchedOcc & UE], na.rm = TRUE),
			  unrateNSA = first(unrateNSA))

with(dirWageChanges, cor(pctPos, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosEE, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosUE, unrateNSA, use = "complete.obs"))


