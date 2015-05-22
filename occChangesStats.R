# April 17, 2015
# Calculate occupational change statistics:
# 1) mean EE, UE, and pooled rate of occupational change over the whole period
# 2) correlation between rates of occupational change and the unemployment rate
# 3) quantile regression of unemployment rate on switching rate
# For more statistics, see summarizeData.R
# Precondition: processData.R has been run.
library(dplyr)
library(stats)
library(ggplot2)
library(xlsx)
library(quantreg)

setwd("~/workspace/CVW/R")

# Use 2 digit occupations from mapped from SOC? (SOC 2d)
useSoc2d <- F

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
                   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
        mutate(month = format(date, "%m"),
               year = format(date, "%Y"),
               date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        select(-year, -month)

# 1996 Panel --------------------------------------------------------------
if(useSoc2d) {
        processed96 <- readRDS("./Data/processed96soc2d.RData")
} else {
        processed96 <- readRDS("./Data/processed96.RData")
}

# Extract data for the switching probit
demoKeepVars <- c("wpfinwgt","race","educ","switchedJob"
				  ,"switchedOcc","unempDur","date","lfStat"
				  ,"age","female","occ","UE","EE")
demoProbit <-  select(processed96,one_of(demoKeepVars))
demoProbit <-  mutate(demoProbit,swOccJob = switchedOcc & switchedJob & !is.na(occ), na.rm =T)


# Calculate probability of switching using raw code
prSwitching <-  group_by(processed96, date) %>%
        summarize(prSwitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE),
                  occObs = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  occEEObs = sum(wpfinwgt[EE], na.rm = TRUE),
                  occUEObs = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = 1996)

# Remove 1996 data from environment
rm(processed96)

# 2001 Panel --------------------------------------------------------------
if(useSoc2d) {
        processed01 <- readRDS("./Data/processed01soc2d.RData")
} else {
        processed01 <- readRDS("./Data/processed01.RData")
}

demoProbit01 <-  select(processed01,one_of(demoKeepVars))
demoProbit01 <-  mutate(demoProbit01,swOccJob = switchedOcc & switchedJob & !is.na(occ), na.rm =T) %>%
	bind_rows(demoProbit)
demoProbit <- demoProbit01
rm(demoProbit01)

# Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed01, date) %>%
        summarize(prSwitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE),
                  occObs = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  occEEObs = sum(wpfinwgt[EE], na.rm = TRUE),
                  occUEObs = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = 2001) %>%
        bind_rows(prSwitching)

# Remove 2001 data from environment
rm(processed01)

# 2004 Panel --------------------------------------------------------------
if(useSoc2d) {
        processed04 <- readRDS("./Data/processed04soc2d.RData")
} else {
        processed04 <- readRDS("./Data/processed04.RData")
}

demoProbit04 <-  select(processed04,one_of(demoKeepVars))
demoProbit04 <-  mutate(demoProbit04,swOccJob = switchedOcc & switchedJob & !is.na(occ), na.rm =T) %>%
	bind_rows(demoProbit)
demoProbit <- demoProbit04
rm(demoProbit04)

#Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed04, date) %>%
        summarize(prSwitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE),
                  occObs = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  occEEObs = sum(wpfinwgt[EE], na.rm = TRUE),
                  occUEObs = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = 2004) %>%
        bind_rows(prSwitching)

# Remove 2004 data from environment
rm(processed04)

# 2008 Panel --------------------------------------------------------------
if(useSoc2d) {
  processed08 <- readRDS("./Data/processed08soc2d.RData")
} else {
  processed08 <- readRDS("./Data/processed08.RData")
}
processed08 <- mutate(processed08, switchedOcc = as.logical(ifelse( is.na(occ),NA,switchedOcc ) ) )

demoProbit08 <-  select(processed08,one_of(demoKeepVars))
demoProbit08 <-  mutate(demoProbit08,swOccJob = switchedOcc & switchedJob & !is.na(occ), na.rm =T) %>%
	bind_rows(demoProbit)
demoProbit <- demoProbit08
rm(demoProbit08)
# save demo probit
saveRDS(demoProbit, "./Data/demoProbit.RData")
rm(demoProbit)

#Calculate probability of switching, add to prSwitching
prSwitching <-  group_by(processed08, date) %>%
        summarize(prSwitchedOcc = weighted.mean(switchedOcc[EE | UE], wpfinwgt[EE | UE],
                                                na.rm = TRUE),
                  prSwitchedOccEE = weighted.mean(switchedOcc[EE], wpfinwgt[EE], na.rm = TRUE),
                  prSwitchedOccUE = weighted.mean(switchedOcc[UE], wpfinwgt[UE], na.rm = TRUE),
                  occObs = sum(wpfinwgt[EE | UE], na.rm = TRUE),
                  occEEObs = sum(wpfinwgt[EE], na.rm = TRUE),
                  occUEObs = sum(wpfinwgt[UE], na.rm = TRUE)) %>%
        mutate(panel = 2008) %>%
        bind_rows(prSwitching)

# Remove 2008 data from environment
rm(processed08)

# Statistics --------------------------------------------------------------

# Drop UE and pooled observations when the number of UE is too low 
# (ad-hoc: less than 10th percentile for panel)
prSwitchingRaw <- prSwitching
prSwitching <- prSwitching %>%
        group_by(panel) %>%
        mutate(cutoffUE = quantile(prSwitchedOccUE, probs = .15, na.rm = TRUE),
               prSwitchedOccUE = ifelse(prSwitchedOccUE < cutoffUE, 
                                        NA, prSwitchedOccUE),
               prSwitchedOcc = ifelse(prSwitchedOccUE < cutoffUE, 
                                        NA, prSwitchedOcc)) %>%
        select(-cutoffUE)
        

# Mean over whole period
with(prSwitching, weighted.mean(prSwitchedOcc, occObs, na.rm = TRUE))
with(prSwitching, weighted.mean(prSwitchedOccEE, occEEObs, na.rm = TRUE))
with(prSwitching, weighted.mean(prSwitchedOccUE, occUEObs, na.rm = TRUE))

# Correlation
prSwitchingAndUnemployment <- prSwitching %>%
        select(date, starts_with("prSwitched"), panel) %>%
        left_join(haver)
with(prSwitchingAndUnemployment, cor(prSwitchedOcc, unrateNSA, use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccEE, unrateNSA, use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccUE, unrateNSA, use = "complete.obs"))

# Quantile regressions ----------------------------------------------------

# Pooled
pooledReg <- rq(unrateNSA ~ prSwitchedOcc, tau = c(0.25, 0.50, 0.75), 
                data = prSwitchingAndUnemployment)

# EE
EEReg <- rq(unrateNSA ~ prSwitchedOccEE, tau = c(0.25, 0.50, 0.75), 
            data = prSwitchingAndUnemployment)

# UE
UEReg <- rq(unrateNSA ~ prSwitchedOccUE, tau = c(0.25, 0.50, 0.75), 
            data = prSwitchingAndUnemployment)

# All together
together <- rq(unrateNSA ~ prSwitchedOcc + prSwitchedOccEE + prSwitchedOccUE,
                 tau = c(0.25, 0.50, 0.75), data = prSwitchingAndUnemployment)

rm(list=c("prSwitchingAndUnemployment","prSwitching"))

# Probit switching ----------------------------------------------------------

demoProbit <- readRDS("./Data/demoProbit.RData")
demoProbit <- subset(demoProbit,!is.na(occ) & switchedJob) %>%
	left_join(haver) %>%
	mutate(nwhite = as.integer(race > 1)) %>%
	mutate(univ = as.integer(educ >= 4)) %>%
	mutate(lths = as.integer(educ == 1)) %>%
	select(-educ,-race)

rec_dates   <- as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2009-06-01"))
demoProbit$Rec <- ((demoProbit$date>rec_dates[1] & demoProbit$date<rec_dates[2] ) | 
						(demoProbit$date>rec_dates[3] & demoProbit$date<rec_dates[4] ))


swDemo.EE <- glm(swOccJob ~ unrateSA + nwhite + univ + lths + female + age , 
			  family=binomial(link="probit"), subset=(EE==1), data= demoProbit, na.action=na.omit)
#swDemo.UE <- glm(swOccJob ~ unrateSA + nwhite + univ + lths + female + age , 
#				 family=binomial(link="probit"), subset=(UE==1), data= demoProbit, na.action=na.omit)
#swDemo <- glm(swOccJob ~ unrateSA + nwhite + univ + lths + female + UE + age , 
#				 family=binomial(link="probit"), data= demoProbit, na.action=na.omit)

swDemo.EE <- glm(swOccJob ~ Rec + nwhite + univ + lths + female + age , 
				 family=binomial(link="probit"), subset=(EE==1), data= demoProbit, na.action=na.omit)
swDemo.UE <- glm(swOccJob ~ Rec + nwhite + univ + lths + female + age , 
				 family=binomial(link="probit"), subset=(UE==1), data= demoProbit, na.action=na.omit)
swDemo <- glm(swOccJob ~ Rec + nwhite + univ + lths + female + UE + age , 
			  family=binomial(link="probit"), data= demoProbit, na.action=na.omit)


rm(list=c("demoProbit","haver"))