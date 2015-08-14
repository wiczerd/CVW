# April 17, 2015
# Calculate occupational change statistics:
# 1) mean EE, UE, and pooled rate of occupational change over the whole period
# 2) correlation between rates of occupational change and the unemployment rate
# 3) quantile regression of unemployment rate on switching rate
# For more statistics, see summarizeData.R
# Precondition: processData.R has been run.
library(weights)
library(dplyr)
library(erer)
library(stats)
library(zoo)
library(ggplot2)
library(xlsx)
library(stargazer)
library(quantreg)

setwd("~/workspace/CVW/R")

# Use 2 digit occupations from mapped from SOC? (SOC 2d)
useSoc2d <- T

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
                   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
        mutate(month = format(date, "%m"),
               year = format(date, "%Y"),
               date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        select(-year, -month)

# push forward and pull backward occupation
nextlastocc <- function(df){
	df <- df %>%
		group_by(id) %>%
		arrange(id,date)
	df$nextOcc = ifelse(df$switchedJob & lead(df$job) != 0, lead(df$occ), NA) 
	df$nextOcc = na.locf(df$nextOcc, na.rm = FALSE, fromLast=TRUE)
	df$lastOcc = ifelse(df$switchedJob & df$job != 0, lag(df$occ), NA)
	df$lastOcc = na.locf(df$lastOcc, na.rm = FALSE) 
	ungroup(df)
	return(df)
}

setwd("./Data")


if(useSoc2d) {
	setwd("./soc2d")
}else {
	setwd("./occ")
}

# Extract data for the switching probit ---------------------------

processed9608 <- readRDS("processed96.RData")
processed01  <- readRDS("processed01.RData")
processed9608<-bind_rows(processed01,processed9608)
rm(processed01)
processed04  <- readRDS("processed04.RData")
processed9608<-bind_rows(processed04,processed9608)
rm(processed04)
processed08  <- readRDS("processed08.RData")
processed9608<-bind_rows(processed08,processed9608)
rm(processed08)


# Extract data for the switching probit
demoKeepVars <- c("wpfinwgt","id","race","educ","switchedJob"
				  ,"switchedOcc","unempDur","date","lfStat"
				  ,"age","female","occ","UE","EE","job","recIndic","waveRec")
demoProbit <-  select(processed9608,one_of(demoKeepVars))
#demoProbit <-  nextlastocc(demoProbit)
#fix SwitchedOcc ?
#demoProbit$switchedOcc <- ifelse(demoProbit$UE, (demoProbit$occ != demoProbit$nextOcc), demoProbit$switchedOcc)
#demoProbit$switchedOcc <- ifelse(demoProbit$UE & (is.na(demoProbit$occ) | is.na(demoProbit$nextOcc)),
#								 NA, demoProbit$switchedOcc)
#demoProbit$switchedOcc <- ifelse(demoProbit$EE & (is.na(demoProbit$occ) ),
#								 NA, demoProbit$switchedOcc)
demoProbit <-  mutate(demoProbit,swOccJob = switchedOcc & (UE|EE) & !is.na(occ), na.rm =T)
# save demo probit
saveRDS(demoProbit, "demoProbit.RData")	

rm(list=c("demoProbit","processed9608"))

# Probit switching estimation ----------------------------------------------

demoProbit <- readRDS("./demoProbit.RData")

demoProbit <- subset(demoProbit,!is.na(occ) & (UE|EE)) %>%
	left_join(haver) %>%
	mutate(nwhite = as.integer(race > 1)) %>%
	mutate(univ = as.integer(educ >= 4)) %>%
	mutate(lths = as.integer(educ == 1)) %>%
	select(-educ,-race)

demoProbit <- mutate(demoProbit,unrateSA = unrateSA/100,
					 unrateNSA = unrateNSA/100)


swDemo.unrate.EE <- glm(swOccJob ~ unrateSA + age + I(age^2/100)  + nwhite + univ + lths + female, 
						family=binomial(link="probit"), subset=(EE==1), data= demoProbit, na.action=na.omit, x=T)
swDemo.unrate.UE <- glm(swOccJob ~ unrateSA + age + I(age^2/100) + unempDur + nwhite + univ + lths + female, 
						family=binomial(link="probit"), subset=(UE==1), data= demoProbit, na.action=na.omit, x=T)
swDemo.unrate <- glm(swOccJob ~ unrateSA  + age + I(age^2/100) +  UE + nwhite + univ + lths + female, 
					 family=binomial(link="probit"), data= demoProbit, na.action=na.omit, x=T)

swDemo.EE <- glm(swOccJob ~ recIndic + nwhite + univ + lths + female + age  + I(age^2), 
				 family=binomial(link="probit"), subset=(EE==1), data= demoProbit, na.action=na.omit)
swDemo.UE <- glm(swOccJob ~ recIndic + nwhite + univ + lths + female + age  + I(age^2) + unempDur, 
				 family=binomial(link="probit"), subset=(UE==1), data= demoProbit, na.action=na.omit)
swDemo <- glm(swOccJob ~ recIndic + nwhite + univ + lths + female + UE + age  + I(age^2) + unempDur, 
			  family=binomial(link="probit"), data= demoProbit, na.action=na.omit)

names(swDemo.unrate$coefficients)<-c("Const","Unemp Rate","Age","Age^2","Unemp Indic", "Non-white","Univ","LT HS","Female")
names(swDemo.unrate.EE$coefficients)<-c("Const","Unemp Rate","Age","Age^2", "Non-white","Univ","LT HS","Female")
names(swDemo.unrate.UE$coefficients)<-c("Const","Unemp Rate","Age","Age^2","Unemp Duration", "Non-white","Univ","LT HS","Female")

#swDemo.Mfx <- maBina(w=swDemo)
swDemo.unrate.Mfx <- maBina(w=swDemo.unrate)
swDemo.unrate.EE.Mfx <- maBina(w=swDemo.unrate.EE)
swDemo.unrate.UE.Mfx <- maBina(w=swDemo.unrate.UE)

setwd("../../")
if(useSoc2d) {
	setwd("./Figures/soc2d")
}else {
	setwd("./Figures/occ")
}

stargazer(swDemo.unrate.Mfx,swDemo.unrate.EE.Mfx,swDemo.unrate.UE.Mfx, 
		  title="Probit for occupational switching",out="swDemounrate.tex",
		  colnames=T,digits=3, column.sep.width= "1pt",df=F,no.space=T,
		  dep.var.labels.include = F)
#texreg( list(swDemo.unrate,swDemo.unrate,swDemo.unrate),label="tab:swDemo",caption="Probit for occupational switching",file="swDemo.tex")
setwd("../../")

if(useSoc2d) {
	setwd("./Data/soc2d")
}else {
	setwd("./Data/occ")
}


# Scroll through 1996-2008 Panels for switching TS---------------------


processed96 <- readRDS("processed96.RData")

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
processed01 <- readRDS("./processed01.RData")


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
processed04 <- readRDS("./processed04.RData")


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
processed08 <- readRDS("./processed08.RData")

  processed08 <- mutate(processed08, switchedOcc = as.logical(ifelse( is.na(occ),NA,switchedOcc ) ) )

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

# save prSwitching TS
saveRDS(prSwitching, "./prSwitching.RData")

# Statistics --------------------------------------------------------------

# Drop UE and pooled observations when the number of UE is too low 
# (ad-hoc: less than 10th percentile for panel)
prSwitchingRaw <- prSwitching
prSwitching <- prSwitching %>%
        group_by(panel) %>%
        mutate(cutoffUE = quantile(occUEObs, probs = .05, na.rm = TRUE),
               prSwitchedOccUE = ifelse(occUEObs < cutoffUE, 
                                        NA, prSwitchedOccUE),
               prSwitchedOcc = ifelse(occUEObs < cutoffUE, 
                                        NA, prSwitchedOcc)) %>%
        select(-cutoffUE)
        

# recession dates
rec_dates   <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01")) #as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2010-06-01")) 
prSwitching$Rec <- ((prSwitching$date>rec_dates[1] & prSwitching$date<rec_dates[2] ) | 
						(prSwitching$date>rec_dates[3] & prSwitching$date<rec_dates[4] ))

# Mean over whole period
with(prSwitching, weighted.mean(prSwitchedOcc, occObs, na.rm = TRUE))
with(prSwitching, weighted.mean(prSwitchedOccEE, occEEObs, na.rm = TRUE))
with(prSwitching, weighted.mean(prSwitchedOccUE, occUEObs, na.rm = TRUE))

# Correlation
prSwitchingAndUnemployment <- prSwitching %>%
        select(date, starts_with("prSwitched"),starts_with("occ"), panel) %>%
        left_join(haver)
with(prSwitchingAndUnemployment, wtd.cors(prSwitchedOcc, unrateSA, weight=occObs )) 
with(prSwitchingAndUnemployment, wtd.cors(prSwitchedOccEE, unrateSA, weight=occEEObs))
with(prSwitchingAndUnemployment, wtd.cors(prSwitchedOccUE, unrateSA, weight=occUEObs))

with(prSwitchingAndUnemployment, cor(prSwitchedOcc, unrateSA, , use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccEE, unrateSA, use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccUE, unrateSA, use = "complete.obs"))

with(prSwitchingAndUnemployment, cor(prSwitchedOcc, as.numeric(date), , use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccEE, as.numeric(date), use = "complete.obs"))
with(prSwitchingAndUnemployment, cor(prSwitchedOccUE, as.numeric(date), use = "complete.obs"))

with(prSwitching, cor(prSwitchedOcc, as.numeric(Rec), , use = "complete.obs"))
with(prSwitching, cor(prSwitchedOccEE, as.numeric(Rec), use = "complete.obs"))
with(prSwitching, cor(prSwitchedOccUE, as.numeric(Rec), use = "complete.obs"))


rm(list=c("prSwitchingAndUnemployment","prSwitching","haver"))
