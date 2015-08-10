# April 2, 2015
# Calculate wage percentiles, probability of switching by wage
# percentile, and mean income changes after switching
# Plot results, save to ./Figures/ directory.
# Precondition: calculateResidualWages.R has been run.
library(dplyr)
library(stats)
library(zoo)
library(reshape2)
library(ggplot2)

setwd("~/workspace/CVW/R")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T
useRegResid <- T

# Function for weighted moving average
# x, wts, and new.name are strings
movingAverage <- function(df, x, wts, new.name) {
        sumProduct <- rollsum(df[x] * df[wts], 13, align = "center", fill = NA)
        sumWts <- rollsum(df[wts], 13, align = "center", fill = NA)
        movingAvg <- sumProduct/sumWts
        df[new.name] <- movingAvg
        return(df)
}

# Create cumuluative distribution function
calculateCumul <- function(df) {
        cumulFunc <- ecdf(df$useWage)
        percentile <- round(cumulFunc(df$useWage)*100)
        data.frame(df, percentile)
}

# 1996 Panel - 2008 Panel --------------------------------------------------------------

setwd("./Data")

if(useSoc2d & useRegResid) {
	setwd("soc2d/RegResid")
} else if(useSoc2d & !useRegResid){
	setwd("soc2d/Raw")
} else if(!useSoc2d & useRegResid){
	setwd("occ/RegResid")
} else if(!useSoc2d & !useRegResid){
	setwd("occ/Raw")
}

analytic9608 <- readRDS("./analytic9608.RData")   

if(useSoc2d){
	analytic9608$soc2d<-analytic9608$occ
}


# Calculate income percentile within two-digit SOC code, append
incomePercentiles <- analytic9608 %>%
        filter(!is.na(soc2d) & lfStat==1) %>%
        group_by(soc2d) %>%
        do(calculateCumul(.)) %>%
        select(wpfinwgt, switchedOcc,switchedJob, percentile, EE, UE, recIndic, waveRec , date)

saveRDS(incomePercentiles,"incomePercentiles.RData")
rm(incomePercentiles)
	

# Remove from environment
rm(analytic9608)

setwd("../../../")

# Plots -------------------------------------------------------------------



if(useSoc2d) {
	if(useRegResid){
		incomePercentiles <- readRDS("./Data/soc2d/RegResid/incomePercentiles.RData")
	}else{
		incomePercentiles <- readRDS("./Data/soc2d/Raw/incomePercentiles.RData")
	}
	setwd("./Figures/soc2d")
} else {
	if(useRegResid){
		incomePercentiles <- readRDS("./Data/occ/RegResid/incomePercentiles.RData")
	}else{
		incomePercentiles <- readRDS("./Data/occ/Raw/incomePercentiles.RData")
	}
    setwd("./Figures/occ")
}


incomePercentiles <- incomePercentiles %>%
	filter(!is.na(percentile)) %>%
	group_by(percentile) %>%
	summarize(prSwitching = weighted.mean(switchedOcc & (UE | EE), wpfinwgt, na.rm = TRUE) /
			  	weighted.mean( (UE | EE), wpfinwgt, na.rm = TRUE),
			  prSwitchingEE = weighted.mean(switchedOcc & EE, wpfinwgt, na.rm = TRUE) /
			  	weighted.mean(EE, wpfinwgt, na.rm = TRUE),
			  prSwitchingUE = weighted.mean(switchedOcc & UE, wpfinwgt, na.rm = TRUE) /
			  	weighted.mean(UE, wpfinwgt, na.rm = TRUE),
			  prSwitchingRec = weighted.mean(switchedOcc & (UE|EE) & recIndic, wpfinwgt, na.rm = TRUE) /
			  	weighted.mean((UE|EE) & recIndic, wpfinwgt, na.rm = TRUE),
			  prSwitchingExp = weighted.mean(switchedOcc & (UE|EE) & !recIndic, wpfinwgt, na.rm = TRUE) /
			  	weighted.mean((UE|EE) & !recIndic, wpfinwgt, na.rm = TRUE),
			  prSwitchingUnc = weighted.mean(switchedOcc & (UE | EE) , wpfinwgt, na.rm = TRUE),
			  prSwitchingEEUnc = weighted.mean(switchedOcc & EE, wpfinwgt, na.rm = TRUE) ,
			  prSwitchingRecUnc = weighted.mean(switchedOcc & (UE|EE) & recIndic, wpfinwgt, na.rm = TRUE) ,
			  prSwitchingExpUnc = weighted.mean(switchedOcc & (UE|EE) & !recIndic, wpfinwgt, na.rm = TRUE)
	)

ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitching)) + geom_point() + geom_smooth() +
	ggtitle("P(Switching|Job change) by Earnings Percentile") + 
	ylab("P(Switching)") + xlab("Earnings percentile within two-digit SOC occupation")
ggsave("prSwitching.eps", width = 6, height = 4)

ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitchingEE)) + geom_point() + geom_smooth() +
	ggtitle("P(Switching| EE job change) by Wage Percentile") + 
	ylab("P(Switching)") + xlab("Wage percentile within two-digit SOC occupation")
ggsave("prSwitchingEE.eps", width = 6, height = 4)

ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitchingRec)) + geom_point() + geom_smooth() +
	ggtitle("P(Switching|Job change & Rec) by Earnings Percentile") + 
	ylab("P(Switching)") + xlab("Earnings percentile within two-digit SOC occupation")
ggsave("prSwitchingRec.eps", width = 6, height = 4)

ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitchingExp)) + geom_point() + geom_smooth() +
	ggtitle("P(Switching|Job change & Expansion) by Earnings Percentile") + 
	ylab("P(Switching)") + xlab("Earnings percentile within two-digit SOC occupation")
ggsave("prSwitchingExp.eps", width = 6, height = 4)



ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitchingUnc)) + geom_point() + geom_smooth() +
	ggtitle("Unconditional P(Switching) by Wage Percentile") + 
	ylab("P(Switching)") + xlab("Wage percentile within two-digit SOC occupation")
ggsave("prSwitchingByIncomeUnc.eps", width = 6, height = 4)

postscript("prSwitchingByIncomeRec.eps", width = 782, height = 569)
ggPr<-ggplot(incomePercentiles, aes(percentile, prSwitching, color = Rec)) + 
	geom_smooth() +
	ggtitle("P(Switching) by Income Percentile") + 
	ylab("P(Switching)") + xlab("Income percentile within occupation")
dev.off()

setwd("../../")

# Income changes by date -------------------------------
# Calculate mean residual income changes by labor force flow

if(useSoc2d) {
	if(useRegResid){
		setwd("./Data/soc2d/RegResid/")
	}else{
		setwd("./Data/soc2d/Raw/")
	}
} else {
	if(useRegResid){
		setwd("./Data/occ/RegResid/")
	}else{
		setwd("./Data/occ/Raw/")
	}
}
analytic9608 <- readRDS("analytic9608.RData")

incomeChanges <- analytic9608 %>%
	filter(switchedOcc) %>%
	group_by(date) %>%
	summarize(meanIncomeChange = weighted.mean(residWageChange[EE | UE],
											   wpfinwgt[EE | UE], na.rm = TRUE),
			  meanIncomeChangeEE = weighted.mean(residWageChange[EE], 
			  								   wpfinwgt[EE], na.rm = TRUE),
			  meanIncomeChangeUE = weighted.mean(residWageChange[UE], 
			  								   wpfinwgt[UE], na.rm = TRUE),
			  numPooled = sum(wpfinwgt[EE | UE], na.rm = TRUE),
			  numEE = sum(wpfinwgt[EE], na.rm = TRUE),
			  numUE = sum(wpfinwgt[UE], na.rm = TRUE)) 


incomeChanges <- incomeChanges %>%
	# why are there NaNs?
	filter(!is.nan(meanIncomeChangeEE)) %>%
	group_by(panel) %>%
	arrange(date) %>%
	do(movingAverage(., "meanIncomeChange", "numPooled", new.name = "maPooled")) %>%
	do(movingAverage(., "meanIncomeChangeEE", "numEE", new.name = "maEE")) %>%
	do(movingAverage(., "meanIncomeChangeUE", "numUE", new.name = "maUE")) %>%
	select(date, starts_with("ma"), panel)


saveRDS(incomeChanges,"incomeChanges.RData")

if(useSoc2d) {
	setwd("./Figures/soc2d")
} else {
	setwd("./Figures/occ")
}

incomeChanges <- melt(incomeChanges, id = c("date", "panel"))
levels(incomeChanges$variable) <- c("Pooled", "Employment-to-Employment", "Unemployment-to-Employment")

png("incomeChangesMA.png", width = 782, height = 569)
ggplot(incomeChanges, aes(date, value, color = panel)) + geom_point() + geom_line() +
        facet_grid(. ~ variable) +
        ggtitle("Moving Average 6(1)6 of Residual Log Income Change by LF Flow") +
        xlab("Date") + ylab("Mean Change in Residual Log Income") +
        labs(color = "Panel")
dev.off()

setwd("../../")