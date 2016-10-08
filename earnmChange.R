# October 6, 2016
library(foreign)
library(data.table)
library(zoo)
library(ggplot2)
library(Hmisc)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

#rootdir <- "G:/Research_Analyst/Eubanks/Occupation Switching/SIPP to go"
rootdir <- "~/workspace/CVW/R/"
datadir <- paste(rootdir, "/InputDataDE", sep = "")
outputdir <- paste(rootdir, "/Results", sep = "")
figuredir <- paste(rootdir, "/Figures", sep = "")
intermed_plots = F
setwd(rootdir)

sipp <- readRDS(paste0(outputdir,"/sipp_2.RData"))

# select relevant variables, exclude recessions
sipp <- sipp[!(recIndic), c("panel","year","month","id","earnm", "wpfinwgt"), with=FALSE]

# only keep those there for a full year (take away all recession years entirely?)
sipp[, monthsInYr := .N, by=list(id,year)]
sipp <- sipp[monthsInYr == 12,]

# collapse to yearly frequency (sum of earnm, mean of weights)
sipp <- sipp[, list(earnm = sum(earnm, na.rm=TRUE),
		    wpfinwgt = mean(wpfinwgt, na.rm=TRUE)), by=list(id,year)]

# get earnings quantile by year, remove 0's first
sipp <- sipp[earnm > 0,]
sipp[, earnmQ := round(rank(earnm)/.N,2), by=year]

# get year to year earnings change
sipp[, next.earnm := shift(earnm, 1, type="lead"), by=id]
sipp[, earnmChg := (next.earnm/earnm - 1)*100]

# get average earnm by earnm quantile, remove 0 earnings
df <- sipp[, .(earnmChg = wtd.quantile(earnmChg, weights=wpfinwgt, probs=c(.1,.5,.9), na.rm=T)), by=earnmQ]
df$qtl <- rep(c(.1,.5,.9),101)

# plot
ggplot(df[earnmQ>=.05 & earnmQ<=.95,], aes(earnmQ, earnmChg,color=qtl,group=qtl)) + 
	geom_point() + 
	geom_line() +
	ylim(c(-25,200)) +
	geom_hline(yintercept=0) +
	ggtitle("Average Annual Earnings Change by Percentile") +
	xlab("Earnm percentile in base year") +
	ylab("Average percent change")
ggsave(paste0(figuredir, "/earnmPctile_trunc2.png"))

# check that distribution is as expected
df <- sipp[, .(earnm = mean(earnm, na.rm=TRUE)), by=list(year,earnmQ)]
ggplot(df, aes(earnmQ, earnm, color=factor(year))) + 
	geom_line() +
	ggtitle("Earnings Distribution by Year") +
	xlab("Earnm percentile") +
	ylab("Average Earnm (within 2 digit percentile)") +
	scale_color_discrete(name="Year")
ggsave(paste0(figuredir, "/earnmDistByYr.png"))




########## There is only garbage below this line! ##########

# get average annual earnm
sipp[, avgEarnm := mean(earnm, na.rm=TRUE), by=list(id,year)]

# try collapsed version
sipp_a <- sipp[!is.na(earnm),]
sipp_a <- sipp_a[, list(avgEarnm=mean(earnm), 
			wpfinwgt=mean(wpfinwgt)), by=list(id, year)]
sipp_a[, earnmQuantile := round(rank(avgEarnm)/.N,2), by=year]
sipp_a[, obsNum := 1:.N, by=id]
sipp_a[, firstObs := (obsNum == 1)]
sipp_a[, lastObs := (obsNum == max(obsNum)), by=id]
sipp_a <- sipp_a[firstObs == TRUE | lastObs == TRUE,]
sipp_a[, final.year := shift(year, 1, type="lead"), by=id]
sipp_a[, final.avgEarnm := shift(avgEarnm, 1, type="lead"), by=id]
sipp_a[, elapsedYears := (final.year - year)]
sipp_a <- sipp_a[firstObs == TRUE,]
sipp_a[, pctRate := ((final.avgEarnm/avgEarnm)^(1/elapsedYears)-1)*100]

df <- sipp_a[, .(avgChg = weighted.mean(pctRate, wpfinwgt, na.rm=TRUE)), by=earnmQuantile]

ggplot(df, aes(earnmQuantile, avgChg)) + 
	geom_point() + 
	geom_line() + 
	# xlim(c(.02,1)) + 
	# ylim(c(-25,100)) +
	geom_hline(yintercept=0) +
	ggtitle("Annual Earnings Growth by Percentile") +
	xlab("Earnm percentile at t=0") +
	ylab("Average annual percentage change")
ggsave(paste0(figuredir, "/earnmPctile_trunc.png"))


# calculate earnings quantiles in each year
sipp[!is.na(earnm), earnmQuantile := round(rank(earnm)/.N,2), by=year]

# identify first and last observations of an id
sipp[, obsNum := 1:.N, by=id]
sipp[, firstObs := (obsNum == 1)]
sipp[, lastObs := (obsNum == max(obsNum)), by=id]
sipp <- sipp[firstObs == TRUE | lastObs == TRUE,]

# find elapsed number of years (fractional)
sipp[, final.year := shift(year, 1, type="lead"), by=id]
sipp[, final.month := shift(month, 1, type="lead"), by=id]
sipp[, final.earnm := shift(earnm, 1, type="lead"), by=id]
sipp[, elapsedYears := (final.year - year) + (final.month - month)/12]
sipp <- sipp[firstObs == TRUE,]

sipp[, c("obsNum", "firstObs", "lastObs") := NULL]

# find annual rate of change in earnm
sipp[, pctRateEarnm := (final.earnm/earnm)^(1/elapsedYears) - 1]

df <- sipp[, .(meanPctRate = mean(pctRateEarnm, na.rm=TRUE)), by=earnmQuantile]





