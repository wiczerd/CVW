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
library(McSpatial)
library(reshape2)

setwd("~/workspace/CVW/R")

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
                   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
        mutate(month = format(date, "%m"),
               year = format(date, "%Y"),
               date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
        select(-year, -month)

toKeep <- c("wpfinwgt", "switchedOcc", "EE", "UE", "switched2d", 
            "residWageChange", "residWageChange_wU", "lfStat", "date",
            "residWageChange_q", "residWageChange_q_wU")

detach("package:xlsx")
detach("package:xlsxjars")
detach("package:rJava")

# Load data --------------------------------------------------------------
analytic96 <- readRDS("./Data/analytic96.RData")
wageChanges <- analytic96 %>%
        select(one_of(toKeep))
rm(analytic96)

analytic01 <- readRDS("./Data/analytic01.RData")
wageChanges <- analytic01 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic01)

analytic04 <- readRDS("./Data/analytic04.RData")
wageChanges <- analytic04 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic04)

analytic08 <- readRDS("./Data/analytic08.RData")
wageChanges <- analytic08 %>%
        select(one_of(toKeep)) %>%
        bind_rows(wageChanges)
rm(analytic08)

#merge in unemployment
wageChanges <- left_join(wageChanges,haver, by="date")


# store full set
saveRDS(wageChanges, "./Data/wageChanges.RData")

# throw out all but job switches
wageChanges <- wageChanges %>%
        mutate(posChange = (residWageChange > 0),
               negChange = (residWageChange < 0)) %>%
	# drop the ones with no wage change (i.e missing values)
        filter(!(is.nan(residWageChange) & is.nan(residWageChange_wU) ) )

# Summary statistics --------------------------------------------------------------

# Mean wage changes
with(wageChanges, wtd.mean(residWageChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE))
with(wageChanges, wtd.mean(residWageChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE))

# Standard deviation of wage changes
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(residWageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], na.rm = TRUE)))

# Median of wage changes
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], probs = .5))
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], probs = .5))
with(wageChanges, wtd.quantile(residWageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], probs = .5))

# Fraction of workers with positive and negative wage changes
# explicitly calculate negative
with(wageChanges, wtd.mean(posChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(posChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE)) 
with(wageChanges, wtd.mean(posChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE))
with(wageChanges, wtd.mean(negChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(negChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE)) 
with(wageChanges, wtd.mean(negChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE)) 

# Correlation
dirWageChanges <- wageChanges %>%
        group_by(date) %>%
        summarize(pctPos = wtd.mean(posChange[switchedOcc & (EE | UE)], 
                                    wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE),
                  pctPosEE = wtd.mean(posChange[switchedOcc & EE], 
                                     wpfinwgt[switchedOcc & EE], na.rm = TRUE),
                  pctPosUE = wtd.mean(posChange[switchedOcc & UE], 
                                      wpfinwgt[switchedOcc & UE], na.rm = TRUE))

with(dirWageChanges, cor(pctPos, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosEE, unrateNSA, use = "complete.obs"))
with(dirWageChanges, cor(pctPosUE, unrateNSA, use = "complete.obs"))


# Quantile regressions ----------------------------------------------------
wageChanges <- readRDS( "./Data/wageChanges.RData")

wageChangesEE <- subset(wageChanges, EE)
wageRegEE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesEE)
EEqr.nSu <-summary(wageRegEE.nSu)
wageRegEE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesEE)
EEqr.wSu <-summary(wageRegEE.wSu)

wageChangesUE <- subset(wageChanges, UE)
wageRegUE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesUE)
UEqr.nSu <-summary(wageRegUE.nSu)
wageRegUE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesUE)
UEqr.wSu <-summary(wageRegUE.wSu)

rm(wageChangesUE)
rm(wageChangesEE)


# Wage change distributions ---------------------------------------------
wageChanges <- readRDS( "./Data/wageChanges.RData")

# recession dates
rec_dates   <- as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2009-06-01"))
wageChanges$Rec <- ((wageChanges$date>rec_dates[1] & wageChanges$date<rec_dates[2] ) | 
                            (wageChanges$date>rec_dates[3] & wageChanges$date<rec_dates[4] ))

wageChangesLong <- melt(wageChanges, id.vars = c("id", "date", "Rec"),
                        measure.vars = c("residWageChange_wU", "residWageChange_q_wU"))
wageChangesLong <- subset(wageChangesLong, !is.na(Rec))

# format variables for better plotting
wageChangesLong$Rec <- as.factor(wageChangesLong$Rec)
levels(wageChangesLong$Rec) <- c("Expansion", "Recession")
levels(wageChangesLong$variable) <- c("Monthly", "Quarterly")

# exp/rec panes
png("./Figures/wageChangeDensityRecExp.png", width = 782, height = 569)
ggplot(wageChangesLong, aes(value, fill = variable)) +
        geom_density(alpha = 0.5) +
        facet_grid(. ~ Rec) +
        xlim(c(-3.75, 3.75)) +
        ggtitle("Wage change distribution in expansion & recession") +
        labs(fill = "Aggregation")
dev.off()

# monthly/quarterly panes
png("./Figures/wageChangeDensityMthQtr.png", width = 782, height = 569)
ggplot(wageChangesLong, aes(value, fill = Rec)) +
        geom_density(alpha = 0.5) +
        facet_grid(. ~ variable) +
        xlim(c(-3.75, 3.75)) +
        ggtitle("Wage change distribution by aggregation") +
        labs(fill = "Business cycle")
dev.off()

# # subset for recessions and expansions
# wageChangesRec <- subset(wageChanges,wageChanges$Rec)
# wageChangesExp <- subset(wageChanges,!wageChanges$Rec)
# 
# wageChangesRecLong <- melt(wageChangesRec, id.vars = c("id", "date"),
#                          measure.vars = c("residWageChange_wU", "residWageChange_q_wU"))
# wageChangesExpLong <- melt(wageChangesExp, id.vars = c("id", "date"),
#                            measure.vars = c("residWageChange_wU", "residWageChange_q_wU"))

resChangeOutlier <- quantile(wageChanges$residWageChange_wU,probs=c(0.025,0.975),na.rm=T)
wageChanges$Out <- (wageChanges$residWageChange_wU<resChangeOutlier[1] |
						wageChanges$residWageChange_wU>resChangeOutlier[2] |
						is.na(wageChanges$residWageChange_wU) )
wageChangesIn <- subset(wageChanges,!Out)
wageChangesIn.kde <- density(wageChangesIn$residWageChange_wU,na.rm=T)
plot(wageChangesIn.kde)

wageChangesIn.rec <- qplot(residWageChange_wU,na.rm=T,data=wageChangesIn, colour = Rec, geom="density")
plot(wageChangesIn.rec)
wageChangesIn.kde <- density(wageChangesIn$residWageChange_wU,na.rm=T)
plot(wageChangesIn.kde)

bp_wU<-boxplot(residWageChange_wU~Rec,data=wageChangesIn,names=c("Expansion","Recession"))
title("Wage change distribution when changing jobs or into unemployment")

bp_wA<-boxplot(residWageChange_wA~Rec,data=wageChangesIn,names=c("Expansion","Recession"))
title("Wage change distribution")

# Machado - Mata Decomposition ----------------------------------------

# do the regressions I'll use
mm_rq.expansion <- rq(residWageChange_wU ~ switchedOcc,tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt,  data=subset(wageChanges,!wageChanges$Rec))

mm_rq.recession <- rq(residWageChange_wU ~ switchedOcc,tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt,  data=subset(wageChanges,wageChanges$Rec))

