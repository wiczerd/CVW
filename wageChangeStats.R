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

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- F
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

toKeep <- c("wpfinwgt", "switchedOcc", "EE", "UE",  
            "residWageChange", "residWageChange_wU", "lfStat", "date",
            "residWageChange_q", "residWageChange_q_wU")

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

# throw out all but job switches
wageChanges <- wageChanges %>%
        mutate(posChange = (residWageChange > 0)) %>%
        mutate(negChange = (residWageChange < 0))

# Summary statistics --------------------------------------------------------------


wageChangesQrtile <-with(wageChanges, wtd.quantile(residWageChange[EE | UE], wpfinwgt[EE | UE], na.rm = TRUE,probs=c(.25,.5,.75) ) )
with(wageChanges, wtd.mean(switchedOcc[(EE | UE) & residWageChange>wageChangesQrtile[3]], 
						   wpfinwgt[(EE | UE) & residWageChange>wageChangesQrtile[3]], na.rm = TRUE))

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


# Quantile regressions ----------------------------------------------------

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

wageChangesEE <- subset(wageChanges, EE)
wageOLSEE.nSu <- lm(residWageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesEE)
wageRegEE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesEE)
wageOLSEE.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, weights= wpfinwgt, data = wageChangesEE)
wageRegEE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesEE)
EEols.nSu<-summary(wageOLSEE.nSu)
EEols.wSu<-summary(wageOLSEE.wSu)
EEqr.nSu <-summary(wageRegEE.nSu)
EEqr.wSu <-summary(wageRegEE.wSu)

wageChangesUE <- subset(wageChanges, UE)
wageOLSUE.nSu <- lm(residWageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesUE)
wageRegUE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesUE)
wageOLSUE.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, weights= wpfinwgt, data = wageChangesUE)
wageRegUE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = c(0.1, 0.25, .5, .75, 0.9), weights= wpfinwgt, data = wageChangesUE)
UEols.nSu <-summary(wageOLSUE.nSu)
UEols.wSu <-summary(wageOLSUE.wSu)
UEqr.nSu <-summary(wageRegUE.nSu)
UEqr.wSu <-summary(wageRegUE.wSu)

rm(wageChangesUE)
rm(wageChangesEE)


# Wage change distributions ---------------------------------------------

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

# recession dates
rec_dates   <- as.Date(c("2001-03-01", "2001-11-01","2007-12-01", "2009-06-01"))
wageChanges$Rec <- ((wageChanges$date>rec_dates[1] & wageChanges$date<rec_dates[2] ) | 
                            (wageChanges$date>rec_dates[3] & wageChanges$date<rec_dates[4] ))
# exclude outliers
resChangeOutlier_rec <- quantile(wageChanges$residWageChange_wU[wageChanges$Rec == 1],probs=c(0.025,0.975),na.rm=T)
resChangeOutlier_exp <- quantile(wageChanges$residWageChange_wU[wageChanges$Rec == 0],probs=c(0.025,0.975),na.rm=T)
#wageChanges$Out <- (wageChanges$residWageChange_wU<resChangeOutlier[1] |
#						wageChanges$residWageChange_wU>resChangeOutlier[2] |
#						is.na(wageChanges$residWageChange_wU) )
wageChanges$Out <- (  (wageChanges$residWageChange_wU<resChangeOutlier_rec[1] |
					   	wageChanges$residWageChange_wU>resChangeOutlier_rec[2]) 
					  & wageChanges$Rec == 1 ) ||
	(  (wageChanges$residWageChange_wU<resChangeOutlier_exp[1] |
			wageChanges$residWageChange_wU>resChangeOutlier_exp[2]) 
			& wageChanges$Rec == 0 )
resChangeOutlier_q_rec <- quantile(wageChanges$residWageChange_q_wU[wageChanges$Rec == 1],probs=c(0.025,0.975),na.rm=T)
resChangeOutlier_q_exp <- quantile(wageChanges$residWageChange_q_wU[wageChanges$Rec == 0],probs=c(0.025,0.975),na.rm=T)
wageChanges$Out_q <- (  (wageChanges$residWageChange_q_wU<resChangeOutlier_q_rec[1] |
					   	wageChanges$residWageChange_q_wU>resChangeOutlier_q_rec[2]) 
					  & wageChanges$Rec == 1 ) ||
	(  (wageChanges$residWageChange_q_wU<resChangeOutlier_q_exp[1] |
			wageChanges$residWageChange_q_wU>resChangeOutlier_q_exp[2]) 
			& wageChanges$Rec == 0 ) 


wageChangesLong <- melt(wageChanges, id.vars = c("id", "date", "Rec"),
                        measure.vars = c("residWageChange_wU", "residWageChange_q_wU"))
wageChangesLong <- subset(wageChangesLong, !is.na(Rec))

# format variables for better plotting
wageChangesLong$Rec <- as.factor(wageChangesLong$Rec)
levels(wageChangesLong$Rec) <- c("Expansion", "Recession")
levels(wageChangesLong$variable) <- c("Monthly", "Quarterly")

# exp/rec panes
#postscript(file = "./Figures/wageChangeDensityRecExp.eps", width = 782, height = 569)
#ggplot(wageChangesLong, aes(value, fill = variable)) +
#        geom_density(alpha = 0.5) +
#        facet_grid(. ~ Rec) +
#        xlim(c(-3.75, 3.75)) +
#        ggtitle("Wage change distribution in expansion & recession") +
#        labs(fill = "Aggregation")
#dev.off()

# monthly/quarterly panes
png(file="./Figures/wageChangeDensityMthQtr.png", width = 782, height = 569)
ggplot(wageChangesLong, aes(value, fill = Rec)) +
        geom_density(alpha = 0.5) +
        facet_grid(. ~ variable) +
        xlim(c(-3.75, 3.75)) +
        ggtitle("Wage change distribution by aggregation") +
        labs(fill = "Business cycle")
dev.off()

# quarterly pane
png(file="./Figures/wageChangeDensityQtr.png",onefile=FALSE, width = 782, height = 569)
ggplot(subset(wageChangesLong,wageChangesLong$variable=="Quarterly"), aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	xlim(c(-3.75, 3.75)) +
	ggtitle("Wage change distribution by aggregation") +
	labs(fill = "Business cycle")
dev.off()

# monthly pane
png(file="./Figures/wageChangeDensityMth.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Monthly"), aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	xlim(c(-3.75, 3.75)) +
	ggtitle("Wage change distribution by aggregation") +
	labs(fill = "Business cycle")
dev.off()

#box plot of the same thing
png(file="./Figures/wageChangeBoxQtr.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Quarterly"), aes(y = value, x = Rec)) +
  geom_boxplot() +
  guides(fill=F) +
  stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
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

bp<-boxplot(residWageChange~Rec,data=wageChangesIn,names=c("Expansion","Recession"))
title("Wage change distribution, excluding unemployment")

# Machado - Mata Decomposition ----------------------------------------

# do the regressions I'll use
wageChanges <- subset(wageChanges, !wageChanges$Out)
wageChangesRec <- subset(wageChanges,wageChanges$Rec & (UE | EE))
wageChangesExp <- subset(wageChanges,!wageChanges$Rec  & (UE | EE))

# these are the quantile regressions we will be decomposing:
qtl_delw <- seq(0.2, .8, .1)
mm_rq.expansion <- rq(residWageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesExp)
mm_rq.recession <- rq(residWageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesRec )
mm_lm.expansion <- lm(residWageChange ~ switchedOcc +UE,weights= wpfinwgt,  data=wageChangesExp)
mm_lm.recession <- lm(residWageChange ~ switchedOcc +UE,weights= wpfinwgt,  data=wageChangesRec)


mm_bmat.exp <-qregbmat(residWageChange ~ switchedOcc + UE , taumat=qtl_delw, data=wageChangesExp, graphb = F)
mm_bmat.rec <-qregbmat(residWageChange ~ switchedOcc + UE , taumat=qtl_delw, data=wageChangesRec, graphb = F)

#use weighted regression? i.e:
#mm_bmat.exp <- t(mm_rq.expansion$coefficients)
#mm_bmat.rec <- t(mm_rq.recession$coefficients)
for(ki in 1:ncol(mm_bmat.exp)){
	for(ti in 1:nrow(mm_bmat.exp)){
		mm_bmat.exp[ti,ki] = mm_rq.expansion$coefficients[ki,ti]
		mm_bmat.rec[ti,ki] = mm_rq.recession$coefficients[ki,ti]
	}
}
qr_coefs <-data.frame(qtl_delw,mm_bmat.exp,mm_bmat.rec)
gg <- ggplot( qr_coefs, aes(y = switchedOccTRUE, x = qtl_delw)) +
	geom_line() + geom_point() 
gg <- gg + geom_line( aes(y = switchedOccTRUE.1, x = qtl_delw)) +
	ggtitle("Marginal Effect of Occupation Switching")
ggsave("./Figures/qtl_swocc.png",width = 5, height = 5)

gg <- ggplot( qr_coefs, aes(y = UETRUE, x = qtl_delw)) +
	geom_line() + geom_point() 
gg <- gg + geom_line( aes(y = UETRUE.1, x = qtl_delw)) +
	ggtitle("Marginal Effect of EUE")
ggsave("./Figures/qtl_eue.png",width = 5, height = 5)


mm_wageChanges <- qregsim2(residWageChange ~ switchedOcc + UE , ~ switchedOcc + UE, wageChangesExp, wageChangesRec,
						   mm_bmat.exp,mm_bmat.rec, timenames=c("Expansion","Recession"), nsim=5000)

