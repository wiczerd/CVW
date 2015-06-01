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
useRegResid <- F

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
            "residWageChange_q", "residWageChange_q_wU","logWage","logEarnm")

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
qtl_delw <- c(0.1, 0.25, .5, .75, 0.9)
wageChangesEE <- subset(wageChanges, EE)
wageOLSEE.nSu <- lm(residWageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesEE)
wageRegEE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau =qtl_delw, weights= wpfinwgt, data = wageChangesEE)
wageOLSEE.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, weights= wpfinwgt, data = wageChangesEE)
wageRegEE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = qtl_delw, weights= wpfinwgt, data = wageChangesEE)
EEols.nSu<-summary(wageOLSEE.nSu)
EEols.wSu<-summary(wageOLSEE.wSu)
EEqr.nSu <-summary(wageRegEE.nSu)
EEqr.wSu <-summary(wageRegEE.wSu)
quantreg::latex(EEqr.nSu,file="./Figures/EEqr_nSu",transpose=T,digits=3)
quantreg::latex(EEqr.wSu,file="./Figures/EEqr_wSu",transpose=T,digits=3)


wageChangesUE <- subset(wageChanges, UE)
wageOLSUE.nSu <- lm(residWageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesUE)
wageRegUE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUE)
wageOLSUE.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, weights= wpfinwgt, data = wageChangesUE)
wageRegUE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUE)
UEols.nSu <-summary(wageOLSUE.nSu)
UEols.wSu <-summary(wageOLSUE.wSu)
UEqr.nSu <-summary(wageRegUE.nSu)
UEqr.wSu <-summary(wageRegUE.wSu)
quantreg::latex(UEqr.nSu,file="./Figures/UEqr_nSu",transpose=T,digits=3)
quantreg::latex(UEqr.wSu,file="./Figures/UEqr_wSu",transpose=T,digits=3)

wageChangesUEEE <- subset(wageChanges, UE | EE)
wageOLS.nSu <- lm(residWageChange ~ switchedOcc + unrateSA + UE, weights= wpfinwgt, data = wageChangesUEEE)
wageReg.nSu <- rq(residWageChange ~ switchedOcc + unrateSA + UE, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
wageOLS.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE, weights= wpfinwgt, data = wageChangesUEEE)
wageReg.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
ols.nSu <-summary(wageOLS.nSu)
ols.wSu <-summary(wageOLS.wSu)
qr.nSu <-summary(wageReg.nSu)
qr.wSu <-summary(wageReg.wSu)
quantreg::latex(qr.nSu,file="./Figures/qr_nSu",transpose=T,digits=3)
quantreg::latex(qr.wSu,file="./Figures/qr_wSu",transpose=T,digits=3)

png("./Figures/qr_nSu.png")
plot(qr.nSu, parm=c("switchedOccTRUE","unrateSA"),xlab="Quantile")
dev.off()

png("./Figures/qr_wSu.png")
plot(qr.wSu, parm=c("switchedOccTRUE","switchedOccTRUE:unrateSA","unrateSA","UETRUE"),xlab="Quantile")
dev.off()

rm(wageChangesUE)
rm(wageChangesEE)

coef_wSu <- cbind(t(wageRegEE.wSu$coefficients),t(wageRegUE.wSu$coefficients))
coef_wSu <- rbind(cbind(t(wageOLSEE.wSu$coefficients),t(wageOLSUE.wSu$coefficients)), coef_wSu )

coef_nSu <- cbind(t(wageRegEE.nSu$coefficients),t(wageRegUE.nSu$coefficients))
coef_nSu <- rbind(cbind(t(wageOLSEE.nSu$coefficients),t(wageOLSUE.nSu$coefficients)), coef_nSu )

qr_nSu <-data.frame(qtl_delw,t(wageRegEE.nSu$coefficients),t(wageRegUE.nSu$coefficients))
ggA <- ggplot( qr_nSu, aes(y = switchedOccTRUE, x = qtl_delw)) +
	geom_line(size = 2, colour = "#66CC99") + geom_hline(aes(yintercept=wageOLSEE.nSu$coefficients[2]), colour = "#66CC99" , size = 2)
ggA <- ggA + geom_line( aes(y = switchedOccTRUE.1, x = qtl_delw), colour = "#9999CC", size = 2) +
	geom_hline(aes(yintercept=wageOLSUE.nSu$coefficients[2]), size = 2, colour="#9999CC") +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of Occupation Switching"))
ggsave("./Figures/qtl_sw_nSu.png",width = 5, height = 5)


ggB <- ggplot( qr_nSu, aes(y = unrateSA, x = qtl_delw)) +
	geom_line(size = 2, colour = "#66CC99") + geom_hline(aes(yintercept=wageOLSEE.nSu$coefficients[3]), size = 2, colour = "#66CC99" )
ggB <- ggB + geom_line( aes(y = unrateSA.1, x = qtl_delw), colour = "#9999CC", size = 2) +
	geom_hline(aes(yintercept=wageOLSUE.nSu$coefficients[3]), colour="#9999CC", size = 2) +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of Unemployment Rate"))
ggsave("./Figures/qtl_u_nSu.png",width = 5, height = 5)

# splitting between young and old
