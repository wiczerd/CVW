# April 17, 2015
# Calculate wage change regressions:
# Requires that wageChangeStats has been run to create wageChanges.RData
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

detach("package:xlsx")
detach("package:xlsxjars")
detach("package:rJava")

## Load data ----------------------------------------------------

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

## Regs on diffs ---------------------------------
wageChangesEE <- subset(wageChanges, EE)
delwageOLSEE <- lm(residWageChange ~ switchedOcc + unrateSA , weights= wpfinwgt, data = wageChangesEE)
delwageQREE <- rq(residWageChange ~ switchedOcc + unrateSA, tau =qtl_delw, weights= wpfinwgt, data = wageChangesEE)
EEdelols<-summary(delwageOLSEE)
EEdelqr <-summary(delwageQREE)
quantreg::latex(EEdelqr,file="./Figures/EEqr_nSu",transpose=T,digits=3)
# with occupation dummies for occupation origin
delwageOLSEE.odum <- lm(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesEE)
delwageQREE.odum <- rq(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesEE)
EEdelols.odum<-summary(delwageOLSEE.odum)
EEdelqr.odum <-summary(delwageQREE.odum)
quantreg::latex(EEdelqr.odum,file="./Figures/EEqr_nSu_odum",transpose=T,digits=3)



wageChangesUE <- subset(wageChanges, UE)
delwageOLSUE <- lm(residWageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesUE)
delwageQRUE <- rq(residWageChange ~ switchedOcc + unrateSA, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUE)
UEdelols <-summary(delwageOLSUE)
UEdelqr <-summary(delwageQRUE)
quantreg::latex(UEdelqr,file="./Figures/UEqr_nSu",transpose=T,digits=3)
# with occupation dummies for occupation origin
delwageOLSUE.odum <- lm(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesUE)
delwageQRUE.odum <- rq(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesUE)
UEdelols.odum<-summary(delwageOLSUE.odum)
UEdelqr.odum <-summary(delwageQRUE.odum)
quantreg::latex(UEdelqr.odum,file="./Figures/UEqr_nSu_odum",transpose=T,digits=3)


wageChangesUEEE <- subset(wageChanges, UE | EE)
delwageOLS <- lm(residWageChange ~ switchedOcc + unrateSA + UE, weights= wpfinwgt, data = wageChangesUEEE)
delwageQR <- rq(residWageChange ~ switchedOcc + unrateSA + UE, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
delols <-summary(delwageOLS)
delqr <-summary(delwageQR)
quantreg::latex(delqr,file="./Figures/qr_nSu",transpose=T,digits=3)
# with occupation dummies for occupation origin
delwageOLS.odum <- lm(residWageChange ~ switchedOcc + unrateSA + UE + factor(soc2d), weights= wpfinwgt, data = wageChangesUEEE)
delwageQR.odum <- rq(residWageChange ~ switchedOcc + unrateSA + UE + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
delols.odum<-summary(delwageOLS.odum)
delqr.odum <-summary(delwageQR.odum)
quantreg::latex(delqr.odum,file="./Figures/qr_nSu_odum",transpose=T,digits=3)


# Regs on lags -----------------------------------
## run regs with lags rather than just differences


wageChangesUEEE <- subset(wageChanges, (UE | EE))
wageOLS <- lm(residWageChange ~ unrateSA + UE + switchedOcc + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE)
wageReg <- rq(residWageChange ~ unrateSA + UE + switchedOcc + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)


# with only-changers
wageChangesUEEE_sw <- subset(wageChanges, (UE | EE) & switchedOcc)
wageOLS.sw <- lm(residWageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE_sw)
wageReg.sw <- rq(residWageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_sw)
ols.sw <-summary(wageOLS.sw)
qr.sw <-summary(wageReg.sw)

# only non-changer
wageChangesUEEE_nsw <- subset(wageChanges, (UE | EE) & !switchedOcc)
wageOLS.nsw <- lm(residWageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE_nsw)
wageReg.nsw <- rq(residWageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_nsw)

png("./Figures/qr_nSu.png")
plot(qr, parm=c("switchedOccTRUE","unrateSA"),xlab="Quantile")
dev.off()

rm(wageChangesUE)
rm(wageChangesEE)

coef_wageReg <- cbind(t(wageRegEE$coefficients),t(wageRegUE$coefficients))
coef_wageReg <- rbind(cbind(t(wageOLSEE$coefficients),t(wageOLSUE$coefficients)), coef_wageReg )

qr_wageReg <-data.frame(qtl_delw,t(wageRegEE$coefficients),t(wageRegUE$coefficients))
ggA <- ggplot( qr_wageReg, aes(y = switchedOccTRUE, x = qtl_delw)) +
	geom_line(size = 2, colour = "#66CC99") + geom_hline(aes(yintercept=wageOLSEE$coefficients[2]), colour = "#66CC99" , size = 2)
ggA <- ggA + geom_line( aes(y = switchedOccTRUE.1, x = qtl_delw), colour = "#9999CC", size = 2) +
	geom_hline(aes(yintercept=wageOLSUE$coefficients[2]), size = 2, colour="#9999CC") +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of Occupation Switching"))
ggsave("./Figures/qtl_sw_nSu.png",width = 5, height = 5)


ggB <- ggplot( qr_wageReg, aes(y = unrateSA, x = qtl_delw)) +
	geom_line(size = 2, colour = "#66CC99") + geom_hline(aes(yintercept=wageOLSEE$coefficients[3]), size = 2, colour = "#66CC99" )
ggB <- ggB + geom_line( aes(y = unrateSA.1, x = qtl_delw), colour = "#9999CC", size = 2) +
	geom_hline(aes(yintercept=wageOLSUE$coefficients[3]), colour="#9999CC", size = 2) +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of Unemployment Rate"))
ggsave("./Figures/qtl_u_nSu.png",width = 5, height = 5)

