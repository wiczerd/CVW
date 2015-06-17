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

toKeep <- c("wpfinwgt", "switchedOcc","soc2d", "EE", "UE",  
			"residWageChange", "residWageChange_wU", "residWageChange_stayer","lfStat", "date",
			"residWageChange_q", "residWageChange_q_wU","resid")

detach("package:xlsx")
detach("package:xlsxjars")
detach("package:rJava")

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
wageOLSEE.nSu <- lm(residWageChange ~ switchedOcc + unrateSA , weights= wpfinwgt, data = wageChangesEE)
wageRegEE.nSu <- rq(residWageChange ~ switchedOcc + unrateSA, tau =qtl_delw, weights= wpfinwgt, data = wageChangesEE)
wageOLSEE.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, weights= wpfinwgt, data = wageChangesEE)
wageRegEE.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA, tau = qtl_delw, weights= wpfinwgt, data = wageChangesEE)
EEols.nSu<-summary(wageOLSEE.nSu)
EEols.wSu<-summary(wageOLSEE.wSu)
EEqr.nSu <-summary(wageRegEE.nSu)
EEqr.wSu <-summary(wageRegEE.wSu)
quantreg::latex(EEqr.nSu,file="./Figures/EEqr_nSu",transpose=T,digits=3)
quantreg::latex(EEqr.wSu,file="./Figures/EEqr_wSu",transpose=T,digits=3)
# with occupation dummies
wageOLSEE.nSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesEE)
wageRegEE.nSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesEE)
wageOLSEE.wSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesEE)
wageRegEE.wSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + factor(soc2d), tau = qtl_delw, weights= wpfinwgt, data = wageChangesEE)
EEols.nSu.odum<-summary(wageOLSEE.nSu.odum)
EEols.wSu.odum<-summary(wageOLSEE.wSu.odum)
EEqr.nSu.odum <-summary(wageRegEE.nSu.odum)
EEqr.wSu.odum <-summary(wageRegEE.wSu.odum)
quantreg::latex(EEqr.nSu.odum,file="./Figures/EEqr_nSu_odum",transpose=T,digits=3)
quantreg::latex(EEqr.wSu.odum,file="./Figures/EEqr_wSu_odum",transpose=T,digits=3)



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
# with occupation dummies
wageOLSUE.nSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesUE)
wageRegUE.nSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesUE)
wageOLSUE.wSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + factor(soc2d), weights= wpfinwgt, data = wageChangesUE)
wageRegUE.wSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + factor(soc2d), tau = qtl_delw, weights= wpfinwgt, data = wageChangesUE)
UEols.nSu.odum<-summary(wageOLSUE.nSu.odum)
UEols.wSu.odum<-summary(wageOLSUE.wSu.odum)
UEqr.nSu.odum <-summary(wageRegUE.nSu.odum)
UEqr.wSu.odum <-summary(wageRegUE.wSu.odum)
quantreg::latex(UEqr.nSu.odum,file="./Figures/UEqr_nSu_odum",transpose=T,digits=3)
quantreg::latex(UEqr.wSu.odum,file="./Figures/UEqr_wSu_odum",transpose=T,digits=3)


wageChangesUEEE <- subset(wageChanges, UE | EE)
wageOLS.nSu <- lm(residWageChange ~ switchedOcc + unrateSA + UE + resid, weights= wpfinwgt, data = wageChangesUEEE)
wageReg.nSu <- rq(residWageChange ~ switchedOcc + unrateSA + UE + resid, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
wageOLS.wSu <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE  + resid, weights= wpfinwgt, data = wageChangesUEEE)
wageReg.wSu <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE + resid, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
ols.nSu <-summary(wageOLS.nSu)
ols.wSu <-summary(wageOLS.wSu)
qr.nSu <-summary(wageReg.nSu)
qr.wSu <-summary(wageReg.wSu)
quantreg::latex(qr.nSu,file="./Figures/qr_nSu",transpose=T,digits=3)
quantreg::latex(qr.wSu,file="./Figures/qr_wSu",transpose=T,digits=3)
# with occupation dummies
wageOLS.nSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + UE + factor(soc2d), weights= wpfinwgt, data = wageChangesUEEE)
wageReg.nSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + UE + factor(soc2d), tau =qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
wageOLS.wSu.odum <- lm(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE + factor(soc2d), weights= wpfinwgt, data = wageChangesUEEE)
wageReg.wSu.odum <- rq(residWageChange ~ switchedOcc + unrateSA + switchedOcc*unrateSA + UE + factor(soc2d), tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)
ols.nSu.odum<-summary(wageOLS.nSu.odum)
ols.wSu.odum<-summary(wageOLS.wSu.odum)
qr.nSu.odum <-summary(wageReg.nSu.odum)
qr.wSu.odum <-summary(wageReg.wSu.odum)
quantreg::latex(qr.nSu.odum,file="./Figures/qr_nSu_odum",transpose=T,digits=3)
quantreg::latex(qr.wSu.odum,file="./Figures/qr_wSu_odum",transpose=T,digits=3)


# with only-changers
wageChangesUEEE_sw <- subset(wageChanges, (UE | EE) & switchedOcc)
wageOLS.nSu.sw <- lm(residWageChange ~ unrateSA + UE + resid, weights= wpfinwgt, data = wageChangesUEEE_all)
wageReg.nSu.sw <- rq(residWageChange ~ unrateSA + UE + resid, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_all)
wageOLS.wSu.sw <- lm(residWageChange ~ unrateSA + switchedOcc*unrateSA + UE  + resid, weights= wpfinwgt, data = wageChangesUEEE_all)
wageReg.wSu.all <- rq(residWageChange ~ unrateSA + switchedOcc*unrateSA + UE + resid, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_all)
ols.nSu.all <-summary(wageOLS.nSu.all)
ols.wSu.all <-summary(wageOLS.wSu.all)
qr.nSu.all <-summary(wageReg.nSu.all)
qr.wSu.all <-summary(wageReg.wSu.all)

# only non-changer
wageChangesUEEE_nsw <- subset(wageChanges, (UE | EE) & !switchedOcc)
wageOLS.nSu.nsw <- lm(residWageChange ~ unrateSA + UE + resid, weights= wpfinwgt, data = wageChangesUEEE_nsw)
wageReg.nSu.nsw <- rq(residWageChange ~ unrateSA + UE + resid, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_nsw)

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

