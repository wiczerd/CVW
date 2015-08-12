# April 17, 2015
# Calculate wage change regressions:
# Requires that wageChangeStats has been run to create wageChanges.RData
library(Hmisc)
library(dplyr)
library(ggplot2)
library(xlsx)
library(quantreg)
library(stargazer)
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
	setwd("soc2d/RegResid")
} else if(useSoc2d & !useRegResid){
	setwd("soc2d/Raw")
} else if(!useSoc2d & useRegResid){
	setwd("occ/RegResid")
} else if(!useSoc2d & !useRegResid){
	setwd("occ/Raw")
}
wageChanges<-readRDS("wageChanges.RData")   

qtl_delw <- c(0.1, 0.25, .5, .75, 0.9)
Nqtl = length(qtl_delw)
## Regs on diffs ---------------------------------
wageChangesEE <- subset(wageChanges, EE)

delwageOLSEE <- lm(wageChange ~ switchedOcc + unrateSA , weights= wpfinwgt, data = wageChangesEE)
delwageQREE <- list()
for(ti in seq(1:Nqtl)){
	delwageQREE[[ti]] <- rq(wageChange ~ switchedOcc + unrateSA, tau =qtl_delw[ti], weights= wpfinwgt, data = wageChangesEE)
}
# with occupation dummies for occupation origin
delwageOLSEE.odum <- lm(wageChange ~ switchedOcc + unrateSA + factor(soc2d) + factor(nextOcc), weights= wpfinwgt, data = wageChangesEE)
delwageQREE.odum <- list()
for(ti in seq(1:Nqtl)){
	delwageQREE.odum[[ti]]  <- rq(wageChange ~ switchedOcc + unrateSA + factor(soc2d) + factor(nextOcc), tau =qtl_delw[ti], weights= wpfinwgt, data = wageChangesEE)
}

wageChangesUE <- subset(wageChanges, UE)
delwageOLSUE <- lm(wageChange ~ switchedOcc + unrateSA, weights= wpfinwgt, data = wageChangesUE)
delwageQRUE <- list()
for(ti in seq(1:Nqtl)){
	delwageQRUE[[ti]] <- rq(wageChange ~ switchedOcc + unrateSA, tau = qtl_delw[ti], weights= wpfinwgt, data = wageChangesUE)
}
# with occupation dummies for occupation origin
delwageOLSUE.odum <- lm(wageChange ~ switchedOcc + unrateSA + factor(soc2d)  + factor(nextOcc), weights= wpfinwgt, data = wageChangesUE)
delwageQRUE.odum <- list()
for(ti in seq(1:Nqtl)){
	delwageQRUE.odum[[ti]] <- rq(wageChange ~ switchedOcc + unrateSA + factor(soc2d) + factor(nextOcc), tau =qtl_delw[ti], weights= wpfinwgt, data = wageChangesUE)
}


wageChangesUEEE <- subset(wageChanges, UE | EE)
delwageOLS <- lm(wageChange ~ switchedOcc + unrateSA + UE, weights= wpfinwgt, data = wageChangesUEEE)
delwageQR <- list()
for(ti in seq(1:Nqtl)){
	delwageQR[[ti]] <- rq(wageChange ~ switchedOcc + unrateSA + UE, tau = qtl_delw[ti], weights= wpfinwgt, data = wageChangesUEEE)
}
# with occupation dummies for occupation origin/destination
delwageOLS.odum <- lm(wageChange ~ switchedOcc + unrateSA + UE + factor(soc2d) + factor(nextOcc), weights= wpfinwgt, data = wageChangesUEEE)
delwageQR.odum <- list()
for(ti in seq(1:Nqtl)){
	delwageQR.odum[[ti]] <- rq(wageChange ~ switchedOcc + unrateSA + UE + factor(soc2d)  + factor(nextOcc), tau =qtl_delw[ti], weights= wpfinwgt, data = wageChangesUEEE)
}

setwd("../../../")
setwd("./Figures")
if(useSoc2d & useRegResid) {
	setwd("soc2d/")
} else if(useSoc2d & !useRegResid){
	setwd("soc2d/")
} else if(!useSoc2d & useRegResid){
	setwd("occ/")
} else if(!useSoc2d & !useRegResid){
	setwd("occ/")
}

# make tables for combined
tab<-stargazer(delwageOLS, delwageQR[1:Nqtl] ,  out="qr_nSu.tex",
			   dep.var.labels.include = FALSE, no.space=T,  omit.stat=c("rsq","f","adj.rsq"),
			   title            = "Earnings change quantile regression",
			   digits=2, column.sep.width= "1pt",df=F,
			   covariate.labels = c("Occup Sw", "Unemp Rate","Unep Indic","Const"))

tabodum<-stargazer(delwageOLS.odum, delwageQR.odum[1:Nqtl] ,  out="qr_nSu_odum.tex",
				   dep.var.labels.include = FALSE, no.space=T, omit.stat=c("rsq","f","adj.rsq"), omit="factor",
				   omit.labels = "Next/Last Occupation Dummies",
				   digits=2, column.sep.width= "1pt",df=F,
				   title            = "Earnings change quantile regression",
				   covariate.labels = c("Occup Sw", "Unemp Rate","Unep Indic","Const"))

#make tables for EE
stargazer(delwageOLSEE, delwageQREE[1:Nqtl] ,  out="EEqr_nSu.tex",
		  dep.var.labels.include = FALSE, no.space=T, omit.stat=c("rsq","f","adj.rsq"),
		  title            = "Earnings change quantile regression, Job-job",
		  digits=2, column.sep.width= "1pt",df=F,
		  covariate.labels = c("Occup Sw", "Unemp Rate","Const"))

stargazer(delwageOLSEE.odum, delwageQREE.odum[1:Nqtl] ,  out="EEqr_nSu_odum.tex",
		  dep.var.labels.include = FALSE, no.space=T, omit.stat=c("rsq","f","adj.rsq"), omit="factor",
		  omit.labels = "Next/Last Occupation Dummies",
		  digits=2, column.sep.width= "1pt",df=F,
		  title            = "Earnings change quantile regression, Job-job",
		  covariate.labels = c("Occup Sw", "Unemp Rate","Const"))


# make tables for UE
tabUE<-stargazer(delwageOLSUE, delwageQRUE[1:Nqtl] ,  out="UEqr_nSu.tex",
				 dep.var.labels.include = FALSE, no.space=T,  omit.stat=c("rsq","f","adj.rsq"),
				 title            = "Earnings change quantile regression, through unemp",
				 digits=2, column.sep.width= "1pt",df=F,
				 covariate.labels = c("Occup Sw", "Unemp Rate","Const"))

tabUEodum<-stargazer(delwageOLSUE.odum, delwageQRUE.odum[1:Nqtl] ,  out="UEqr_nSu_odum.tex",
					 dep.var.labels.include = FALSE, no.space=T, omit.stat=c("rsq","f","adj.rsq"), omit="factor",
					 omit.labels = "Next/Last Occupation Dummies",
					 digits=2, column.sep.width= "1pt",df=F,
					 title            = "Earnings change quantile regression, through unemp",
					 covariate.labels = c("Occup Sw", "Unemp Rate","Const"))


setwd("../../")

# Regs on lags -----------------------------------
## run regs with lags rather than just differences


wageChangesUEEE <- subset(wageChanges, (UE | EE))
wageOLS <- lm(wageChange ~ unrateSA + UE + switchedOcc + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE)
wageReg <- rq(wageChange ~ unrateSA + UE + switchedOcc + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE)


# with only-changers
wageChangesUEEE_sw <- subset(wageChanges, (UE | EE) & switchedOcc)
wageOLS.sw <- lm(wageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE_sw)
wageReg.sw <- rq(wageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_sw)
ols.sw <-summary(wageOLS.sw)
qr.sw <-summary(wageReg.sw)

# only non-changer
wageChangesUEEE_nsw <- subset(wageChanges, (UE | EE) & !switchedOcc)
wageOLS.nsw <- lm(wageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, weights= wpfinwgt, data = wageChangesUEEE_nsw)
wageReg.nsw <- rq(wageChange ~ unrateSA + UE + lastResidWage + lastOccWage + occWageChange, tau = qtl_delw, weights= wpfinwgt, data = wageChangesUEEE_nsw)

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

