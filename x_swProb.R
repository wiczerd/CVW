library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)
library(xtable)
library(quantreg)
library(ggplot2)
library(texreg)

wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

wagechanges <- readRDS("./Data/balancedwagechanges.RData")
CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
CPSunempRt$unrt <- CPSunempRt$unrt/100

wagechanges <- merge(wagechanges, CPSunempRt, by = "date", all.x = TRUE)


toKeep <- c("switchedOcc",
			"age","female","black","hispanic","union",			
			"HSCol",
			"recIndic",
			"wagechange",
			"wagechange_EUE", 
			"wagechange_all", 
			"balanceweight", 
			"EE","EU","UE",
			"unrt","maxunempdur",
			"wave","id")

# select toKeep columns only
wagechanges <- wagechanges[, toKeep, with = FALSE]

DTall <- readRDS("./Data/DTall_6.RData")
DTall <- merge(DTall, CPSunempRt, by = "date", all.x = TRUE)

toKeep <- c(toKeep,"wpfinwgt","switchedJob")


# select toKeep columns only
DTall <- DTall[, toKeep, with = FALSE]
DTall <- subset(DTall, is.finite(wpfinwgt) & is.finite(wagechange_all))

DTall[, allwt := wpfinwgt]
DTall[EU==T|UE==T|EE==T, allwt := balanceweight]
DTall[, wagechange_allEUE := ifelse(EU==T, wagechange_EUE,wagechange_all)]
DTall[UE==T, wagechange_allEUE := NA_real_]
DTall[, allwtEUE := allwt]
DTall[EU==T, allwtEUE := allwtEUE*2.]
DTall[UE==T, allwtEUE := 0.]

DTall<-DTall[ is.finite(EE)&is.finite(EU)&is.finite(UE),]

# wage change probit -------------------

# balanceweights, mean1:
wagechanges[, balanceweightEUE:=  balanceweight]
wagechanges[EU==T, balanceweightEUE:=  2*balanceweight]
wagechanges[UE==T, balanceweightEUE:=  0.]

swIndiclm <- lm( switchedOcc ~ recIndic + EU+ age + I(age^2/100) + female + black + hispanic + factor(HSCol, levels=c(1,0,2))
				, weights=balanceweight ,data=subset(wagechanges,EU|EE))
swIndic <- glm( switchedOcc ~ recIndic + EU+ age + I(age^2/100) + female + black + hispanic + factor(HSCol, levels=c(1,0,2))
	 , family=binomial(link=probit), weights=balanceweight ,data=subset(wagechanges,EU|EE))
swIndic$estcoef <- swIndic$coefficients
prSw <- wagechanges[EU|EE, wtd.mean(switchedOcc,weights=balanceweight,na.rm=T)]
for(ci in seq(1,length(swIndic$coefficients))){
	swIndic$coefficients[ci] <- swIndic$estcoef[ci] * pnorm( mean(swIndic$linear.predictors) )
}

swUnrtlm <- lm( switchedOcc ~ unrt + EU+ age + I(age^2/100) + female + black + hispanic + factor(HSCol, levels=c(1,0,2))
			   , weights=balanceweight ,data=subset(wagechanges,EU|EE))
swUnrt <- glm( switchedOcc ~ unrt + EU+ age + I(age^2/100) + female + black + hispanic + factor(HSCol, levels=c(1,0,2))
				, family=binomial(link=probit), weights=balanceweight ,data=subset(wagechanges,EU|EE))
swUnrt$estcoef <- swUnrt$coefficients

for(ci in seq(1,length(swUnrt$coefficients))){
	swUnrt$coefficients[ci] <- swUnrt$estcoef[ci] * pnorm( mean(swUnrt$linear.predictors) )
}
swUnrtlmEE <- lm( switchedOcc ~ unrt + age + I(age^2/100) + female + black + hispanic + factor(HSCol, levels=c(1,0,2))
				, weights=balanceweight ,data=subset(wagechanges,EE))
swUnrtlmEU <- lm( switchedOcc ~ unrt + age + I(age^2/100) + female + black + hispanic +maxunempdur + factor(HSCol, levels=c(1,0,2))
				  , weights=balanceweight ,data=subset(wagechanges,EU))

Nc<- length(swUnrtlm$coefficients)
texreg(list(swUnrtlm,swIndiclm,swUnrtlmEE,swUnrtlmEU),custom.model.names=c("All","All","EE","EUE"), reorder.coef=c(2,Nc+1,3,Nc+2,seq(4,Nc),1) ,
	   custom.coef.names=c("Const","Unemp Rt","Unemp Indic","Age","Age$^2$/100","Female","Black","Hispanic","<HS","Col+","Rec Indic","Unemp Dur"),
	   file="./Figures/swReg.tex")


