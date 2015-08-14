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
library(quantreg)
library(reshape2)
library(xtable)

setwd("~/workspace/CVW/R")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T
useRegResid <- T


calculateCumulUse <- function(df) {
	cumulFunc <- ecdf(df$useWage)
	usepctile <- round(cumulFunc(df$useWage)*100)
	data.frame(df, usepctile)
}
calculateCumulNext <- function(df) {
	cumulFunc <- ecdf(df$nextWage)
	nextpctile <- round(cumulFunc(df$nextWage)*100)
	data.frame(df, nextpctile)
}
# Load data --------------------------------------------------------------


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
analytic9608<-readRDS( "analytic9608.RData")

wageChangesQrtileAll <- wtd.quantile(analytic9608$wageChange_all, analytic9608$wpfinwgt, na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileAll_rec <- with(analytic9608, wtd.quantile(wageChange_all[recIndic], wpfinwgt[recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) )
wageChangesQrtileAll_exp <- with(analytic9608, wtd.quantile(wageChange_all[!recIndic], wpfinwgt[!recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) )
rm(analytic9608)



# summarize changes changes of changers --------------------------------------------

wageChanges <- readRDS("wageChanges.RData")
attach(wageChanges)
wageChangesQrtile <- wtd.quantile(wageChange[EE | UE], wpfinwgt[EE | UE], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileUE <-wtd.quantile(wageChange[UE], wpfinwgt[UE], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileEE <-wtd.quantile(wageChange[EE], wpfinwgt[EE], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtile_sw <-wtd.quantile(wageChange[switchedOcc & (EE | UE)], wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 
wageChangesQrtile_sty <-wtd.quantile(wageChange[!switchedOcc & (EE | UE)], wpfinwgt[!switchedOcc & (EE | UE)], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 
wageChangesQrtileUE_sw <-wtd.quantile(wageChange[switchedOcc & UE], wpfinwgt[switchedOcc & UE], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 
wageChangesQrtileEE_sw <-wtd.quantile(wageChange[switchedOcc & EE], wpfinwgt[switchedOcc & EE], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 


wageChangesQrtile_rec <- wtd.quantile(wageChange[(EE | UE) & recIndic], wpfinwgt[(EE | UE) & recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileUE_rec <-wtd.quantile(wageChange[UE & recIndic], wpfinwgt[UE & recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileEE_rec <-wtd.quantile(wageChange[EE & recIndic], wpfinwgt[EE & recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtile_sw_rec <-wtd.quantile(wageChange[switchedOcc & (EE | UE) & recIndic], wpfinwgt[switchedOcc & (EE | UE) & recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 
wageChangesQrtileUE_sw_rec <-wtd.quantile(wageChange[switchedOcc & UE & recIndic], wpfinwgt[switchedOcc & UE & recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 

wageChangesQrtile_exp <- wtd.quantile(wageChange[(EE | UE) & !recIndic], wpfinwgt[(EE | UE) & !recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtileUE_exp <-wtd.quantile(wageChange[UE  & !recIndic], wpfinwgt[UE & !recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9)) 
wageChangesQrtile_sw_exp <-wtd.quantile(wageChange[switchedOcc & (EE | UE) & !recIndic], wpfinwgt[switchedOcc & (EE | UE) & !recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 
wageChangesQrtileUE_sw_exp <-wtd.quantile(wageChange[switchedOcc & UE & !recIndic], wpfinwgt[switchedOcc & UE & !recIndic], na.rm = TRUE,probs=c(.1,.25,.5,.75,.9) ) 

setwd("../../../")
if(useSoc2d) {
	setwd("./Figures/soc2d")
}else{
	setwd("./Figures/occ")
	}
wgCh <- rbind(wageChangesQrtileAll,wageChangesQrtile,wageChangesQrtile_sty,wageChangesQrtile_sw)
rownames(wgCh) <- c("Full sample","Job Changers", "Occupation Stayers", "Occuption Switchers")
wgCh.xt <- xtable(wgCh,label="tab:wgCh",digits=3,caption="Month-to-Month Earnings Changes")
print(wgCh.xt,file="wgCh.tex",hline.after=c(-1,-1,0,2,nrow(wgCh.xt)) )

wgCh_recexp <- rbind(wageChangesQrtileAll_rec,wageChangesQrtileAll_exp,wageChangesQrtile_rec,wageChangesQrtile_exp,wageChangesQrtile_sw_rec,wageChangesQrtile_sw_exp)
rownames(wgCh_recexp) <- c("Full sample, Rec","Full sample, Exp","Job Changers, Rec","Job Changers, Exp", "Occupation Switchers, Exp", "Occuption Switchers, Rec")
wgCh_recexp.xt <- xtable(wgCh_recexp,label="tab:wgCh_recexp",digits=3,caption="Month-to-Month Earnings Changes Over the Cycle")
print(wgCh_recexp.xt,file="wgCh_recexp.tex",hline.after=c(-1,-1,0,2,4,nrow(wgCh_recexp.xt)) )

wgCh_UE <- rbind(wageChangesQrtileEE,wageChangesQrtileUE,wageChangesQrtileEE_sw,wageChangesQrtileUE_sw,wageChangesQrtileEE_rec,wageChangesQrtileUE_rec)
rownames(wgCh_UE) <- c("E->E","E->U->E, Exp","E->E, Occupation Switchers","E->U->E, Occupation Switchers", "E->E, Rec", "E->U->E, Rec")
wgCh_UE.xt <- xtable(wgCh_UE,label="tab:wgCh_UE",digits=3,caption="Month-to-Month Earnings Changes, by Empoyment Status")
print(wgCh_UE.xt,file="wgCh_UE.tex",hline.after=c(-1,-1,0,2,4,nrow(wgCh_UE.xt)) )


# Wage change ladder -----------------------------------------------

wageChanges <- wageChanges %>%
	filter(!is.na(soc2d) & lfStat==1) %>%
	group_by(soc2d) %>%
	do(calculateCumulUse(.)) %>%
	ungroup()
if(useSoc2d){
	wageChanges <- wageChanges %>%
		group_by(id) %>%
		arrange(id,date) %>%
		mutate(nextOcc = lead(soc2d)) %>%
		ungroup() %>%
		group_by(nextOcc) %>%
		do(calculateCumulNext(.)) 
}


ggPr<-ggplot(wageChanges, aes(usepctile, nextpctile)) + geom_smooth(size=2) +
	ggtitle("Dynamics of Earnings") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggPr <- ggPr + geom_line( aes(y=usepctile,x=usepctile), size=1, linetype="dotted")
ggsave("DynPctile.eps", width = 6, height = 4)
ggsave("DynPctile.png", width = 6, height = 4)

ggPr<-ggplot(subset(wageChanges,UE), aes(usepctile, nextpctile)) + geom_smooth(size=2) +
	ggtitle("Dynamics of Earnings, Through Unemployment") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggPr <- ggPr + geom_line( aes(y=usepctile,x=usepctile), size=1, linetype="dotted")
ggPr <- ggPr + geom_smooth(data = wageChanges, aes(usepctile, nextpctile),size=1, colour="coral", linetype="dashed")
ggsave("DynPctileUE.eps", width = 6, height = 4)
ggsave("DynPctileUE.png", width = 6, height = 4)

ggPr<-ggplot(subset(wageChanges,EE), aes(usepctile, nextpctile)) + geom_smooth(size=2) +
	ggtitle("Dynamics of Earnings, Job-Job") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggPr <- ggPr + geom_line( aes(y=usepctile,x=usepctile), size=1, linetype="dotted")
ggPr <- ggPr + geom_smooth(data = wageChanges, aes(usepctile, nextpctile),size=1, colour="coral", linetype="dashed")
ggsave("DynPctileEE.eps", width = 6, height = 4)
ggsave("DynPctileEE.png", width = 6, height = 4)

ggPr<-ggplot(subset(wageChanges,recIndic), aes(usepctile, nextpctile)) + geom_smooth(size=2) +
	ggtitle("Dynamics of Earnings, Recession") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggPr <- ggPr + geom_line( aes(y=usepctile,x=usepctile), size=1, linetype="dotted")
ggPr <- ggPr + geom_smooth(data = wageChanges, aes(usepctile, nextpctile),size=1, colour="coral", linetype="dashed")
ggsave("DynPctileRec.eps", width = 6, height = 4)
ggsave("DynPctileRec.png", width = 6, height = 4)

ggPr<-ggplot(subset(wageChanges,!recIndic), aes(usepctile, nextpctile)) + geom_smooth() +
	ggtitle("Dynamics of Earnings, Expansion") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggsave("DynPctileExp.eps", width = 6, height = 4)
ggsave("DynPctileExp.png", width = 6, height = 4)


ggPr<-ggplot(subset(wageChanges,switchedOcc), aes(usepctile, nextpctile)) + geom_smooth() +
	ggtitle("Dynamics of Earnings, Occupation Switchers") + 
	ylab("Earnings percentile within next occupation") + xlab("Earnings percentile within two-digit SOC occupation")
ggPr <- ggPr + geom_line( aes(y=usepctile,x=usepctile), size=1, linetype="dotted")
ggPr <- ggPr + geom_smooth(data = wageChanges, aes(usepctile, nextpctile),size=1, colour="coral", linetype="dashed")
ggsave("DynPctileSw.eps", width = 6, height = 4)
ggsave("DynPctileSw.png", width = 6, height = 4)


# OTHER MOMENTS -----------------------------------------------------
# Mean wage changes
with(wageChanges, wtd.mean(wageChange[switchedOcc & (EE | UE)], 
                           wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[switchedOcc & EE], 
                           wpfinwgt[switchedOcc & EE], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[switchedOcc & UE], 
                           wpfinwgt[switchedOcc & UE], na.rm = TRUE))
# no change
with(wageChanges, wtd.mean(wageChange[!switchedOcc & (EE | UE)], 
						   wpfinwgt[!switchedOcc & (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[!switchedOcc & EE], 
						   wpfinwgt[!switchedOcc & EE], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[!switchedOcc & UE], 
						   wpfinwgt[!switchedOcc & UE], na.rm = TRUE))
# tot
with(wageChanges, wtd.mean(wageChange[ (EE | UE)], 
						   wpfinwgt[ (EE | UE)], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[ EE], 
						   wpfinwgt[ EE], na.rm = TRUE))
with(wageChanges, wtd.mean(wageChange[ UE], 
						   wpfinwgt[ UE], na.rm = TRUE))


# Standard deviation of wage changes
with(wageChanges, sqrt(wtd.var(wageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(wageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], na.rm = TRUE)))
with(wageChanges, sqrt(wtd.var(wageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], na.rm = TRUE)))

# Median of wage changes
with(wageChanges, wtd.quantile(wageChange[switchedOcc & (EE | UE)], 
                               wpfinwgt[switchedOcc & (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(wageChange[switchedOcc & EE], 
                               wpfinwgt[switchedOcc & EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(wageChange[switchedOcc & UE], 
                               wpfinwgt[switchedOcc & UE], probs = c(.25,.5,0.75 ) ))
#no switch
with(wageChanges, wtd.quantile(wageChange[!switchedOcc & (EE | UE)], 
							   wpfinwgt[!switchedOcc & (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(wageChange[!switchedOcc & EE], 
							   wpfinwgt[!switchedOcc & EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(wageChange[!switchedOcc & UE], 
							   wpfinwgt[!switchedOcc & UE], probs = c(.25,.5,0.75 ) ))
#tot
with(wageChanges, wtd.quantile(wageChange[ (EE | UE)], 
							   wpfinwgt[ (EE | UE)], probs = c(.25,.5,0.75 ) ) )
with(wageChanges, wtd.quantile(wageChange[ EE], 
							   wpfinwgt[ EE], probs = c(.25,.5,0.75 ) ))
with(wageChanges, wtd.quantile(wageChange[ UE], 
							   wpfinwgt[ UE], probs = c(.25,.5,0.75 ) ))

# Positive and negative wage changes -------------------
wageChanges <- wageChanges %>%
	mutate(posChange = (wageChange > 0)) %>%
	mutate(negChange = (wageChange < 0))


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


