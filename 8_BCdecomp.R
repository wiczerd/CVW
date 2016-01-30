# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

NS = 6 #number of subgroups


DHLdecomp <- function(wcDF, rec,wc,wt){

	wcRec <- subset(wcDF, rec==T)
	wcExp <- subset(wcDF, rec==F)
	
	
	distptsRec <-with(wcRec,wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExp <-with(wcExp,wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	
	invdistRec <- approxfun(seq(0.01,0.99,by=0.005),distptsRec)
	invdistExp <- approxfun(seq(0.01,0.99,by=0.005),distptsExp)
	distRec <- approxfun(distptsRec,seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	distExp <- approxfun(distptsExp,seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	
	
	# 6 subgroups, Sw X (EE UE EU), sets up conditional distributions.
	distptsRecS <- matrix(0.,197,NS)
	distptsRecS[,1] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==T & wcRec$UE==F & wcRec$EU ==F), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsRecS[,2] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==F & wcRec$UE==T & wcRec$EU ==F), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsRecS[,3] <- with(subset(wcRec,wcRec$switchedOcc==T & wcRec$EE==F & wcRec$UE==F & wcRec$EU ==T), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsRecS[,4] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==T & wcRec$UE==F & wcRec$EU ==F), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsRecS[,5] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==F & wcRec$UE==T & wcRec$EU ==F), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsRecS[,6] <- with(subset(wcRec,wcRec$switchedOcc==F & wcRec$EE==F & wcRec$UE==F & wcRec$EU ==T), 
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS <- matrix(0.,197,6)
	distptsExpS[,1] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==T & wcExp$UE==F & wcExp$EU ==F),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS[,2] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==F & wcExp$UE==T & wcExp$EU ==F),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS[,3] <- with(subset(wcExp,wcExp$switchedOcc==T & wcExp$EE==F & wcExp$UE==F & wcExp$EU ==T),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS[,4] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==T & wcExp$UE==F & wcExp$EU ==F),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS[,5] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==F & wcExp$UE==T & wcExp$EU ==F),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExpS[,6] <- with(subset(wcExp,wcExp$switchedOcc==F & wcExp$EE==F & wcExp$UE==F & wcExp$EU ==T),
							wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	
	# create list of interp functions to make conditional distributions
	distRecS <-list(distExp,distExp,distExp,distExp,distExp,distExp)
	invdistRecS <-list(invdistExp,invdistExp,invdistExp,invdistExp,invdistExp,invdistExp)
	distExpS <-list(distExp,distExp,distExp,distExp,distExp,distExp)
	invdistExpS <-list(invdistExp,invdistExp,invdistExp,invdistExp,invdistExp,invdistExp)
	for (si in seq(1,NS)){
		invdistRecS[[si]] <- approxfun(seq(0.01,0.99,by=0.005),distptsRecS[,si])
		invdistExpS[[si]] <- approxfun(seq(0.01,0.99,by=0.005),distptsExpS[,si])
		distRecS[[si]] <- approxfun(distptsRecS[,si],seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
		distExpS[[si]] <- approxfun(distptsExpS[,si],seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	}
	
	PsExp <- rep(0,NS)
	PsRec <- rep(0,NS)
	for (swi in c(0,1)){
		PsExp[1+3*swi] = wcExp[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==T & UE==F & EU ==F,wt)]
		PsExp[2+3*swi] = wcExp[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==F & UE==T & EU ==F,wt)]
		PsExp[3+3*swi] = wcExp[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==F & UE==F & EU ==T,wt)]
		PsRec[1+3*swi] = wcRec[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==T & UE==F & EU ==F,wt)]
		PsRec[2+3*swi] = wcRec[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==F & UE==T & EU ==F,wt)]
		PsRec[3+3*swi] = wcRec[!is.na(wc), wtd.mean(switchedOcc==as.logical(1-swi) & EE==F & UE==F & EU ==T,wt)]
	}
	
	
	
	qtls<- seq(0.1,0.9,by=0.1)
	share <- array(0, dim=length(qtls))
	shift <- array(0, dim= c(length(qtls),NS))
	qi =1
	for (q in qtls){
		# set wc^E, wc^R with the two distributions
		wcE <- invdistExp(q)
		wcR <- invdistRec(q)
		fstar = 0. #integrate this
		for(si in seq(1,NS)){
			fstar = PsExp[si]*(distExpS[[si]](wcR) - distExpS[[si]](wcE)) + PsRec[si]*( distExpS[[si]](wcR) - distExpS[[si]](wcE) ) +fstar
		}
		fstar = fstar/(wcR - wcE)
		share[qi] = 0.
		for(si in seq(1,NS)){
			share[qi] = PsExp[si]*(distExpS[[si]](wcE)-distRecS[[si]](wcE) ) + PsRec[si]*(distExpS[[si]](wcR)-distRecS[[si]](wcR)) + share[qi]
		}
		share[qi] = share[qi]/fstar
		for(si in seq(1,NS)){
			shift[qi,si] = ((distExpS[[si]](wcR)-q) + (distRecS[[si]](wcE)-q))*(PsRec[si]-PsExp[si])
			shift[qi,si] = shift[qi,si]/fstar
		}
		
		qi = qi+1
	}
	
	shift_share <-list(shift,Share)
	return(shift_share)
}


wagechangesfull <- readRDS("./Data/balancedwagechanges.RData")

shift_share <- DHLdecomp(wagechangesfull,wagechangesfull$recIndic,wagechangesfull4wagechange,wagechangesfull$balanceweight)