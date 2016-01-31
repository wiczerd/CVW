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


wagechangesfull <- readRDS("./Data/balancedwagechanges.RData")


DHLdecomp <- function(wcDF,NS, recname,wcname,wtname){
# wcDF : data set
# NS   : Number of subgroups
# recname	: name of recession indicator
# wcname	: name of earnings change variable
# wtname	: name of weighting variable

	wcDF$wc  <- wcDF[[wcname]]
	wcDF$wt  <- wcDF[[wtname]]
	wcDF$rec <- wcDF[[recname]]
	wcDF$s   <- 0
	# setup subgroup indices
	if(NS ==6){
		# 6 subgroups, Sw X (EE UE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 2]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 4]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 5]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 6]
	}else if(NS==4){
		# 4 subgroups, Sw X (EE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$EU ==F, s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$EU ==T, s := 2]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$EU ==F, s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$EU ==T, s := 4]		
	}
	wcRec <- subset(wcDF, rec==T)
	wcExp <- subset(wcDF, rec==F)
	
	
	distptsRec <-with(wcRec,wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	distptsExp <-with(wcExp,wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	
	invdistRec <- approxfun(seq(0.01,0.99,by=0.005),distptsRec)
	invdistExp <- approxfun(seq(0.01,0.99,by=0.005),distptsExp)
	distRec <- approxfun(distptsRec,seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	distExp <- approxfun(distptsExp,seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	
	
	
	distptsRecS <- matrix(0.,197,NS)
	distptsExpS <- matrix(0.,197,NS)
	
	for(si in seq(1,NS)){
		distptsRecS[,si] <- with(subset(wcRec,wcRec$s == si), 
								wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
		distptsExpS[,si] <- with(subset(wcExp,wcExp$s == si),
								wtd.quantile(wc, wt, probs = seq(0.01,0.99,by=0.005)))
	}
	if(NS == 6){
		# create list of interp functions to make conditional distributions
		distRecS <-list(distExp,distExp,distExp,distExp,distExp,distExp)
		invdistRecS <-list(invdistExp,invdistExp,invdistExp,invdistExp,invdistExp,invdistExp)
		distExpS <-list(distExp,distExp,distExp,distExp,distExp,distExp)
		invdistExpS <-list(invdistExp,invdistExp,invdistExp,invdistExp,invdistExp,invdistExp)
	}else if(NS ==4){
		# create list of interp functions to make conditional distributions
		distRecS <-list(distExp,distExp,distExp,distExp)
		invdistRecS <-list(invdistExp,invdistExp,invdistExp,invdistExp)
		distExpS <-list(distExp,distExp,distExp,distExp)
		invdistExpS <-list(invdistExp,invdistExp,invdistExp,invdistExp)
	}
	
	
	for (si in seq(1,NS)){
		invdistRecS[[si]] <- approxfun(seq(0.01,0.99,by=0.005),distptsRecS[,si])
		invdistExpS[[si]] <- approxfun(seq(0.01,0.99,by=0.005),distptsExpS[,si])
		distRecS[[si]] <- approxfun(distptsRecS[,si],seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
		distExpS[[si]] <- approxfun(distptsExpS[,si],seq(0.01,0.99,by=0.005),yleft=0.,yright=1.)
	}
	
	PsExp <- rep(0,NS)
	PsRec <- rep(0,NS)
	for (si in seq(1,NS)){
		PsExp[si] = wcExp[!is.na(wc), wtd.mean(s == si,wt)]
		PsRec[si] = wcRec[!is.na(wc), wtd.mean(s == si,wt)]
	}
	
	qtls<- seq(0.1,0.9,by=0.1)
	share <- array(0, dim= length(qtls))
	wchng <- array(0, dim= length(qtls))
	shift <- array(0, dim= c(length(qtls),NS))
	qi =1
	for (q in qtls){
		# set wc^E, wc^R with the two distributions
		wcE <- invdistExp(q)
		wcR <- invdistRec(q)
		wchng[qi] = wcR - wcE
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
	shift_share <-list(wchng,shift,share)
	return(shift_share)
}

shift_share <- DHLdecomp(wagechangesfull,6,"recIndic","wagechange","balanceweight")
pct_share <- shift_share[[3]]/shift_share[[1]]
pct_shift <- shift_share[[2]]/matrix(shift_share[[1]],nrow=9,ncol=6 )

shift_share_EUE <- DHLdecomp(wagechangesfull,4,"recIndic","wagechange_EUE","balanceweight")
pct_share_EUE <- shift_share_EUE[[3]]/shift_share_EUE[[1]]
pct_shift_EUE <- shift_share_EUE[[2]]/matrix(shift_share_EUE[[1]],nrow=9,ncol=4 )

