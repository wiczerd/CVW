# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(xtable)
library(Hmisc)
library(quantreg)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Results"

setwd(wd0)

keep <- c("wagechange","wagechange_EUE","EU","UE","EE","recIndic","switchedOcc","switchedInd","balanceweight","date")

mmtabchngqtls <- seq(0.1,0.9,0.1)
mmtaballqtls  <- c(seq(.01,.1,.03),seq(0.25,0.75,0.25),seq(.9,.99,.03))
Nsims = 10

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
	}else if(NS ==7){
		# 6 subgroups, Sw X (EE UE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 2]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 4]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 5]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 6]
		wcDF[                      wcDF$EE==F & wcDF$UE==F & wcDF$EU ==F, s := 7]
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
	# allocate lists, and will grow them later
	distRecS <-list(distRec) 
	invdistRecS <-list(invdistRec)
	distExpS <-list(distExp)
	invdistExpS <-list(invdistExp)
	
	
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

MMdecomp <- function(wcDF,NS,recname,wcname,wtname){
	# wcDF : data set
	# NS   : Number of subgroups
	# recname	: name of recession indicator
	# wcname	: name of earnings change variable
	# wtname	: name of weighting variable
	
	wcDF$wc  <- wcDF[[wcname]]
	wcDF$wt  <- wcDF[[wtname]]
	wcDF$rec <- wcDF[[recname]]
	wcDF$s   <- NA_integer_
	
	# setup subgroup indices
	if(NS ==6){
		# 6 subgroups, Sw X (EE UE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 2]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F, s := 4]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F, s := 5]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T, s := 6]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(wcDF$s), s5 := ifelse(s==5,1,0)]
		wcDF[!is.na(wcDF$s), s6 := ifelse(s==6,1,0)]		
	}else if(NS==7){
		# 7 subgroups, Sw X (EE UE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F , s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F , s := 2]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T , s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F , s := 4]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F , s := 5]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T , s := 6]
		wcDF[                    !(wcDF$EE==T | wcDF$UE==T | wcDF$EU ==T), s := 7]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(wcDF$s), s5 := ifelse(s==5,1,0)]
		wcDF[!is.na(wcDF$s), s6 := ifelse(s==6,1,0)]		
		wcDF[!is.na(wcDF$s), s7 := ifelse(s==7,1,0)]
	}else if(NS==4){
		# 4 subgroups, Sw X (EE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$EU ==F, s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$EU ==T, s := 2]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$EU ==F, s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$EU ==T, s := 4]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
	}else if(NS==5){
		# 4 subgroups, Sw X (EE EU), sets up conditional distributions.
		wcDF[wcDF$switchedOcc==T & wcDF$EE==T & wcDF$EU ==F , s := 1]
		wcDF[wcDF$switchedOcc==T & wcDF$EE==F & wcDF$EU ==T , s := 2]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==T & wcDF$EU ==F , s := 3]
		wcDF[wcDF$switchedOcc==F & wcDF$EE==F & wcDF$EU ==T , s := 4]
		wcDF[                    !(wcDF$EE==T | wcDF$EU ==T), s := 5]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(wcDF$s), s5 := ifelse(s==4,1,0)]		
	}
	wcDF<- subset(wcDF,is.finite(wcDF$s))
	wcRec <- subset(wcDF, rec==T)
	wcExp <- subset(wcDF, rec==F)
	
	# run quantile regressions on a relatively coarse grid (have to run 8 regressions)
	if(NS == 5 | NS==7){
		qtlgrid <- mmtaballqtls
	}else{
		qtlgrid <- c(0.01,0.05,mmtabchngqtls,0.95,0.99)
	}
	#qtlgrid <- seq(0.01,0.99,0.01)
	
	betaptsR <- matrix(0.,nrow = length(qtlgrid),ncol=NS)
	betaptsE <- matrix(0.,nrow = length(qtlgrid),ncol=NS)

	#rhere <- rq( as.formula(paste0("wc~",paste(grep("s[1-9]", names(wcDF),value=T), collapse="+"),"+0") ),
	#			 tau= q, data=wcRec, weights = wt)
	#rhere <- rq( wc ~ factor(s)+0 ,tau= qtlgrid, wcRec, weights = wt, method="pfn", Mm.factor = 1.6)
	rhere <- rq( wc ~ factor(s)+0 ,tau= qtlgrid, wcRec, weights = wt, method="br")
	betaptsR = t(rhere$coefficients)
	#rhere <- rq( as.formula(paste0("wc~",paste(grep("s[1-9]", names(wcDF),value=T), collapse="+"),"+0") ), 
	#			 tau= q, data=wcExp, weights = wt)
	rhere <- rq( wc ~factor(s) +0,tau= qtlgrid, data=wcExp, weights = wt, method="br")
	betaptsE = t(rhere$coefficients)
	rm(rhere)

	if(NS==6){
		colnames(betaptsE) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw")
		colnames(betaptsR) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw")
	}else if(NS==7){
		colnames(betaptsE) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay")
		colnames(betaptsR) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay")
	}
	
	qtlgridOut <- seq(0.01,0.99,0.01)
	betaE <- array(0.,dim=c(NS,length(qtlgridOut)) )
	betaR <- array(0.,dim=c(NS,length(qtlgridOut)) )
	for(si in seq(1,NS)){
		betaE[si,] <- spline(x=qtlgrid, y=betaptsE[,si], method="hyman", xout=qtlgridOut)$y
		betaR[si,] <- spline(x=qtlgrid, y=betaptsR[,si], method="hyman", xout=qtlgridOut)$y
	}
	
	qtlgrid <- qtlgridOut
	nsamp = floor(nrow(wcExp)/length(qtlgrid))

	seedint = 941987
	set.seed(seedint)
	# initialize distributions:
	wc_cf_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_cf_sw_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_cf_un_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_rec_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_exp_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	
	for( simiter in seq(1,Nsims)){	
		set.seed(seedint + 13*simiter)
		qi = 1
		wc_cf <- matrix(NA, nrow=nsamp*length(qtlgrid),ncol=1) #storing the counter-factual distribution
		for(q in qtlgrid){
			if(NS == 6){
				wc_cf[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
					  	betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
					  	betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
			}else if(NS==7){
				wc_cf[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
						betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 + 
						betaE[7,qi]*s7]
			}else if(NS==4){
				wc_cf[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
						betaE[3,qi]*s3 + betaE[4,qi]*s4]
			}else if(NS==5){
				wc_cf[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
						betaE[3,qi]*s3 + betaE[4,qi]*s4 +
						betaE[5,qi]*s5	]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_cf_pctile[,simiter] <- quantile(wc_cf,probs=qtlgrid,na.rm=T)
		rm(wc_cf)
		set.seed(seedint + 13*simiter)
		qi=1
		wc_cf_sw <- matrix(NA, nrow=nsamp*length(qtlgrid),ncol=1) #storing the counter-factual distribution - only switch
		for(q in qtlgrid){
			if(NS == 6){
				wc_cf_sw[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
							betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
							betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6] 
			}else if(NS == 7){
				wc_cf_sw[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
							betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
							betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 + 
							betaE[7,qi]*s7]
			}else if(NS==4){
				wc_cf_sw[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
							betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
							betaR[3,qi]*s3 + betaR[4,qi]*s4]
			}else if(NS==5){
				wc_cf_sw[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
							betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
							betaR[3,qi]*s3 + betaR[4,qi]*s4 + 
							betaE[5,qi]*s5]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_cf_sw_pctile[,simiter] <- quantile(wc_cf_sw,probs=qtlgrid,na.rm=T)
		rm(wc_cf_sw)
		
		set.seed(seedint + 13*simiter)
		qi=1
		wc_cf_un <- matrix(NA, nrow=nsamp*length(qtlgrid),ncol=1) #storing the counter-factual distribution - only unemployment
		for(q in qtlgrid){
			if(NS == 6){
				wc_cf_un[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
							betaR[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
							betaR[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
			}else if(NS == 7){
				wc_cf_un[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
							betaR[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
							betaR[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 +
							betaR[7,qi]*s7]
			}else if(NS==4){
				wc_cf_un[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
							betaR[1,qi]*s1 + betaE[2,qi]*s2 + 
							betaR[3,qi]*s3 + betaE[4,qi]*s4]
			}else if(NS==5){
				wc_cf_un[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
							betaR[1,qi]*s1 + betaE[2,qi]*s2 + 
							betaR[3,qi]*s3 + betaE[4,qi]*s4 +
							betaR[5,qi]*s5]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_cf_un_pctile[,simiter] <- quantile(wc_cf_un,probs=qtlgrid,na.rm=T)
		rm(wc_cf_un)
		
		set.seed(seedint + 13*simiter)
		qi=1
		wc_rec <- matrix(NA, nrow=nsamp*length(qtlgrid),ncol=1) #storing the counter-factual distribution - only unemployment
		for(q in qtlgrid){
			if(NS == 6){
				wc_rec[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
					betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
					betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6] 
			}else if(NS == 7){
				wc_rec[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt),
					betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
					betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 +
					betaR[7,qi]*s7]
			}else if(NS==4){
				wc_rec[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
					betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
					betaR[3,qi]*s3 + betaR[4,qi]*s4]
			}else if(NS==5){
				wc_rec[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcRec[sample(nrow(wcRec), nsamp ,replace =T, prob=wt), 
					betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
					betaR[3,qi]*s3 + betaR[4,qi]*s4 +
					betaR[5,qi]*s5]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_rec_pctile[,simiter] <- quantile(wc_rec,probs=qtlgrid,na.rm=T)
		rm(wc_rec)
	
		set.seed(seedint + 13*simiter)
		qi=1
		wc_exp <- matrix(NA, nrow=nsamp*length(qtlgrid),ncol=1) #storing the counter-factual distribution - only unemployment
		for(q in qtlgrid){
			if(NS == 6){
				wc_exp[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcExp[sample(nrow(wcExp), nsamp ,replace =T, prob=wt),
															   betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
															   betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
			}else if(NS == 7){
				wc_exp[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcExp[sample(nrow(wcExp), nsamp ,replace =T, prob=wt),
															   betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
															   betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 +
															   betaE[7,qi]*s7]
			}else if(NS==4){
				wc_exp[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcExp[sample(nrow(wcExp), nsamp ,replace =T, prob=wt), 
															   betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
															   betaE[3,qi]*s3 + betaE[4,qi]*s4]
			}else if(NS==5){
				wc_exp[ ((qi-1)*nsamp+1):(qi*nsamp) ] <- wcExp[sample(nrow(wcExp), nsamp ,replace =T, prob=wt), 
															   betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
															   betaE[3,qi]*s3 + betaE[4,qi]*s4 +
															   betaE[5,qi]*s5]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_exp_pctile[,simiter] <- quantile(wc_exp,probs=qtlgrid,na.rm=T)
		rm(wc_exp)
	}	
	
		
	return(list(betaptsE = betaptsE,betaptsR=betaptsR,wc_cf = wc_cf_pctile, wc_cf_sw= wc_cf_sw_pctile, wc_cf_un= wc_cf_un_pctile,
				wc_rec = wc_rec_pctile,wc_exp = wc_exp_pctile))
}


########################################################################
########################################################################
## Seams version -------------------------------

library(ggplot2)
DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))
DTseam[wagechange_wave_bad==F &  !(EU_wave==T|UE_wave==T|EE_wave==T) & lfstat_wave==1, stayer:= T]
DTseam[wagechange_wave_bad==F &  (EU_wave==T|UE_wave==T|EE_wave==T)  , changer:= T]
DTseam[changer==T, stayer:= F]
DTseam[stayer==T , changer:=F]
DTseam <- DTseam[(changer|stayer)]

toKeep <- c("waveweight","wavetruncweight","wpfinwgt","EU_wave","UE_wave","EE_wave","switchedOcc_wave","wagechange_wave","wagechangeEUE_wave","recIndic_wave","date")


# select toKeep columns only
DTseam <- DTseam[, toKeep, with = FALSE]
DTseam <- subset(DTseam, is.finite(waveweight) & is.finite(wagechange_wave) & is.finite(EU_wave) & is.finite(UE_wave)& is.finite(EE_wave))
DTseam[ , EE := EE_wave]
DTseam[ , EU := EU_wave]
DTseam[ , UE := UE_wave]
DTseam[ , switchedOcc := switchedOcc_wave]
DTseamchng <- subset(DTseam, EU==T|UE==T|EE==T)


MM_wavechng_betaE_betaR_cf    <- MMdecomp(DTseamchng,6,"recIndic_wave","wagechange_wave","wavetruncweight")
MM_waveall_betaE_betaR_cf <- MMdecomp(DTseam,7,"recIndic_wave","wagechange_wave","wavetruncweight")
MM_waveallEUE_betaE_betaR_cf <- MMdecomp(DTseam,5,"recIndic_wave","wagechangeEUE_wave","wavetruncweight")

#all observations (7 categories)
wcExp <- subset(DTseam,recIndic_wave==F)
wcRec <- subset(DTseam,recIndic_wave==T)
dist_exp      <- wcExp[ , wtd.quantile(wagechange_wave,probs=seq(0.01,.99,0.01),weights=wavetruncweight, na.rm=T)]
dist_rec      <- wcRec[ , wtd.quantile(wagechange_wave,probs=seq(0.01,.99,0.01),weights=wavetruncweight, na.rm=T)]
dist_cf       <- rowMeans(MM_waveall_betaE_betaR_cf$wc_cf   )
dist_cf_sw    <- rowMeans(MM_waveall_betaE_betaR_cf$wc_cf_un)
dist_cf_un    <- rowMeans(MM_waveall_betaE_betaR_cf$wc_cf_sw)
dist_sim_rec  <- rowMeans(MM_waveall_betaE_betaR_cf$wc_rec  )
dist_sim_exp  <- rowMeans(MM_waveall_betaE_betaR_cf$wc_exp  )

#moving window to summarize the quantiles
win_cf      <- array(0.,dim=89)
win_cf_un   <- array(0.,dim=89)
win_cf_sw   <- array(0.,dim=89)
win_exp     <- array(0.,dim=89)
win_rec     <- array(0.,dim=89)
win_rec_dat <- array(0.,dim=89)
win_exp_dat <- array(0.,dim=89)
di = 1
for(lb in seq(1,89,1)){
	ub = lb+10
	win_cf[di]      <- mean(dist_cf     [ lb:ub ])
	win_rec[di]     <- mean(dist_sim_rec[ lb:ub ])
	win_exp[di]     <- mean(dist_sim_exp[ lb:ub ])
	win_cf_un[di]   <- mean(dist_cf_un  [ lb:ub ])
	win_cf_sw[di]   <- mean(dist_cf_sw  [ lb:ub ])	
	win_rec_dat[di] <- mean(dist_rec    [ lb:ub ])
	win_exp_dat[di] <- mean(dist_exp    [ lb:ub ])	
	di= di+1
}
win_dif_rec_exp_dat <- win_rec_dat - win_exp_dat
win_dif_rec_exp     <- win_rec     - win_exp
win_dif_cf_exp      <- win_cf      - win_exp

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EUE wage changes  observations (5 categories)
distEUE_exp      <- wcExp[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.01,.99,0.01),weights=wavetruncweight, na.rm=T)]
distEUE_rec      <- wcRec[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.01,.99,0.01),weights=wavetruncweight, na.rm=T)]
distEUE_cf       <- rowMeans(MM_waveallEUE_betaE_betaR_cf$wc_cf   )
distEUE_cf_sw    <- rowMeans(MM_waveallEUE_betaE_betaR_cf$wc_cf_un)
distEUE_cf_un    <- rowMeans(MM_waveallEUE_betaE_betaR_cf$wc_cf_sw)
distEUE_sim_rec  <- rowMeans(MM_waveallEUE_betaE_betaR_cf$wc_rec  )
distEUE_sim_exp  <- rowMeans(MM_waveallEUE_betaE_betaR_cf$wc_exp  )

#moving window to summarize the quantiles
winEUE_cf      <- array(0.,dim=89)
winEUE_cf_un   <- array(0.,dim=89)
winEUE_cf_sw   <- array(0.,dim=89)
winEUE_exp     <- array(0.,dim=89)
winEUE_rec     <- array(0.,dim=89)
winEUE_rec_dat <- array(0.,dim=89)
winEUE_exp_dat <- array(0.,dim=89)
di = 1
for(lb in seq(1,89,1)){
	ub = lb+10
	winEUE_cf[di]      <- mean(distEUE_cf     [ lb:ub ])
	winEUE_rec[di]     <- mean(distEUE_sim_rec[ lb:ub ])
	winEUE_exp[di]     <- mean(distEUE_sim_exp[ lb:ub ])
	winEUE_cf_un[di]   <- mean(distEUE_cf_un  [ lb:ub ])
	winEUE_cf_sw[di]   <- mean(distEUE_cf_sw  [ lb:ub ])	
	winEUE_rec_dat[di] <- mean(distEUE_rec    [ lb:ub ])
	winEUE_exp_dat[di] <- mean(distEUE_exp    [ lb:ub ])	
	di= di+1
}
winEUE_dif_rec_exp_dat <- winEUE_rec_dat - winEUE_exp_dat
winEUE_dif_rec_exp     <- winEUE_rec     - winEUE_exp
winEUE_dif_cf_exp      <- winEUE_cf      - winEUE_exp


##################################################
### plot stuff
##################################################
dt_mm <- data.table(cbind( seq(0.06,0.94,0.01),win_dif_rec_exp_dat ))
names(dt_mm) <- c("Quantile","Data")
dt_mm_melted <- melt(dt_mm, id= "Quantile")
ggplot( dt_mm_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	scale_color_manual(values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[1])) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_data.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_data.png",height=5,width=10)


dt_mm <- data.table(cbind( seq(0.06,0.94,0.01),win_dif_rec_exp_dat,win_dif_cf_exp ))
names(dt_mm) <- c("Quantile","Data","Counter Factual")
dt_mm_melted <- melt(dt_mm, id= "Quantile")
ggplot( dt_mm_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[1:2]) ) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all.eps",height=5,width=10)
ggsave("./Figures/MMwave_all.png",height=5,width=10)

dt_mm_un <- data.table(cbind( seq(0.05,0.95,0.01),spn_dif_rec_exp_dat,spn_dif_cf_exp,spn_dif_un_exp ))
names(dt_mm_un) <- c("Quantile","Data","Counter Factual","Recession Switching Coefficients")
dt_mm_un_melted <- melt(dt_mm_un, id= "Quantile")
ggplot( dt_mm_un_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	scale_color_manual(values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[1:3])) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.6,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_un.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_un.png",height=5,width=10)

dt_mm_sw <- data.table(cbind( seq(0.05,0.95,0.01),spn_dif_rec_exp_dat,spn_dif_rec_exp,spn_dif_sw_exp ))
names(dt_mm_sw) <- c("Quantile","Counter Factual","Data","Only Unemployment Coefficients")
dt_mm_sw_melted <- melt(dt_mm_sw, id= "Quantile")
ggplot( dt_mm_sw_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	#scale_color_manual(values = c(hcl(h=seq(15, 375, length=Ndemo+1), l=50, c=100)[1:Ndemo], "black")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_sw.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_sw.png",height=5,width=10)



wcExp <- subset(DTseamchng,recIndic_wave==F)
wcRec <- subset(DTseamchng,recIndic_wave==T)
distchng_exp    <- wcExp[, wtd.quantile(wagechange_wave,probs=mmtabchngqtls,weights=waveweight, na.rm=T)]
distchng_rec    <- wcRec[, wtd.quantile(wagechange_wave,probs=mmtabchngqtls,weights=wcRec$allwt, na.rm=T)]
distchng_sim_rec<- rowMeans(MM_wavechng_betaE_betaR_cf$wc_rec)
distchng_sim_exp<- rowMeans(MM_wavechng_betaE_betaR_cf$wc_exp)
distchng_cf     <- rowMeans(MM_wavechng_betaE_betaR_cf$wc_cf)
distchng_cf_sw  <- rowMeans(MM_wavechng_betaE_betaR_cf$wc_cf_sw)
distchng_cf_un  <- rowMeans(MM_wavechng_betaE_betaR_cf$wc_cf_un)
distchng_pct    <- (distchng_cf - distchng_exp)/(distchng_rec - distchng_exp)
distchng_pct    <- (distchng_cf - distchng_sim_exp)/(distchng_sim_rec - distchng_sim_exp)
distchng_pct_sw <- (distchng_cf_sw - distchng_exp)/(distchng_rec - distchng_exp)
distchng_pct_un <- (distchng_cf_un - distchng_exp)/(distchng_rec - distchng_exp)


#take decile averages
decchng_cf  <- array(0.,dim = 10)
decchng_exp <- array(0.,dim = 10)
decchng_rec <- array(0.,dim = 10)
decchng_exp_dat <- array(0.,dim = 10)
decchng_rec_dat <- array(0.,dim = 10)
for(di in seq(1,10,1)){
	decchng_cf [di] <- mean(distchng_cf     [ ((di-1)*10+1):min(di*10,99)])
	decchng_rec[di] <- mean(distchng_sim_rec[ ((di-1)*10+1):min(di*10,99)])
	decchng_exp[di] <- mean(distchng_sim_exp[ ((di-1)*10+1):min(di*10,99)])
	decchng_rec_dat[di] <- mean(distchng_rec[ ((di-1)*10+1):min(di*10,99)])
	decchng_exp_dat[di] <- mean(distchng_exp[ ((di-1)*10+1):min(di*10,99)])	
}

spn_difchng_cf_exp <- spline(seq(0.05,0.95,0.1), (decchng_cf-decchng_exp),xout=seq(0.05,0.95,0.01))$y
spn_difchng_rec_exp <- spline(seq(0.05,0.95,0.1), (decchng_rec-decchng_exp),xout=seq(0.05,0.95,0.01))$y

dt_mmchng <- data.table(cbind( seq(0.05,0.95,0.01),spn_difchng_cf_exp,spn_difchng_rec_exp ))
names(dt_mmchng) <- c("Quantile","Counter Factual","Data")
dt_mmchng_melted <- melt(dt_mmchng, id= "Quantile")
ggplot( dt_mmchng_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	#scale_color_manual(values = c(hcl(h=seq(15, 375, length=Ndemo+1), l=50, c=100)[1:Ndemo], "black")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",
		  								 color = "black"))
ggsave("./Figures/MMwave_chng.eps",height=5,width=10)
ggsave("./Figures/MMwave_chng.png",height=5,width=10)


waveshift_share <- DHLdecomp(DTseam,7,"recIndic_wave","wagechange_wave","waveweight")
wave_pct_share <- waveshift_share[[3]]/waveshift_share[[1]]
wave_pct_shift <- waveshift_share[[2]]/matrix(waveshift_share[[1]],nrow=9,ncol=7 )

#table out:
MM_all_tab <- data.table(cbind( mmtabqtls,(dist_cf),dist_rec,dist_exp,dist_rec- dist_exp, 
								dist_pct,dist_pct_sw,dist_pct_un))
names(MM_all_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MM_all_tab) <- mmtabqtls
MM_all_tab <- xtable(MM_all_tab, digits=2, 
				 align="ll|lll|l|lll", caption="Machado-Mata, including unemployment \\label{tab:MM_all_tab}")
print(MM_all_tab,include.rownames=T, hline.after= c(0,nrow(MM_all_tab)), file="./Figures/MM_all.tex")

#table out:
MM_allEUE_tab <- data.table(cbind( mmtabqtls,(distEUE_cf),distEUE_rec,distEUE_exp,distEUE_rec- distEUE_exp, 
								   distEUE_pct,distEUE_pct_sw,distEUE_pct_un))
names(MM_allEUE_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MM_allEUE_tab) <- mmtabqtls
MM_allEUE_tab <- xtable(MM_allEUE_tab, digits=2, 
					 align="ll|lll|l|lll", caption="Machado-Mata, connecting unemployment spells \\label{tab:MM_allEUE_tab}")
print(MM_allEUE_tab,include.rownames=T, hline.after= c(0,nrow(MM_allEUE_tab)), file="./Figures/MM_allEUE.tex")



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# only month changes data set
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
wagechanges <- readRDS(paste0(datadir,"/balancedwagechanges.RData"))
wagechanges <- subset(wagechanges, select=keep)


# create vector of recession dates : above 7%
recDates3 <- as.Date(c("1991-10-01", "1993-07-01","2008-12-01", "2013-11-01"))
wagechanges[, recIndic3 := (date > recDates3[1] & date < recDates3[2]) | 
				(date > recDates3[3] & date < recDates3[4])]


shift_share <- DHLdecomp(wagechanges,6,"recIndic","wagechange","balanceweight")
pct_share <- shift_share[[3]]/shift_share[[1]]
pct_shift <- shift_share[[2]]/matrix(shift_share[[1]],nrow=9,ncol=6 )

shift_share_EUE <- DHLdecomp(wagechanges,4,"recIndic","wagechange_EUE","balanceweightEUE")
pct_share_EUE <- shift_share_EUE[[3]]/shift_share_EUE[[1]]
pct_shift_EUE <- shift_share_EUE[[2]]/matrix(shift_share_EUE[[1]],nrow=9,ncol=4 )

MM_betaE_betaR_cf <- MMdecomp(wagechanges,6,"recIndic","wagechange","balanceweight")

MMEUE_betaE_betaR_cf <- MMdecomp(subset(wagechanges,EE==T|EU==T),4,"recIndic","wagechange_EUE","balanceweightEUE")

mmtabqtls <- seq(0.1,0.9,0.1)
wcExp <- subset(wagechanges,recIndic==F)
wcRec <- subset(wagechanges,recIndic==T)
distEUE_exp <- wcExp[EE==T | EU==T, wtd.quantile(wagechange_EUE,probs=mmtabqtls,weights=balanceweightEUE, na.rm=T)]
distEUE_rec <- wcRec[EE==T | EU==T, wtd.quantile(wagechange_EUE,probs=mmtabqtls,weights=balanceweightEUE, na.rm=T)]
distEUE_cf  <- MMEUE_betaE_betaR_cf$wc_cf[mmtabqtls*100] 
distEUE_cf_sw  <- MMEUE_betaE_betaR_cf$wc_cf_un[mmtabqtls*100]
distEUE_cf_un  <- MMEUE_betaE_betaR_cf$wc_cf_sw[mmtabqtls*100]
distEUE_pct <- (distEUE_cf - distEUE_exp)/(distEUE_rec - distEUE_exp)
distEUE_pct_sw <- (distEUE_cf_sw - distEUE_exp)/(distEUE_rec - distEUE_exp)
distEUE_pct_un <- (distEUE_cf_un - distEUE_exp)/(distEUE_rec - distEUE_exp)

dist_exp <- wtd.quantile(wcExp$wagechange,probs=mmtabqtls,weights=wcExp$balanceweight, na.rm=T)
dist_rec <- wtd.quantile(wcRec$wagechange,probs=mmtabqtls,weights=wcRec$balanceweight, na.rm=T)
dist_cf  <- MM_betaE_betaR_cf$wc_cf[mmtabqtls*100]
dist_cf_sw  <- MM_betaE_betaR_cf$wc_cf_sw[mmtabqtls*100]
dist_cf_un  <- MM_betaE_betaR_cf$wc_cf_un[mmtabqtls*100] 
dist_pct <- (dist_cf - dist_exp)/(dist_rec - dist_exp)
dist_pct_sw <- (dist_cf_sw - dist_exp)/(dist_rec - dist_exp)
dist_pct_un <- (dist_cf_un - dist_exp)/(dist_rec - dist_exp)


#-- spit these out into tables

#DHW: EE, EU, UE
pct_share_shift <- data.table(cbind(seq(0.1,0.9,0.1),pct_share,pct_shift))
names(pct_share_shift) <- c("Decile","Share","EE","UE","EU","EE","UE","EU")
addswitched <-list(pos = list(-1),command="\\hline\\hline& & \\multicolumn{3}{c}{Switched Occ} & \\multicolumn{3}{c}{Not Switched Occ} \\\\ ")
pct_share_shift <- xtable(pct_share_shift, label="tab:pct_share_shift", digits=2, 
						  align="ll|l|lll|lll", caption="Shift-Share (DHL) decomposition, including unemployment")
print(pct_share_shift,add.to.row=addswitched, include.rownames=F, hline.after= c(0,nrow(pct_shift)), file="./Figures/pct_share_shift.tex")

#DHW: EE, EUE
pct_share_shift_EUE <- data.table(cbind(seq(0.1,0.9,0.1),pct_share_EUE,pct_shift_EUE))
names(pct_share_shift_EUE) <- c("Decile","Share","EE","EUE","EE","EUE")
addswitched <-list(pos = list(-1),command="\\hline\\hline& & \\multicolumn{2}{c}{Switched Occ} & \\multicolumn{2}{c}{Not Switched Occ} \\\\ ")
pct_share_shift_EUE <- xtable(pct_share_shift_EUE, label="tab:pct_share_shift_EUE", digits=2, 
							  align="ll|l|ll|ll", caption="Shift-Share (DHL) decomposition, connecting across unemployment")
print(pct_share_shift_EUE,add.to.row=addswitched, include.rownames=F, hline.after= c(0,nrow(pct_shift)), file="./Figures/pct_share_shift_EUE.tex")

#MM : EE, EU, UE
MM_tab <- data.table(cbind( mmtabqtls,(dist_cf),dist_rec,dist_exp,dist_rec- dist_exp, 
							dist_pct,dist_pct_sw,dist_pct_un))
names(MM_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MM_tab) <- mmtabqtls
MM_tab <- xtable(MM_tab, digits=2, 
				 align="ll|lll|l|lll", caption="Machado-Mata, including unemployment \\label{tab:MM_tab}")
print(MM_tab,include.rownames=T, hline.after= c(0,nrow(MM_tab)), file="./Figures/MM.tex")
ks.test(wcRec$wagechange,wcExp$wagechange,alternative = "greater")
# plot the coefficients
MMcoef <- data.table(cbind( mmtabqtls,MM_betaE_betaR_cf$betaptsE,MM_betaE_betaR_cf$betaptsE))
names(MMcoef) <- c("Quantiles","ExpSwEE","ExpSwUE","ExpSwEU","ExpNSwEE","ExpNSwUE","ExpNSwEU","RecSwEE","RecSwUE","RecSwEU","RecNSwEE","RecNSwUE","RecNSwEU")
EUfrac <- wagechanges[recIndic==F, wtd.mean(EU,na.rm=T,weights=balanceweight)]
UEfrac <- wagechanges[recIndic==F, wtd.mean(UE,na.rm=T,weights=balanceweight)]
EEfrac <- wagechanges[recIndic==F, wtd.mean(EE,na.rm=T,weights=balanceweight)]
Swfrac <- wagechanges[recIndic==F, wtd.mean(switchedOcc,na.rm=T,weights=balanceweight)]


#MM : EE,EUE
MMEUE_tab <- data.table(cbind( mmtabqtls,(distEUE_cf),distEUE_rec,distEUE_exp,distEUE_rec- distEUE_exp, 
							   distEUE_pct,distEUE_pct_sw,distEUE_pct_un ))
names(MMEUE_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MMEUE_tab) <- mmtabqtls
MMEUE_tab <- xtable(MMEUE_tab, digits=2, 
					align="ll|lll|l|lll", caption="Machado-Mata, connecting across unemployment \\label{tab:MMEUE_tab}")
print(MMEUE_tab,include.rownames=F, hline.after= c(0,nrow(MMEUE_tab)), file="./Figures/MMEUE.tex")

ks.test(wcRec$wagechange_EUE,wcExp$wagechange_EUE,alternative = "greater")

rm(list = c("wcRec","wcExp"))


rm(list=ls())