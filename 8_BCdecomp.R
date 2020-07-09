# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(xtable)
library(Hmisc)
library(quantreg)
library(ggplot2)
library(foreign)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Figures"

setwd(wd0)


wt <- "truncweight"
wc <- "wagechange_anan"
recDef <- "recIndic2_wave"

wdur = F #will include duration in the variables in the regression

minEarn = 1040 
minLEarn = log(2*minEarn)

qtlgridEst  <- c(seq(0.005,0.03,0.005),seq(0.04,0.06,0.01),seq(.07,.13,.02),seq(0.15,0.85,0.05),seq(.87,.93,.02),seq(0.94,0.96,0.01),seq(0.97,0.995,.005))
Nestpt = length(qtlgridEst)
qtlgridOut <- c(seq(.01,0.99,0.005))
qtlgridOutTruncCor <- (qtlgridOut-min(qtlgridEst))/(max(qtlgridEst)-min(qtlgridEst)) 
MMstd_errs = F
Nmoments = 5
moment_names <- c("Mean","Median","Med Abs Dev","Groenv-Meeden","Moors")
#recession counter-factual returns beta^E * recession inidcators and beta^R * expansion indicators

wtd.mad <- function(xt,wt,md=-Inf){
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	if(is.finite(md)==F){
		md  <- wtd.quantile(xt,weights=wt,probs = 0.5,na.rm = T)
	}
	wtd.quantile(abs(xt - md),weights=wt,probs = 0.5,na.rm = T)
}
wtd.skewness <- function(xt, wt){  
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	(sum( wt*(xt-wtd.mean(xt,weights=wt))^3 ) / sum(wt) ) / wtd.var(xt,weights=wt)^(1.5)
}
wtd.GroenveldMeeden <- function(xt, wt,md=-Inf){  
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	if(is.finite(md)==F){
		md  <- wtd.quantile(xt,weights=wt,probs = 0.5,na.rm = T)
	}
	mn  <- wtd.mean(xt,weights=wt)
	(mn-md)/wtd.mean(abs(xt - md),weights=wt)
}
wtd.Moors <- function(xt, wt,samp=T){  
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	if(samp==T){
		sampidx <- sample(length(xt), length(xt)/4,replace = T)
		xt <- xt[sampidx]
		wt <- wt[sampidx]
	}
	octls<-wtd.quantile(xt,weights=wt,probs = c(seq(1/8,3/8,1/8),seq(5/8,7/8,1/8)),na.rm = T) #excludes the median
	((octls[6]-octls[4]) + (octls[3]-octls[1]))/(octls[5]-octls[2]) 
}
wtd.intMoors <- function(xt, wt,samp=T){  
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	if(samp==T){
		sampidx <- sample(length(xt), length(xt)/4, replace = T)
		xt <- xt[sampidx]
		wt <- wt[sampidx]
	}
	qtls_hr<-wtd.quantile(xt,weights=wt,probs = seq(.01,.99,.01),na.rm = T)
	kurt <- array(NA,dim=13)
	for( i in seq(1,13)){
		kurt[i] <- ((qtls_hr[100-i] - qtls_hr[50+i]) + (qtls_hr[50-i] - qtls_hr[i]))/(qtls_hr[100-2*i] - qtls_hr[2*i])
	}
	return(mean(kurt,na.rm = T))
}
wtd.kurtosis <- function(xt, wt){  
	wt  <- wt[is.na(xt)==F]
	xt      <- xt[is.na(xt)==F]
	(sum( wt*(xt-wtd.mean(xt,weights=wt))^4 ) / sum(wt) ) / wtd.var(xt,weights=wt)^(2) - 3.
}
wtd.4qtlmoments <- function(xt,wt){
	#wraps 4 quantile-based moments
	median <- wtd.quantile(xt,na.rm=T,weights=wt,probs=0.5)
	mad    <- wtd.mad(xt,wt,median)
	GrnMd  <- wtd.GroenveldMeeden(xt,wt,median)
	set.seed(12281951)
	Moors  <- wtd.Moors(xt,wt)
	return(c(median,mad,GrnMd,Moors))
}

MMdecomp <- function(wcDF,NS,recname,wcname,wtname, std_errs=F,no_occ=F,durEU=F){
	# wcDF : data set
	# NS   : Number of subgroups
	# recname	: name of recession indicator
	# wcname	: name of earnings change variable
	# wtname	: name of weighting variable
	
	wcDF$wc  <- wcDF[[wcname]]
	wcDF$wt  <- wcDF[[wtname]]
	wcDF$rec <- wcDF[[recname]]
	wcDF$s   <- NA_integer_
	if(missing(std_errs)){
		std_errs = F
	}
	if(missing(no_occ)){
		no_occ = F
	}
	if(missing(durEU)){
		durEU = F
	}
	
	wcDF <- subset(wcDF,wt>0 & is.finite(wc))
	
	# setup subgroup indices
	if(NS==7){
		# 7 subgroups, Sw X (EE UE EU) + stay, sets up conditional distributions.
		wcDF[sw==T & EEfrq==T & UEfrq==F & EUfrq ==F , s := 1]
		wcDF[sw==T & EEfrq==F & UEfrq==T & EUfrq ==F , s := 2]
		wcDF[sw==T & EEfrq==F & UEfrq==F & EUfrq ==T , s := 3]
		wcDF[sw==F & EEfrq==T & UEfrq==F & EUfrq ==F , s := 4]
		wcDF[sw==F & EEfrq==F & UEfrq==T & EUfrq ==F , s := 5]
		wcDF[sw==F & EEfrq==F & UEfrq==F & EUfrq ==T , s := 6]
		wcDF[st==T                                   , s := 7]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(s), s5 := ifelse(s==5,1,0)]
		wcDF[!is.na(s), s6 := ifelse(s==6,1,0)]		
		wcDF[!is.na(s), s7 := ifelse(s==7,1,0)]
	}
	if(NS==8){
		# 8 subgroups, Sw X (EE UE EU stay), sets up conditional distributions.
		wcDF[sw==T &  EEfrq==T                                  , s := 1]
		wcDF[sw==T &                  UEfrq==T                  , s := 2]
		wcDF[sw==T &                                  EUfrq ==T , s := 3]
		
		wcDF[sw==F &  EEfrq==T                                  , s := 4]
		wcDF[sw==F &                  UEfrq==T                  , s := 5]
		wcDF[sw==F &                                  EUfrq ==T , s := 6]
		
		wcDF[sw==T & st ==T , s := 7]
		wcDF[sw==F & st ==T , s := 8]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(s), s5 := ifelse(s==5,1,0)]
		wcDF[!is.na(s), s6 := ifelse(s==6,1,0)]		
		wcDF[!is.na(s), s7 := ifelse(s==7,1,0)]
		wcDF[!is.na(s), s8 := ifelse(s==8,1,0)]
	}else if(NS==4){
		# 4 subgroups, Sw X (EE EU), sets up conditional distributions.
		wcDF[  EEfrq==T & EUfrq ==F & UEfrq ==F , s := 1]
		wcDF[  EEfrq==F & EUfrq ==T & UEfrq ==F , s := 2]
		wcDF[  EEfrq==F & EUfrq ==F & UEfrq ==T , s := 3]
		wcDF[!(EEfrq==T | EUfrq ==T | UEfrq ==T), s := 4]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
	}else if(NS==5){
		# 4 subgroups, Sw X (EE EU) + stay, sets up conditional distributions.
		wcDF[sw==T & EEfrq==T & EUfrq ==F , s := 1]
		wcDF[sw==T & EEfrq==F & EUfrq ==T , s := 2]
		wcDF[sw==F & EEfrq==T & EUfrq ==F , s := 3]
		wcDF[sw==F & EEfrq==F & EUfrq ==T , s := 4]
		wcDF[st==T                                  , s := 5]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(s), s5 := ifelse(s==5,1,0)]		
	}else if(NS==6){
		# 6 subgroups, Sw X (EE EU Stay) , sets up conditional distributions.
		wcDF[sw==T & EEfrq==T & EUfrq ==F & UEfrq==F,  s := 1]
		wcDF[sw==T & EEfrq==F & EUfrq ==T & UEfrq==F , s := 2]
		wcDF[sw==F & EEfrq==T & EUfrq ==F & UEfrq==F , s := 3]
		wcDF[sw==F & EEfrq==F & EUfrq ==T & UEfrq==F , s := 4]
		wcDF[sw==F & st   ==T                               , s := 5]
		wcDF[sw==T & st   ==T                               , s := 6]
		
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(s), s5 := ifelse(s==5,1,0)]		
		wcDF[!is.na(s), s6 := ifelse(s==6,1,0)]
	}else if(NS==4){
		# 4 subgroups, (EE EU UE) + stay, sets up conditional distributions.
		wcDF[  EEfrq==T & EUfrq ==F & UEfrq ==F , s := 1]
		wcDF[  EEfrq==F & EUfrq ==T & UEfrq ==F , s := 2]
		wcDF[  EEfrq==F & EUfrq ==F & UEfrq ==T , s := 3]
		wcDF[!(EEfrq==T | EUfrq ==T | UEfrq ==T), s := 4]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(s), s4 := ifelse(s==4,1,0)]
	}else if(NS==3){
		# 4 subgroups, (EE EU) + stay, sets up conditional distributions.
		wcDF[  EEfrq==T & EUfrq ==F & UEfrq ==F , s := 1]
		wcDF[  EEfrq==F & EUfrq ==T , s := 2]
		wcDF[!(EEfrq==T | EUfrq ==T), s := 3]
		wcDF[!is.na(s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(s), s3 := ifelse(s==3,1,0)]
	}
	if(durEU==T){
		wcDF[ !(wcDF$EUfrq==T|wcDF$UEfrq==T), dur := 0]
		wcDF[ is.na(wcDF$dur)==T, dur := 0]
	}
	wcDF<- subset(wcDF,is.finite(wcDF$s))
	wcRec <- subset(wcDF, rec==T)
	wcExp <- subset(wcDF, rec==F)
	
	# run quantile regressions on a relatively coarse grid (have to run 8 regressions)
	betaptsR <- matrix(0.,nrow = length(qtlgridEst),ncol=NS)
	betaptsE <- matrix(0.,nrow = length(qtlgridEst),ncol=NS)

	if(std_errs ==T){
		Nsims = 50
	}else{
		Nsims = 1
	}

	qtlgridSamp <- seq(min(qtlgridOut),max(qtlgridOut),0.001)
	seedint = 941987
	set.seed(seedint)
	
	#draw the sample
	nsampE = floor(nrow(wcExp)/length(qtlgridSamp))
	nsampR = floor(nrow(wcRec)/length(qtlgridSamp))
	sampE <-matrix(0.,nrow=nsampE,ncol=length(qtlgridSamp))
	sampR <-matrix(0.,nrow=nsampR,ncol=length(qtlgridSamp))
	for(qi in seq(1,length(qtlgridSamp))){
		set.seed(seedint+ qi*nsampE)
		sampE[,qi] <- sample(nrow(wcExp),nsampE,prob=wcExp$wt,replace=T) # weight the sampling. I will also weight the regression 
		sampR[,qi] <- sample(nrow(wcRec),nsampR,prob=wcRec$wt,replace=T) #,prob=wcRec$wt
	}
	# initialize output distributions:
	wc_IR_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_IR_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_BR_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_IR_sw_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_sw_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_IR_un_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_un_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_rec_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_rec_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_exp_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_exp_moments <- array(0.,dim=c(Nmoments,Nsims))
	
	wc_BR_exstay_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_exstay_moments <- array(0.,dim=c(Nmoments,Nsims))
	
	if( durEU ==F ){
		wc_BR_exi_moments <- array(0.,dim=c(Nmoments,NS,Nsims))
		wc_BR_exi_pctile <- array(0.,dim=c(length(qtlgridOut),NS,Nsims))
		wc_BR_onlyi_pctile <- array(0.,dim=c(length(qtlgridOut),NS,Nsims))
		wc_BR_onlyi_moments <- array(0.,dim=c(Nmoments,NS,Nsims))
	}else{
		wc_BR_exi_moments <- array(0.,dim=c(Nmoments,NS+1,Nsims))
		wc_BR_exi_pctile <- array(0.,dim=c(length(qtlgridOut),NS+1,Nsims))
		wc_BR_onlyi_pctile <- array(0.,dim=c(length(qtlgridOut),NS+1,Nsims))
		wc_BR_onlyi_moments <- array(0.,dim=c(Nmoments,NS+1,Nsims))
	}
	
	for( simiter in seq(1,Nsims)){	
		
		if(std_errs == T){
			set.seed(seedint+simiter*Nsims)		
			wcRec <- subset(wcDF, rec==T)
			wcExp <- subset(wcDF, rec==F)
			datsampR <- sample(nrow(wcRec),nrow(wcRec),replace=T) # <-,prob=wcRec$wt do not weight the sample. I will weight the regression and stuff
			datsampE <- sample(nrow(wcExp),nrow(wcExp),replace=T) #   ,prob=wcExp$wt
			wcRec <- wcRec[datsampR,]
			wcExp <- wcExp[datsampE,]
		}
		# here is where the quantile regression happens
		if(durEU==T){ 
			#do we include duration in the regression?
			regform <- formula(paste(c("wc~factor(s)","dur","0"),collapse=" + ") )
		}else{
			#this form is just the dummies:
			regform <- formula(paste(c("wc~factor(s)","0"),collapse=" + ") )
		}
		rhere <- rq( regform ,tau= qtlgridEst, data=wcRec, weights = wt, method="sfn")
		betaptsR = t(rhere$coefficients)
		
	
		rhere <- rq( regform ,tau= qtlgridEst, data=wcExp, weights = wt, method="sfn") #sfn
		betaptsE = t(rhere$coefficients)
		rm(rhere)
	
		if(NS==7){
			names <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay")
		}else if(NS==8){
			names <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay_nosw","stay_sw")
		}else if(NS==5){
			names <-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay")
		}else if(NS==6){
		  names <- c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay_nosw","stay_sw")
		}else if(NS==4){
			names <- c("EE","UE","EU","stay")
		}else if(NS==3){
			names <-c("EE","EU","stay")
		}
		if(durEU==T){
  			colnames(betaptsE) <-c(names,"dur")
  			colnames(betaptsR) <-c(names,"dur")
		}else{
			colnames(betaptsE) <-names
			colnames(betaptsR) <-names
		}
		
		if(durEU==T){
		  betaE <- array(0.,dim=c(NS+1,length(qtlgridSamp)) )
		  betaR <- array(0.,dim=c(NS+1,length(qtlgridSamp)) )
		}else{
		  betaE <- array(0.,dim=c(NS,length(qtlgridSamp)) )
		  betaR <- array(0.,dim=c(NS,length(qtlgridSamp)) )
		}
		#we interpolate over the coefficients so that we have the quatile regresison coefficeients defined over the whole space, not just where we estimated
		for(si in seq(1,ncol(betaptsE))){
			#first make sure that coefficients on dummies are monotone in quantile (this might not be a needed check):
			if(si<=NS){
				gpE <- c(T,betaptsE[2:length(qtlgridEst),si]>=betaptsE[1:length(qtlgridEst)-1,si]) #ensure monotonicity
				gpR <- c(T,betaptsR[2:length(qtlgridEst),si]>=betaptsR[1:length(qtlgridEst)-1,si]) #ensure monotonicity
				methodhr = "hyman"
			}else{
				#doesn't have to be monotone in duration
				gpE <- gpE==gpE
				gpR <- gpR==gpR #just setting them all to true
				methodhr = "natural"
			}
			
			#interpolate the coefficients
			if(min(qtlgridSamp)<min(qtlgridEst[gpR]) | max(qtlgridSamp)>max(qtlgridEst[gpR])){
				betaE[si,] <- spline(x=c(min(qtlgridSamp) ,qtlgridEst[gpE], max(qtlgridSamp)), 
									 y=c(min(betaptsE[gpE,si]), betaptsE[gpE,si],max(betaptsE[gpE,si])), method=methodhr, xout=qtlgridSamp)$y
			}
			if(min(qtlgridSamp)<min(qtlgridEst[gpR]) | max(qtlgridSamp)>max(qtlgridEst[gpR])){
				betaR[si,] <- spline(x=c(min(qtlgridSamp) ,qtlgridEst[gpR], max(qtlgridSamp)), 
									 y=c(min(betaptsR[gpR,si]),betaptsR[gpR,si],max(betaptsR[gpR,si])), method=methodhr, xout=qtlgridSamp)$y	
			}
			if( min(qtlgridSamp)>=min(qtlgridEst) & max(qtlgridSamp)<=max(qtlgridEst[gpR]) & max(qtlgridSamp)<=max(qtlgridEst[gpE])){
				betaE[si,] <- spline(x=qtlgridEst[gpE], y=betaptsE[gpE,si], method=methodhr, xout=qtlgridSamp)$y
				betaR[si,] <- spline(x=qtlgridEst[gpR], y=betaptsR[gpR,si], method=methodhr, xout=qtlgridSamp)$y
			}
			if( any(gpE ==F) | any(gpR==F) ){
				warning(c('non-monotone regression coefficients, ',as.character(sum(gpE==F) + sum(gpR==F)), ' time(s)') )
			}
		}

		wc_IR <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution
		wc_BR <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution
		
		qi = 1
		#sets up to create distributions of earnings change (wc) based on coefficients:
		for(q in qtlgridSamp){
			sumBetaE <- ""
			sumBetaR <- ""
			for(si in seq(1,NS) ){ #loop over the coefficients and create a string that says to sum them:
				sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
				sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
			}
			if( durEU == T){
				sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
				sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
			}else{
				sumBetaE <- paste(sumBetaE,"0" )
				sumBetaR <- paste(sumBetaR,"0" )
			}
			
			#creates the simulation-based distribution:      sample group       evaluate command that sums estimated coefficients
			wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[ sampR[,qi] , eval(parse(text=sumBetaE))]
			wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(text=sumBetaR))]

			qi = qi+1
		}
		#clean up the space:
		wc_IR_pctile[,simiter] <- quantile(wc_IR,probs=qtlgridOutTruncCor,na.rm=T)
		wc_BR_pctile[,simiter] <- quantile(wc_BR,probs=qtlgridOutTruncCor,na.rm=T)
		#est Pareto tails
		#wc_IR_pctile[1,simiter]                  <- ptail(qtlgridOut[2:(length(qtlgridOut)-1)] ,wc_IR_pctile[2:(length(qtlgridOut)-1),simiter],0.04,0.01,"L")
		#wc_IR_pctile[length(qtlgridOut),simiter] <- ptail(qtlgridOut[2:(length(qtlgridOut)-1)] ,wc_IR_pctile[2:(length(qtlgridOut)-1),simiter],0.96,0.99,"R")
		
		wc_IR_moments[,simiter]<- moments_compute_qtls(qtlpts=qtlgridOut , distpts=wc_IR_pctile[,simiter])
		wc_IR_moments[1,simiter] <- mean(wc_IR,na.rm = T)
		wc_BR_moments[,simiter]<- moments_compute_qtls(qtlpts=qtlgridOut , distpts=wc_BR_pctile[,simiter])
		wc_BR_moments[1,simiter] <- mean(wc_BR)

		rm(wc_IR)
		rm(wc_BR)

		qi=1
		wc_rec <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the recession distribution 
		for(q in qtlgridSamp){
			sumBetaR <- ""
			for(si in seq(1,NS) ){
				sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
			}
			if( durEU == T){
				sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
			}else{
				sumBetaR <- paste(sumBetaR,"0")
			}

			wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], eval(parse(text=sumBetaR)) ] 
		
			qi = qi+1
		}
		#clean up the space:
		wc_rec_pctile[,simiter] <- quantile(wc_rec,probs=qtlgridOutTruncCor,na.rm=T)
		wc_rec_moments[ ,simiter] <- moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_rec_pctile[,simiter])
		wc_rec_moments[1,simiter] <- mean( wc_rec,na.rm = T )
		rm(wc_rec)
		
		qi=1
		wc_exp <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1) #storing the expansion distribution 
		for(q in qtlgridSamp){
			sumBetaE <- ""
			for(si in seq(1,NS) ){
				sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
			}
			if( durEU == T){
				sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
			}else{
				sumBetaE <- paste(sumBetaE,"0")
			}
			wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(text=sumBetaE))]
			qi = qi+1
		}
		#clean up the space:
		wc_exp_pctile[,simiter] <- quantile(wc_exp,probs=qtlgridOutTruncCor,na.rm=T)
		wc_exp_moments[,simiter] <-moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_exp_pctile[,simiter])
		wc_exp_moments[1,simiter] <- mean( wc_exp,na.rm = T )
		rm(wc_exp)
		
		#some additional counter-factuals
		# try replacing stayers:
		qi=1
		wc_BR_exstay <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1)
		idxstay <- grepl("stay", colnames(betaptsR)) #picks up any "stay"
		for(q in qtlgridSamp){
			sumBetaR <- ""
			for(si in seq(1,NS) ){
				if(idxstay[si] ==F){
					sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
				}
				else{# use the expansions coefficients
					sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
				}
			}
			if( durEU == T){
				sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
			}else{
				sumBetaR <- paste(sumBetaR,"0")
			}
			wc_BR_exstay[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(text=sumBetaR))]
			qi = qi+1
		}
		wc_BR_exstay_pctile[,simiter] <- quantile(wc_BR_exstay,probs=qtlgridOutTruncCor,na.rm=T)
		wc_BR_exstay_moments[ ,simiter] <-moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_BR_exstay_pctile[,simiter])
		wc_BR_exstay_moments[1,simiter] <- mean( wc_BR_exstay,na.rm = T )
		
		#turn off each individually
		wc_BR_exi <-matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1)
		NNS = ifelse(durEU==T,NS+1,NS)
		for(ni in seq(1,NNS)){
			qi = 1
			for(q in qtlgridSamp){
				sumBetaR <- ""
				for(si in seq(1,NS) ){
					if(si == ni){
						sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
					else{# use the expansions coefficients
						sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
				}
				if( durEU == T){
					if(si==NS+1){
						sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
					}else{
						sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
					}
				}else{
					sumBetaR <- paste(sumBetaR,"0")
				}
				wc_BR_exi[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(text=sumBetaR))]
				qi = qi+1
			}
			wc_BR_exi_pctile[,ni,simiter] <- quantile(wc_BR_exi,probs=qtlgridOutTruncCor,na.rm=T)
			wc_BR_exi_moments[ ,ni,simiter] <- moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_BR_exi_pctile[,ni,simiter])
			wc_BR_exi_moments[1,ni,simiter] <- mean( wc_BR_exi,na.rm = T )
		}
		#turn on each individually
		wc_BR_onlyi <-matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1)
		NNS = ifelse(durEU==T,NS+1, NS)
		for(ni in seq(1,NNS)){
			qi = 1
			for(q in qtlgridSamp){
				sumBetaE <- ""
				for(si in seq(1,NS) ){
					if(si == ni){
						sumBetaE <- paste(sumBetaE,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
					else{# use the expansions coefficients
						sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
				}
				if( durEU == T){
					if(ni==NS+1){
						sumBetaE <- paste(sumBetaE,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
					}else{
						sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
					}
				}else{
					sumBetaE <- paste(sumBetaE,"0")
				}
				wc_BR_onlyi[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(text=sumBetaE))]
				qi = qi+1
			}
			wc_BR_onlyi_pctile[,ni,simiter] <- quantile(wc_BR_onlyi,probs=qtlgridOutTruncCor,na.rm=T)
			wc_BR_onlyi_moments[ ,ni,simiter] <- moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_BR_onlyi_pctile[,ni,simiter])
			wc_BR_onlyi_moments[1,ni,simiter] <- mean( wc_BR_onlyi,na.rm = T )
		}
		
		if(no_occ ==F){
			#contribution of occupation switchers
			idxsw <- grepl("_sw", colnames(betaptsR))
			wc_IR_sw <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only switch
			wc_BR_sw <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only switch
			wc_IR_sw_moments <- matrix(NA,Nmoments,Nsims)
			qi=1
			for(q in qtlgridSamp){
				sumBetaE <- ""
				sumBetaR <- ""
				for(si in seq(1,NS) ){
					if(idxsw[si] ==F){ #compare to the rest of the counter-factuals
						sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
						sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
					else{# holding the switching terms as their base
						sumBetaE <- paste(sumBetaE,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
						sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
					}
				}
				if( durEU == T){
					sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
					sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(NS+1),",qi]*dur" ))
				}else{
					sumBetaE <- paste(sumBetaE,"0")
					sumBetaR <- paste(sumBetaR,"0")
				}
				
				wc_IR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],eval(parse(text=sumBetaE))] 
				wc_BR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcExp[sampR[,qi],eval(parse(text=sumBetaR))] 
				qi = qi+1
			}
			#clean up the space:
			wc_IR_sw_pctile[,simiter] <- quantile(wc_IR_sw,probs=qtlgridOutTruncCor,na.rm=T)
			wc_BR_sw_pctile[,simiter] <- quantile(wc_BR_sw,probs=qtlgridOutTruncCor,na.rm=T)
			wc_IR_sw_moments[,simiter] <-moments_compute_qtls(qtlpts=qtlgridOut,distpts=wc_IR_sw_pctile[,simiter])
			wc_IR_sw_moments[1,simiter] <- mean( wc_IR_sw,na.rm = T )
			rm(wc_IR_sw)
			
		}			
		wc_IR_un <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
		wc_BR_un <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
		idxun <- grepl("EU_|UE_", colnames(betaptsR))
		qi=1
		for(q in qtlgridSamp){
			sumBetaE <- ""
			sumBetaR <- ""
			for(si in seq(1,NS) ){
				if(idxun[si] ==F){
					sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
					sumBetaR <- paste(sumBetaR,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
				}
				else{# use the expansions coefficients
					sumBetaE <- paste(sumBetaE,paste0("betaR[",as.character(si),",qi]*s",as.character(si)," + ") )
					sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(si),",qi]*s",as.character(si)," + ") )
				}
			}
			if( durEU == T){
				sumBetaE <- paste(sumBetaE,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
				sumBetaR <- paste(sumBetaR,paste0("betaE[",as.character(NS+1),",qi]*dur" ))
			}else{
				sumBetaE <- paste(sumBetaE,"0")
				sumBetaR <- paste(sumBetaR,"0")
			}
			wc_IR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],eval(parse(text=sumBetaE))] 
			wc_BR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcExp[sampR[,qi],eval(parse(text=sumBetaR))] 
			qi = qi+1
		}
		#clean up the space:
		wc_IR_un_pctile[,simiter] <- quantile(wc_IR_un,probs=qtlgridOutTruncCor,na.rm=T)
		wc_BR_un_pctile[,simiter] <- quantile(wc_BR_un,probs=qtlgridOutTruncCor,na.rm=T)
		rm(wc_IR_un)

	}#simiter	
	
	if(no_occ==F){	
		return(list(betaptsE = betaptsE,betaptsR=betaptsR,wc_IR = wc_IR_pctile, wc_IR_sw= wc_IR_sw_pctile, wc_IR_un= wc_IR_un_pctile,wc_BR_un= wc_BR_un_pctile,
				wc_rec = wc_rec_pctile,wc_exp = wc_exp_pctile, wc_BR = wc_BR_pctile,wc_BR_mmts = wc_BR_moments,wc_IR_mmts = wc_IR_moments, 
				wc_exp_mmts = wc_exp_moments, wc_rec_mmts = wc_rec_moments, wc_IR_sw_mmts = wc_IR_sw_moments,wc_BR_sw = wc_BR_sw_pctile, wc_BR_exi = wc_BR_exi_pctile,
				wc_BR_exi_mmts = wc_BR_exi_moments,wc_BR_onlyi = wc_BR_onlyi_pctile))
	}else{
		return(list(betaptsE = betaptsE,betaptsR=betaptsR,wc_IR = wc_IR_pctile, wc_IR_un= wc_IR_un_pctile,wc_BR_un= wc_BR_un_pctile,
					wc_rec = wc_rec_pctile,wc_exp = wc_exp_pctile, wc_BR = wc_BR_pctile, wc_BR_exi = wc_BR_exi_pctile,wc_BR_onlyi = wc_BR_onlyi_pctile,
					wc_BR_mmts = wc_BR_moments,wc_IR_mmts = wc_IR_moments,wc_exp_mmts = wc_exp_moments, wc_rec_mmts = wc_rec_moments,
					wc_BR_exi_mmts = wc_BR_exi_moments, wc_BR_onlyi_mmts = wc_BR_onlyi_moments))
	}
}
ptail <- function(qtlpts,distpts,Mqtl,Tqtl,LR = "R"){
	if(missing(LR)){
		LR = "R"
	}
	dist_fun <- approxfun(qtlpts,distpts)
	if(LR=="R"){
		alphahat = (Tqtl - Mqtl)/( (1-Tqtl)*log(abs(dist_fun(Tqtl))) + sum( 0.01*log(abs(dist_fun(seq(Mqtl+0.01,Tqtl,0.01)))) )- (1-Mqtl)*log(abs(dist_fun(Mqtl))) )
	}else{
		alphahat = (Mqtl - Tqtl)/( Tqtl*log(abs(dist_fun(Tqtl)))     + sum( 0.01*log(abs(dist_fun(seq(Tqtl,Mqtl-0.01,0.01)))) ) -    Mqtl*log(abs(dist_fun(Mqtl))) )
	}
	return(alphahat)
}

moments_compute_qtls <- function(qtlpts, distpts){
	# take the distributions points evaluated at the quantile points and crank out some moments
	npts <- length(qtlpts)
	dist_fun <- splinefun(qtlpts, distpts)
	qtl_grid <- c(0.5*qtlpts[1] , 0.5*(qtlpts[seq(2,npts)]+qtlpts[seq(1,npts-1)]), 0.5*qtlpts[npts]+0.5)
	qtl_step <- qtl_grid[seq(2,npts+1)]-qtl_grid[seq(1,npts)]
	mean_hr <- qtl_step%*%distpts
	median_hr <- dist_fun(0.5)
	mad_hr <- wtd.quantile( x=abs(distpts- median_hr),probs = 0.5) #need to figure out how weights work here
	GroenMeed_hr <- (mean_hr-median_hr)/(sum(abs(distpts - median_hr)*qtl_step))
	Moors_hr <- (dist_fun(7/8)-dist_fun(5/8)+dist_fun(3/8)-dist_fun(1/8))/(dist_fun(3/4)-dist_fun(1/4))
	return( c(mean_hr,median_hr,mad_hr,GroenMeed_hr,Moors_hr))
}

########################################################################
########################################################################
## Seams version -------------------------------

DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))
#create some vars to handle net flows
DTseam[ , last.occ :=shift(occ), by=id]
DTseam[ !is.finite(occL), occL := occ]
DTseam[ EE_wave==T , occL := last.occ]
DTseam[ , next.occ :=shift(occ,type="lead"), by=id]
DTseam[ !is.finite(occD), occD := next.occ]
#++++++++++++++++++++
# export to stata:
#++++++++++++++++++++
export_loc = "~/Dropbox/Carrillo_Visschers_Wiczer/SIPP/sipp_wavefreq.dta"
write.dta(subset(DTseam,select=c("ageGrp","date","EE_wave","EU_wave","HSCol","id","lfstat_wave","max.unempdur_wave",
 	"occ_1d","occ_wave","occ90","occD","occL","rawwgchange_wave","recIndic_wave","recIndic2_wave","switched_wave",
 	"UE_wave","unrt","wagechange_anan","wagechange_wave","wave","wavewage","wpfinwgt","perwt","truncweight")),
 	export_loc)
#
toKeep <- c("switchedOcc_wave","switched_wave","switched_anan","esr_max",
				 "ageGrp","HSCol","next.stable_emp","last.stable_emp","wavewage",
				 "recIndic","recIndic_wave","recIndic2_wave","recIndic_stint","levwage","max.unempdur_wave","max.unempdur",
				 "wagechange_wave","wagechangeEUE_wave","rawwgchange_wave","rawwgchangeEUE_wave","wagechange_anan",
				 "wagechange_notransbad","wagechange_wave_low","wagechange_wave_high","wagechange_wave_jcbad","pctmaxmis",
				 "EE_wave","EU_wave","UE_wave","changer","stayer","EE_anan","EU_anan","UE_anan","changer_anan","stayer_anan","valid_anan",
				 "unrt","wpfinwgt","perwt","truncweight","cleaningtruncweight","lastann.wavewage","matched_EUUE_anan",
				 "occwage_wave","occwagechange_wave","occwagechangeEUE_wave","occ_wave","occL", "occD",
				 "lfstat_wave","next.lfstat_wave","wave","id","date","panel")
# select toKeep columns only
DTseam <- DTseam[, toKeep, with = FALSE]

DTseam[ , last2.stable_emp:= shift(last.stable_emp), by=id]
DTseam[ , last3.stable_emp:= shift(last2.stable_emp), by=id]
DTseam[ wave-2 != shift(wave,2), last2.stable_emp := NA]
DTseam[ wave-3 != shift(wave,3), last3.stable_emp := NA]

DTseam[ EU_anan& changer_anan & is.na(matched_EUUE_anan),changer_anan:= NA]
DTseam[ UE_anan& changer_anan & is.na(matched_EUUE_anan),changer_anan:= NA]
DTseam[ EU_anan& changer_anan & is.na(matched_EUUE_anan),EU_anan:= NA]
DTseam[ UE_anan& changer_anan & is.na(matched_EUUE_anan),UE_anan:= NA]

DTseam[ , nextann.wavewage:=NULL]
DTseam[ , nextann.wavewage := levwage + shift(levwage,type="lead") + shift(levwage,2,type="lead")  , by=id]
DTseam[ wave+1!=shift(wave,type = "lead") | wave+2!=shift(wave,2,type = "lead") , nextann.wavewage := NA]
DTseam[ , nextann.wavewage := log(nextann.wavewage + (1+nextann.wavewage^2)^.5) ]


if(wc == "wagechangeEUE_wave"|wc == "reswgchangeEUE_wave"){
	DTseam[ ,wtEUE:= eval(as.name(wt))]
	DTseam[ UE_wave==T,wtEUE:= 0.]
	DTseam[ EU_wave==T,wtEUE:= eval(as.name(wt))]
	origwt = wt
	wt = "wtEUE"
}else{
	wt = "truncweight"
	#upweight the EE_wave to EE_anan and EU_wave and UE_wave to EU_anan and UE_anan
	DTseam[ , wtanan:= eval(as.name(wt))]
	
	#could upweight the EE_wave to EE_anan and EU_wave and UE_wave to EU_anan and UE_anan
	fracEE = DTseam[ changer_anan & !(EU_wave|UE_wave|EE_wave) , mean(EE_anan)]  #DTseam[ EE_anan==T & !(EU_anan|UE_anan)& is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]/DTseam[ EE_wave==T & is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]
	fracEU = DTseam[ changer_anan & !(EU_wave|UE_wave|EE_wave) , mean(EU_anan)]  #DTseam[ EU_anan==T &                     is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]/DTseam[ EU_wave==T & is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]
	fracUE = DTseam[ changer_anan & !(EU_wave|UE_wave|EE_wave) , mean(UE_anan)]  #DTseam[ UE_anan==T &                     is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]/DTseam[ UE_wave==T & is.finite(wagechange_anan) & is.finite(switched_wave), sum(eval(as.name(wt)))]
	fracEUUE = DTseam[ changer_anan & !(EU_wave|UE_wave|EE_wave) , mean(UE_anan|EU_anan)] 
	fracEE = fracEE/(fracEE+fracEUUE)
	fracEUUE = fracEUUE/(fracEE+fracEUUE)
	
	scaleEE = fracEE*DTseam[ changer_anan==T ,sum(eval(as.name(wt))) ]/DTseam[ changer_anan==T & EE_wave==T ,sum(eval(as.name(wt))) ]
	scaleEU = fracEU*DTseam[ changer_anan==T ,sum(eval(as.name(wt))) ]/DTseam[ changer_anan==T & EU_wave==T ,sum(eval(as.name(wt))) ]
	scaleUE = fracUE*DTseam[ changer_anan==T ,sum(eval(as.name(wt))) ]/DTseam[ changer_anan==T & UE_wave==T ,sum(eval(as.name(wt))) ]
	scaleEUUE = fracEUUE*DTseam[ changer_anan==T ,sum(eval(as.name(wt))) ]/DTseam[ changer_anan==T & (EU_wave|UE_wave) ,sum(eval(as.name(wt))) ]
	
	DTseam[ EE_wave==T, wtanan:= eval(as.name(wt))*scaleEE]
	DTseam[ EU_wave==T, wtanan:= eval(as.name(wt))*scaleEUUE]
	DTseam[ UE_wave==T, wtanan:= eval(as.name(wt))*scaleEUUE]
	origwt = wt
	wt = "wtanan"
}

# setup labels
if(wc == "wagechangeEUE_wave"){
	wclab = "reswgEUE"
}else if( wc == "wagechange_wave"){
	wclab = "reswg"
}else if(wc == "rawwgchangeEUE_wave"){
	wclab = "rawwgEUE"
}else if( wc == "rawwgchange_wave"){
	wclab = "rawwg"
}else if( wc == "wagechange_anan" ){
	wclab = "reswgANAN"
}

if(recDef == "recIndic_wave"){
	reclab = "NBER"
}else if( recDef == "recIndic2_wave"){
	reclab = "urt"
}else if(recDef == "recIndic_stint"){
	reclab = "recstint"
}
if( wc == "wagechange_anan"|wc == "rawwgchange_anan"){
	NS=8
	freq = "annual"
}else{
	NS=6
	freq = "wave"	
}


if(freq == "wave"){
	DTseam[, ch := changer]
	DTseam[, st := stayer]
	DTseam[, sw := switched_wave]
	DTseam[, EUfrq := EU_wave]
	DTseam[, EEfrq := EE_wave]
	DTseam[, UEfrq := UE_wave]
	DTseam <- DTseam[ (sw|!sw) & (ch|st), ]
}else{
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0),ch := changer_anan]
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0),st := stayer_anan]
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0),sw := switched_wave]
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0), EUfrq := EU_wave]
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0), EEfrq := EE_wave]
	DTseam[lastann.wavewage>minLEarn & is.finite(lastann.wavewage) & (EU_wave==T|nextann.wavewage>0), UEfrq := UE_wave]
	DTseam[!(EUfrq|EEfrq|UEfrq) & ch==T, ch := NA]
	DTseam <- DTseam[ valid_anan==T, ]
}
DTseam[ EUfrq==T|UEfrq==T, dur := max.unempdur_wave]
DTseam[ !(EUfrq==T|UEfrq==T) | !is.finite(dur), dur := 0.]
#++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------- The Decomposition ----------------
#++++++++++++++++++++++++++++++++++++++++++++++++++++
if( freq == "wave"){
	if(wdur ==T){
		MM_betaE_betaR_IR <- MMdecomp(DTseam,NS,recDef,wcname=wc,wtname=wt,std_errs = MMstd_errs, no_occ = F,durEU = T)
		saveRDS(MM_betaE_betaR_IR,paste0(outputdir,"/MM_waveallEUE.RData"))
	}else{
		MM_betaE_betaR_IR <- MMdecomp(DTseam,NS,recDef,wcname=wc,wtname=wt,std_errs = MMstd_errs, no_occ = F,durEU = F)
		saveRDS(MM_betaE_betaR_IR,paste0(outputdir,"/MM_waveallEUE_nodur.RData"))
	}
}else{
	if(wdur==T){
		MM_betaE_betaR_IR <- MMdecomp(DTseam,NS,recDef,wcname=wc,wtname=wt,std_errs = MMstd_errs, no_occ = F,durEU = T)
		saveRDS(MM_betaE_betaR_IR,paste0(outputdir,"/MM_ANAN_dur.RData"))
	}else{
		MM_betaE_betaR_IR <- MMdecomp(DTseam,NS,recDef,wcname=wc,wtname=wt,std_errs = MMstd_errs, no_occ = F,durEU = F)
		saveRDS(MM_betaE_betaR_IR,paste0(outputdir,"/MM_ANAN.RData"))
	}
	# MM_betaE_betaR_IR <-readRDS(paste0(outputdir,"/MM_ANAN.RData"))
}

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# The calibration targets ------------------------------------------------
#J2J
DTseam[ lfstat_wave ==1  & (st|ch) , wtd.mean(EE_wave,na.rm = T, weights = truncweight )]
#find rate
DTseam[ lfstat_wave >1  , wtd.mean(UE_wave,na.rm = T, weights = truncweight )]
#sep rate
DTseam[ lfstat_wave ==1  & (st|ch), wtd.mean(EU_wave,na.rm = T, weights = truncweight )]
#sw_U rate
DTseam[ changer_anan & EU_wave==T , wtd.mean(switched_wave,na.rm = T, weights = truncweight )]
#sw_E rate
DTseam[ changer_anan & EE_wave==T , wtd.mean(switched_wave,na.rm = T, weights = truncweight )]
#sw_st rate
DTseam[ stayer_anan==T  , wtd.mean(switched_wave,na.rm = T, weights = truncweight )]

DTseam[ switched_wave==T & changer_anan & UE_wave==T, wtd.mean(max.unempdur_wave, na.rm=T, weights=truncweight)]

DTseam[ switched_wave==F & changer_anan & UE_wave==T, wtd.mean(max.unempdur_wave, na.rm=T, weights=truncweight)]

DTseam[ lfstat_wave ==1 & is.finite(occ_wave), wtd.mean(occwage_wave,na.rm=T, weights=truncweight)/2 , by=occ_wave] - max(DTseam[ lfstat_wave ==1  & is.finite(occ_wave), wtd.mean(occwage_wave,na.rm=T, weights=truncweight)/2 , by=occ_wave])

#earnings change distributions:
DTseam[ stayer_anan==T & switched_wave==T, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]
DTseam[ stayer_anan==T & switched_wave==F, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]

DTseam[ changer_anan==T & EU_wave==T & switched_wave==T, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]
DTseam[ changer_anan==T & EU_wave==T & switched_wave==F, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]

DTseam[ changer_anan==T & UE_wave==T & switched_wave==T, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]
DTseam[ changer_anan==T & UE_wave==T & switched_wave==F, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]

DTseam[ changer_anan==T & EE_wave==T & switched_wave==T, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]
DTseam[ changer_anan==T & EE_wave==T & switched_wave==F, wtd.quantile(wagechange_anan, probs = c(.1,.25,.5,.75,.9),na.rm = T, weights = truncweight )]

#cyclicality of flows
DTseam[ lfstat_wave == 1 & recIndic_wave==F & (st|ch), wtd.mean(EE_wave, na.rm=T, weights=truncweight)]/
DTseam[ lfstat_wave == 1 & recIndic_wave==T & (st|ch), wtd.mean(EE_wave, na.rm=T, weights=truncweight)]

DTseam[ lfstat_wave > 1  & recIndic_wave==F , wtd.mean(UE_wave, na.rm=T, weights=truncweight)]/
DTseam[ lfstat_wave > 1  & recIndic_wave==T , wtd.mean(UE_wave, na.rm=T, weights=truncweight)]

DTseam[ lfstat_wave == 1 & recIndic_wave==F & (st|ch), wtd.mean(EU_wave, na.rm=T, weights=truncweight)]/
DTseam[ lfstat_wave == 1 & recIndic_wave==T & (st|ch), wtd.mean(EU_wave, na.rm=T, weights=truncweight)]

DTseam[ EE_wave == 1 & recIndic_wave==F & (st|ch), wtd.mean(switched_anan, na.rm=T, weights=truncweight)]/
	DTseam[ EE_wave == 1 & recIndic_wave==T & (st|ch), wtd.mean(switched_anan, na.rm=T, weights=truncweight)]

DTseam[ UE_wave== 1  & recIndic_wave==F, wtd.mean(switched_anan, na.rm=T, weights=truncweight)]/
	DTseam[ UE_wave== 1  & recIndic_wave==T, wtd.mean(switched_anan, na.rm=T, weights=truncweight)]

DTseam[ EU_wave == 1 & recIndic_wave==F & (st|ch), wtd.mean(switched_anan, na.rm=T, weights=truncweight)]/
	DTseam[ EU_wave == 1 & recIndic_wave==T & (st|ch), wtd.mean(switched_anan, na.rm=T, weights=truncweight)]

#cyclicality of returns (old measure)
MVswrec1_qtls =  DTseam[ switched_wave==T & ch==T & recIndic_wave==T, wtd.quantile(wagechange_anan, probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=T,weights=truncweight) ]
MVswrec0_qtls =  DTseam[ switched_wave==T & ch==T & recIndic_wave==F, wtd.quantile(wagechange_anan, probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=T,weights=truncweight) ]

MVnsrec1_qtls =  DTseam[ switched_wave==F & ch==T & recIndic_wave==T, wtd.quantile(wagechange_anan, probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=T,weights=truncweight) ]
MVnsrec0_qtls =  DTseam[ switched_wave==F & ch==T & recIndic_wave==F, wtd.quantile(wagechange_anan, probs=c(0.1,0.25,0.5,0.75,0.9), na.rm=T,weights=truncweight) ]

MVswrec1_qtls/MVswrec0_qtls

MVnsrec1_qtls/MVnsrec0_qtls

#cyclicality of returns (new measure)
MVswrec1_qtls =  DTseam[ switched_wave==T & ch==T & recIndic_wave==T, wtd.quantile(wagechange_anan, probs=c(0.025,0.05,0.5,0.95,0.975), na.rm=T,weights=truncweight) ]
MVswrec0_qtls =  DTseam[ switched_wave==T & ch==T & recIndic_wave==F, wtd.quantile(wagechange_anan, probs=c(0.025,0.05,0.5,0.95,0.975), na.rm=T,weights=truncweight) ]

MVnsrec1_qtls =  DTseam[ switched_wave==F & ch==T & recIndic_wave==T, wtd.quantile(wagechange_anan, probs=c(0.025,0.05,0.5,0.95,0.975), na.rm=T,weights=truncweight) ]
MVnsrec0_qtls =  DTseam[ switched_wave==F & ch==T & recIndic_wave==F, wtd.quantile(wagechange_anan, probs=c(0.025,0.05,0.5,0.95,0.975), na.rm=T,weights=truncweight) ]

MVswrec1_qtls- MVswrec0_qtls

MVnsrec1_qtls- MVnsrec0_qtls


# overall net flow matrix
DTseam[ sw==T & ch==T & occL ==occD, occLDNA := T ]
DTseam[ sw==T & ch==T & occLDNA ==T , occL := NA ]
DTseam[ sw==T & ch==T & occLDNA ==T , occD := NA ]
DTseam[ , occLDNA := NULL ]

grossL <- array(0,dim = c(4,4)) 
noccsw <- array(0,dim = c(4)) 
noccsw <- DTseam[ sw==T & ch==T , ftable(occL)/sum(is.finite(occL)& is.finite(occD))]
grossL[1,c(2,3,4)] <- DTseam[ occL==1 & sw==T & ch==T , ftable(occD)/sum(is.finite(occD))]*noccsw[1]
grossL[2,c(1,3,4)] <- DTseam[ occL==2 & sw==T & ch==T , ftable(occD)/sum(is.finite(occD))]*noccsw[2]
grossL[3,c(1,2,4)] <- DTseam[ occL==3 & sw==T & ch==T , ftable(occD)/sum(is.finite(occD))]*noccsw[3]
grossL[4,c(1,2,3)] <- DTseam[ occL==4 & sw==T & ch==T , ftable(occD)/sum(is.finite(occD))]*noccsw[4]

netL <- array(0, dim=c(4,4))
netLUi <- array(0, dim=c(4,4,2))
grossVar <- array(dim=c(5,2))
for( li in seq(1,4)){
	for( di in seq(1,4)){
		netL[li,di] = grossL[li,di] - grossL[di,li]
	}
}

gross_marg             <- DTseam[ sw==T & ch==T , ftable(occD)/sum(is.finite(occD))]


for(Ui in c(T,F)){
	ii = ifelse(Ui, 1,2)
	grossL <- array(0,dim = c(4,4)) 
	noccsw <- array(0,dim = c(4)) 
	noccsw <- DTseam[ sw==T & ch==T & EUfrq==Ui, ftable(occL)/sum(is.finite(occL)& is.finite(occD))]
	grossL[1,c(2,3,4)] <- DTseam[ occL==1 & sw==T & ch==T & EUfrq==Ui , ftable(occD)/sum(is.finite(occD))]*noccsw[1]
	grossL[2,c(1,3,4)] <- DTseam[ occL==2 & sw==T & ch==T & EUfrq==Ui , ftable(occD)/sum(is.finite(occD))]*noccsw[2]
	grossL[3,c(1,2,4)] <- DTseam[ occL==3 & sw==T & ch==T & EUfrq==Ui , ftable(occD)/sum(is.finite(occD))]*noccsw[3]
	grossL[4,c(1,2,3)] <- DTseam[ occL==4 & sw==T & ch==T & EUfrq==Ui , ftable(occD)/sum(is.finite(occD))]*noccsw[4]
	grossVar[1,ii] = var(grossL[1,c(2,3,4)]/noccsw[1])
	grossVar[2,ii] = var(grossL[2,c(1,3,4)]/noccsw[2])
	grossVar[3,ii] = var(grossL[3,c(1,2,4)]/noccsw[3])
	grossVar[4,ii] = var(grossL[4,c(1,2,3)]/noccsw[4])
	# weighted average across these, the calibration target!
	grossVar[5,ii] = grossVar[1,ii]*noccsw[1]+grossVar[2,ii]*noccsw[2]+grossVar[3,ii]*noccsw[3]+grossVar[4,ii]*noccsw[4]
	
	netL <- array(0, dim=c(4,4))
	for( li in seq(1,4)){
		for( di in seq(1,4)){
			netLUi[li,di,ii] = grossL[li,di] - grossL[di,li]
		}
	}
}
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# End moments
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#Extracting the dist and moments -----------------------------------
 # wcExp <- subset(DTseam,eval(as.name(recDef))==F)
 # wcRec <- subset(DTseam,eval(as.name(recDef))==T)
 # dist_exp      <- wcExp[ , wtd.quantile(eval(as.name(wc)),probs=qtlgridOut,weights=eval(as.name(wt)), na.rm=T)]
 # dist_rec      <- wcRec[ , wtd.quantile(eval(as.name(wc)),probs=qtlgridOut,weights=eval(as.name(wt)), na.rm=T)]
# wc_Exp_mmts   <- wcExp[  , c(wtd.mean(eval(as.name(wc)),weights=eval(as.name(wt)),na.rm=T), 
# 							 wtd.4qtlmoments( xt= eval(as.name(wc)) ,wt= eval(as.name(wt))))]
# wc_Rec_mmts   <- wcRec[  , c(wtd.mean(eval(as.name(wc)),weights=eval(as.name(wt)),na.rm=T), 
# 							 wtd.4qtlmoments( xt= eval(as.name(wc)) ,wt= eval(as.name(wt))))]

Dist_dif = data.table(cbind(qtlgridOut,MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_rec,MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR))#,
							#MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR_un,MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR_sw))
names(Dist_dif) <- c("Qtls", "Dif", "CF Dif")#,"CF Dif, ex un","CF Dif, ex sw")
Dist_melt <- melt(Dist_dif,id.vars = "Qtls")
ggplot(Dist_melt, aes(x=Qtls,y=value,color=variable))+geom_smooth(se=F,span=0.1)+theme_bw()+
	xlab("Annual-Annual Log Earnings Change Quantile")+ylab("Difference Expansion - Recession")+
	scale_color_manual(labels=c("Data","Recession Returns"),#,"Recession Returns, Expansion Unemployment Returns","Recession Returns, Expansion Switching Returns"),
					   values=c("black",hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:3)]))+ 
	scale_x_continuous(breaks=seq(0,1,1/4),minor_breaks = seq(0,1,1/8) )+
	theme(legend.title = element_blank(),text = element_text(size=20),
		  legend.position = c(0.5,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))	
if(wdur==T){
	nametab = "cf_qtls_wdur"
}else{
	nametab = "cf_qtls"
}
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".eps"),device=cairo_ps,height=5,width=10)
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".png"),height=5,width=10)
out_wavedist <- xtable(Dist_dif, digits=3, 
					   caption=paste0("Decomposition of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
print(out_wavedist,include.rownames=T, include.colnames=T,
	  file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
Dist_dif_qtl <- Dist_dif 

#with levels on the x-axis
Dist_dif = data.table(cbind(MM_betaE_betaR_IR$wc_exp*.75+MM_betaE_betaR_IR$wc_rec*.25,
							MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_rec,MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR))#,
						#	MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR_un,MM_betaE_betaR_IR$wc_exp-MM_betaE_betaR_IR$wc_BR_sw))
names(Dist_dif) <- c("WChng", "Dif", "CF Dif")#,"CF Dif, ex un","CF Dif, ex sw")
Dist_melt <- melt(Dist_dif,id.vars = "WChng")
ggplot(Dist_melt, aes(x=WChng,y=value,color=variable))+geom_smooth(se=F,span=0.1,method="loess",n=20)+theme_bw()+
	xlab("Annual-Annual Log Earnings Change")+ylab("Difference Expansion - Recession")+
	scale_color_manual(labels=c("Data","Recession Returns"),#,"Recession Returns, Expansion Unemployment Returns","Recession Returns, Expansion Switching Returns"),
					   values=c("black",hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:3)]))+ 
	theme(legend.title = element_blank(),text = element_text(size=20),
		  legend.position = c(0.65,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))	
if(wdur==T){
	nametab = "cf_wchng_wdur"
}else{
	nametab = "cf_wchng"
}
nametab = "cf_wchng_wdur"
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".eps"),device=cairo_ps,height=5,width=10)
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".png"),height=5,width=10)
out_wavedist <- xtable(Dist_dif, digits=3, 
					   caption=paste0("Decomposition of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
print(out_wavedist,include.rownames=T, include.colnames=T,
	  file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))

#figure out what does what: decomp top/bottom 5%, 10%:
if(wdur==T){
	NNS=NS+1
}else{
	NNS= NS
}
Decomp_tails = array(0.,dim=c(4,NNS+1))
Nqtl_5 = sum(qtlgridOut<=0.05)
Nqtl_10 = sum(qtlgridOut<=0.025)
for(si in seq(1,NNS)){
	Decomp_tails[1, si] = sum(-0.5*(MM_betaE_betaR_IR$wc_BR[qtlgridOut<=.025] - MM_betaE_betaR_IR$wc_BR_exi[qtlgridOut<=.025,si,1]) + 0.5*(MM_betaE_betaR_IR$wc_exp[qtlgridOut<=.025] - MM_betaE_betaR_IR$wc_BR_onlyi[qtlgridOut<=.025,si,1]))/sum(MM_betaE_betaR_IR$wc_exp[qtlgridOut<=.025] - MM_betaE_betaR_IR$wc_rec[qtlgridOut<=.025])
	Decomp_tails[2, si] = sum(-0.5*(MM_betaE_betaR_IR$wc_BR[qtlgridOut<=0.05] - MM_betaE_betaR_IR$wc_BR_exi[qtlgridOut<=0.05,si,1]) + 0.5*(MM_betaE_betaR_IR$wc_exp[qtlgridOut<=0.05] - MM_betaE_betaR_IR$wc_BR_onlyi[qtlgridOut<=0.05,si,1]))/sum(MM_betaE_betaR_IR$wc_exp[qtlgridOut<=0.05] - MM_betaE_betaR_IR$wc_rec[qtlgridOut<=0.05])
	Decomp_tails[3, si] = sum(-0.5*(MM_betaE_betaR_IR$wc_BR[qtlgridOut>=0.95] - MM_betaE_betaR_IR$wc_BR_exi[qtlgridOut>=0.95,si,1]) + 0.5*(MM_betaE_betaR_IR$wc_exp[qtlgridOut>=0.95] - MM_betaE_betaR_IR$wc_BR_onlyi[qtlgridOut>=0.95,si,1]))/sum(MM_betaE_betaR_IR$wc_exp[qtlgridOut>=0.95] - MM_betaE_betaR_IR$wc_rec[qtlgridOut>=0.95])
	Decomp_tails[4, si] = sum(-0.5*(MM_betaE_betaR_IR$wc_BR[qtlgridOut>=.975] - MM_betaE_betaR_IR$wc_BR_exi[qtlgridOut>=.975,si,1]) + 0.5*(MM_betaE_betaR_IR$wc_exp[qtlgridOut>=.975] - MM_betaE_betaR_IR$wc_BR_onlyi[qtlgridOut>=.975,si,1]))/sum(MM_betaE_betaR_IR$wc_exp[qtlgridOut>=.975] - MM_betaE_betaR_IR$wc_rec[qtlgridOut>=.975])
}
Decomp_tails[, NNS+1] = rowSums(Decomp_tails)
if(wdur==T){
	Decomp_tails <- Decomp_tails[, c(seq(1,NS-2),NS+1,NS-1,NS,NNS+1)]
}

Decomp_tails<- data.table(Decomp_tails)
rownames(Decomp_tails) <- c("$\\leq$ 2.5\\%","$\\leq$ 5.0\\%","$\\geq$ 95.0\\%","$\\geq$ 97.5\\%")
if(NS==8){
	if(wdur==F){ names(Decomp_tails)<-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay_nosw","stay_sw","total")
	}else{  names(Decomp_tails)<-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","dur","stay_nosw","stay_sw","total")}
}else{
	names(Decomp_tails)<-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay_nosw","stay_sw","total")
}

if(wdur==T){
nametab <- "TailMM_wdur"
}else{
	nametab<-"tailMM"
}
aligntxt =ifelse(wdur==T, "l|cccccc|c|cc|c","l|cccccc|cc|c")

Decomp_tails <- xtable(Decomp_tails, digits=3, 
							align=aligntxt, caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
print(Decomp_tails,include.rownames=T, hline.after= c(nrow(Decomp_tails)), 
	  file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
	

#!!!!!!!!!!!
# Compare distributions ---------------------------------------------------
tabqtls <- c(.05,.10,.25,.5,.75,.90,.95)
tN <- (length(tabqtls)+1)
ann_wavedist <- array(0., dim=c(3,length(tabqtls)+2))
# ann_wavedist[1,1]   <- DTseam[!(eval(as.name(recDef)) ) & (sw|!sw) 
# 							  &(st|ch),  wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
# ann_wavedist[1,2:tN]<- DTseam[!(eval(as.name(recDef)) ) & (sw|!sw) 
# 							  &(st|ch),  wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
ann_wavedist[1,1+tN]<- DTseam[                             (sw|!sw) 
							   &(st|ch), wtd.mean(!eval(as.name(recDef)),na.rm=T,weights=eval(as.name(wt)))]
# ann_wavedist[3,1]   <- DTseam[ (eval(as.name(recDef))  ) & (sw|!sw) 
# 							   &(st|ch), wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
# ann_wavedist[3,2:tN]<- DTseam[ (eval(as.name(recDef))  ) & (sw|!sw) 
# 							   &(st|ch), wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
ann_wavedist[3,1+tN]<- DTseam[                             (sw|!sw) 
							   &(st|ch), wtd.mean(eval(as.name(recDef)),na.rm=T,weights=eval(as.name(wt)))]
ann_wavedist[2,1]   <- MM_betaE_betaR_IR$wc_BR_mmts[1]
ann_wavedist[2,2:tN]<- approx(qtlgridOut, MM_betaE_betaR_IR$wc_BR, xout=tabqtls)$y
ann_wavedist[2,1+tN]<- ann_wavedist[3,1+tN]

ann_wavedist[3,2:tN]<- approx(qtlgridOut, MM_betaE_betaR_IR$wc_rec, xout=tabqtls)$y
ann_wavedist[1,2:tN]<- approx(qtlgridOut, MM_betaE_betaR_IR$wc_exp, xout=tabqtls)$y
ann_wavedist[1,1]<- MM_betaE_betaR_IR$wc_exp_mmts[1]
ann_wavedist[3,1]<- MM_betaE_betaR_IR$wc_rec_mmts[1]


plt_wavedist <- data.table(ann_wavedist)
names(plt_wavedist) <- c("Mean","P5","P10","P25","P50","P75","P90","P95","Frac")
plt_wavedist[ , Cycle := as.factor(c("Exp","CF","Rec"))]
plt_melt <- melt(plt_wavedist, id.vars = "Cycle")

ggplot( plt_wavedist , aes(Cycle)) + theme_bw()+
	geom_boxplot( aes( ymin=P10,lower=P25,middle=Mean,upper=P75,ymax=P90 , color=Cycle, width=Frac),stat="identity",varwidth = T )+
	scale_x_discrete(labels=c("Counter-Factual","Expansion","Recession"))+ xlab("")+ylab("Log earnings change")+
	scale_color_manual(values=c("purple","blue","red"))  #+ylim(c(-0.51,0.51))
nametab = "cf_box"
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".eps"),device=cairo_ps,height=5,width=10)
ggsave(file=paste0(outputdir,"/",nametab,"_",reclab,"_",wclab,".png"),height=5,width=10)
out_wavedist <- xtable(plt_wavedist, digits=3, 
					   caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
print(out_wavedist,include.rownames=T, include.colnames=T,
	  file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))





##################################################
### plot stuff
##################################################
dt_mm <- data.table(cbind( seq(0.07,0.93,0.01),win_dif_rec_exp_sim ))
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


dt_mm <- data.table(cbind( seq(0.07,0.93,0.01),win_dif_rec_exp_sim,win_dif_IR_exp ))
names(dt_mm) <- c("Quantile","Data","Counter Factual - Rec Flows")
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

dt_mm <- data.table(cbind( seq(0.07,0.93,0.01),win_dif_rec_exp_sim,win_dif_BR_exp ))
names(dt_mm) <- c("Quantile","Data","Counter Factual - Rec Returns")
dt_mm_melted <- melt(dt_mm, id= "Quantile")
ggplot( dt_mm_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Recession - Expansion Distribution") + 
	xlab("Earnings Growth Quantile") +
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[1:2]) ) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_BR_all.eps",height=5,width=10)
ggsave("./Figures/MMwave_BR_all.png",height=5,width=10)

dt_mm_un <- data.table(cbind( seq(0.07,0.93,0.01),win_dif_rec_exp_sim,win_dif_IR_exp,win_dif_IR_un_exp ))
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

#some coefficients:
mmtaballqtls <- qtlgridEst
dt_co_stay <- data.table(cbind( qtlgridEst,MM_waveall_betaE_betaR_IR$betaptsE[,7], MM_waveall_betaE_betaR_IR$betaptsR[,7]))
names(dt_co_stay) <- c("Quantile","Expansion","Recession")
dt_co_stay_melted<-melt(dt_co_stay,id="Quantile")
ggplot( dt_co_stay_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Coefficient Value") + 
	xlab("Earnings Growth Quantile") +
	scale_color_manual(values = c("blue","red")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.6,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_costay.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_costay.png",height=5,width=10)

dt_co_EEnos <- data.table(cbind( qtlgridEst,MM_waveall_betaE_betaR_IR$betaptsE[,4], MM_waveall_betaE_betaR_IR$betaptsR[,4]))
names(dt_co_EEnos) <- c("Quantile","Expansion","Recession")
dt_co_EEnos_melted<-melt(dt_co_EEnos,id="Quantile")
ggplot( dt_co_EEnos_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Coefficient Value") + 
	xlab("Earnings Growth Quantile") +ylim(c(-2.5,2.5))+
	scale_color_manual(values = c("blue","red")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.6,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_coEEnosw.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_coEEnosw.png",height=5,width=10)

dt_co_EEs <- data.table(cbind( qtlgridEst,MM_waveall_betaE_betaR_IR$betaptsE[,1], MM_waveall_betaE_betaR_IR$betaptsR[,1]))
names(dt_co_EEs) <- c("Quantile","Expansion","Recession")
dt_co_EEs_melted<-melt(dt_co_EEs,id="Quantile")
ggplot( dt_co_EEs_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Coefficient Value") + 
	xlab("Earnings Growth Quantile")  +ylim(c(-2.5,2.5))+
	scale_color_manual(values = c("blue","red")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.6,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_coEEsw.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_coEEsw.png",height=5,width=10)

dt_co_EUs <- data.table(cbind( qtlgridEst,MM_waveall_betaE_betaR_IR$betaptsE[,3], MM_waveall_betaE_betaR_IR$betaptsR[,3]))
names(dt_co_EUs) <- c("Quantile","Expansion","Recession")
dt_co_EUs_melted<-melt(dt_co_EUs,id="Quantile")
ggplot( dt_co_EUs_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Coefficient Value") + 
	xlab("Earnings Growth Quantile")  +
	scale_color_manual(values = c("blue","red")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.6,0.2),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_coEUsw.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_coEUsw.png",height=5,width=10)

dt_co_UEs <- data.table(cbind( qtlgridEst,MM_waveall_betaE_betaR_IR$betaptsE[,2], MM_waveall_betaE_betaR_IR$betaptsR[,2]))
names(dt_co_UEs) <- c("Quantile","Expansion","Recession")
dt_co_UEs_melted<-melt(dt_co_UEs,id="Quantile")
ggplot( dt_co_UEs_melted, aes(x=Quantile,y=value,colour=variable) ) + 
	geom_line(size=2) + 
	ylab("Coefficient Value") + 
	xlab("Earnings Growth Quantile")  +
	scale_color_manual(values = c("blue","red")) +
	theme(legend.title = element_blank(),
		  legend.position = c(0.2,0.6),
		  legend.background = element_rect(linetype = "solid",color = "black"))
ggsave("./Figures/MMwave_all_coUEsw.eps",height=5,width=10)
ggsave("./Figures/MMwave_all_coUEsw.png",height=5,width=10)


#Tables Out---------------------


#Moments:
MM_moments <- cbind(MM_waveall_betaE_betaR_IR$wc_exp_mmts,MM_waveall_betaE_betaR_IR$wc_rec_mmts,
					MM_waveall_betaE_betaR_IR$wc_IR_mmts,MM_waveall_betaE_betaR_IR$wc_IR_sw_mmts)
#MM_moments <- data.table(MM_moments)

MM_EUE_moments <- cbind(MM_waveallEUE_betaE_betaR_IR$wc_exp_mmts,MM_waveallEUE_betaE_betaR_IR$wc_rec_mmts,
						MM_waveallEUE_betaE_betaR_IR$wc_IR_mmts,MM_waveallEUE_betaE_betaR_IR$wc_IR_sw_mmts)
#MM_EUE_moments <- data.table(MM_EUE_moments)



colnames(MM_moments) <- c("Exp", "Rec" , "CF: $\\beta^E$, $\\indic^R$", 
						  "\\parbox{1\\columnwidth}{\\vspace{3px} CF: $\\beta^E$, $\\indic^R$  \\\\ + Switch Rtn\\vspace{3px} }") 
					#c("Exp","Rec","CF-Rec Flows",
					 #  "CF-Rec Flows + Sw Rets")
rownames(MM_moments) <- moment_names
MM_moments <- xtable(MM_moments, digits=2, 
					 align="l|ll|ll", caption="Moments of Machado-Mata Deocmpostion (EU,UE-view)\\label{tab:MM_moments}")
print(MM_moments,include.rownames=T, include.colnames = T ,  sanitize.colnames.function = function(x) {x} ,sanitize.text.function=function(x) {x},
	  hline.after= c(-1,-1,0,nrow(MM_moments)), file="./Figures/MM_moments.tex")
##
colnames(MM_EUE_moments) <- c("Exp", "Rec" , "CF: $\\beta^E$, $\\indic^R$", 
						   "\\parbox{1\\columnwidth}{\\vspace{3px} CF: $\\beta^E$, $\\indic^R$  \\\\ + Switch Rtn\\vspace{3px} }") 
						#c("Exp","Rec","CF-Rec Flows",
					    #"CF-Rec Flows + Sw Rets")
rownames(MM_EUE_moments) <- moment_names
MM_EUE_moments <- xtable(MM_EUE_moments, digits=2, 
					 align="l|ll|ll", caption="Moments of Machado-Mata Deocmpostion (EUE-view)\\label{tab:MM_EUE_moments}")
print(MM_EUE_moments,include.rownames=T, include.colnames = T , sanitize.colnames.function = function(x) {x} ,sanitize.text.function=function(x) {x},
	  hline.after= c(-1,-1,0,nrow(MM_EUE_moments)), file="./Figures/MM_EUE_moments.tex")

#quantile tables out:
MM_all_tab <- data.table(cbind( mmtabqtls,(dist_IR),dist_rec,dist_exp,dist_rec- dist_exp, 
								dist_pct,dist_pct_sw,dist_pct_un))
names(MM_all_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MM_all_tab) <- mmtabqtls
MM_all_tab <- xtable(MM_all_tab, digits=2, 
				 align="ll|lll|l|lll", caption="Machado-Mata, including unemployment \\label{tab:MM_all_tab}")
print(MM_all_tab,include.rownames=T, hline.after= c(0,nrow(MM_all_tab)), file="./Figures/MM_all.tex")

#table out:
MM_allEUE_tab <- data.table(cbind( mmtabqtls,(distEUE_IR),distEUE_rec,distEUE_exp,distEUE_rec- distEUE_exp, 
								   distEUE_pct,distEUE_pct_sw,distEUE_pct_un))
names(MM_allEUE_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
rownames(MM_allEUE_tab) <- mmtabqtls
MM_allEUE_tab <- xtable(MM_allEUE_tab, digits=2, 
					 align="ll|lll|l|lll", caption="Machado-Mata, connecting unemployment spells \\label{tab:MM_allEUE_tab}")
print(MM_allEUE_tab,include.rownames=T, hline.after= c(0,nrow(MM_allEUE_tab)), file="./Figures/MM_allEUE.tex")

