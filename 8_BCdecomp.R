# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(xtable)
library(Hmisc)
library(quantreg)
library(ggplot2)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Results"

setwd(wd0)


wt <- "truncweight"
wc <- "wagechangeEUE_wave"

qtlgridEst  <- c(seq(.02,.1,.02),seq(0.15,0.85,0.05),seq(.9,.98,.02))
qtlgridOut <- seq(.02,0.98,0.01)
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
	
	wcDF <- subset(wcDF,wt>0)
	
	# setup subgroup indices
	if(NS==7){
		# 7 subgroups, Sw X (EE UE EU) + stay, sets up conditional distributions.
		wcDF[wcDF$switch==T & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F , s := 1]
		wcDF[wcDF$switch==T & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F , s := 2]
		wcDF[wcDF$switch==T & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T , s := 3]
		wcDF[wcDF$switch==F & wcDF$EE==T & wcDF$UE==F & wcDF$EU ==F , s := 4]
		wcDF[wcDF$switch==F & wcDF$EE==F & wcDF$UE==T & wcDF$EU ==F , s := 5]
		wcDF[wcDF$switch==F & wcDF$EE==F & wcDF$UE==F & wcDF$EU ==T , s := 6]
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
		wcDF[  wcDF$EE==T & wcDF$EU ==F & wcDF$UE ==F , s := 1]
		wcDF[  wcDF$EE==F & wcDF$EU ==T & wcDF$UE ==F , s := 2]
		wcDF[  wcDF$EE==F & wcDF$EU ==F & wcDF$UE ==T , s := 3]
		wcDF[!(wcDF$EE==T | wcDF$EU ==T | wcDF$UE ==T), s := 4]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
	}else if(NS==5){
		# 4 subgroups, Sw X (EE EU) + stay, sets up conditional distributions.
		wcDF[wcDF$switch==T & wcDF$EE==T & wcDF$EU ==F , s := 1]
		wcDF[wcDF$switch==T & wcDF$EE==F & wcDF$EU ==T , s := 2]
		wcDF[wcDF$switch==F & wcDF$EE==T & wcDF$EU ==F , s := 3]
		wcDF[wcDF$switch==F & wcDF$EE==F & wcDF$EU ==T , s := 4]
		wcDF[               !(wcDF$EE==T | wcDF$EU ==T), s := 5]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(wcDF$s), s5 := ifelse(s==5,1,0)]		
	}else if(NS==6){
		# 6 subgroups, Sw X (EE EU Stay) , sets up conditional distributions.
		wcDF[wcDF$switch==T & wcDF$EE==T & wcDF$EU ==F , s := 1]
		wcDF[wcDF$switch==T & wcDF$EE==F & wcDF$EU ==T , s := 2]
		wcDF[wcDF$switch==F & wcDF$EE==T & wcDF$EU ==F , s := 3]
		wcDF[wcDF$switch==F & wcDF$EE==F & wcDF$EU ==T , s := 4]
		wcDF[wcDF$switch==F & !(wcDF$EE==T | wcDF$EU ==T), s := 5]
		wcDF[wcDF$switch==T & !(wcDF$EE==T | wcDF$EU ==T), s := 6]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
		wcDF[!is.na(wcDF$s), s5 := ifelse(s==5,1,0)]		
		wcDF[!is.na(wcDF$s), s6 := ifelse(s==6,1,0)]
	}else if(NS==4){
		# 4 subgroups, (EE EU UE) + stay, sets up conditional distributions.
		wcDF[  wcDF$EE==T & wcDF$EU ==F & wcDF$UE ==F , s := 1]
		wcDF[  wcDF$EE==F & wcDF$EU ==T & wcDF$UE ==F , s := 2]
		wcDF[  wcDF$EE==F & wcDF$EU ==F & wcDF$UE ==T , s := 3]
		wcDF[!(wcDF$EE==T | wcDF$EU ==T | wcDF$UE ==T), s := 4]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
		wcDF[!is.na(wcDF$s), s4 := ifelse(s==4,1,0)]
	}else if(NS==3){
		# 4 subgroups, (EE EU) + stay, sets up conditional distributions.
		wcDF[  wcDF$EE==T & wcDF$EU ==F & wcDF$UE ==F , s := 1]
		wcDF[  wcDF$EE==F & wcDF$EU ==T , s := 2]
		wcDF[!(wcDF$EE==T | wcDF$EU ==T), s := 3]
		wcDF[!is.na(wcDF$s), s1 := ifelse(s==1,1,0)]
		wcDF[!is.na(wcDF$s), s2 := ifelse(s==2,1,0)]
		wcDF[!is.na(wcDF$s), s3 := ifelse(s==3,1,0)]
	}
	if(durEU==T){
		wcDF[ !(wcDF$EU==T), dur := 0]
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
		sampE[,qi] <- sample(nrow(wcExp),nsampE,prob=wcExp$wt,replace=T) # weight the sampling. I will also weight the regression 
		sampR[,qi] <- sample(nrow(wcRec),nsampR,prob=wcRec$wt,replace=T) #,prob=wcRec$wt
	}
	
	# initialize output distributions:
	wc_IR_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_IR_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_BR_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_IR_sw_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_IR_un_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_rec_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_rec_moments <- array(0.,dim=c(Nmoments,Nsims))
	wc_exp_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_exp_moments <- array(0.,dim=c(Nmoments,Nsims))
	
	wc_BR_exstay_pctile <- array(0.,dim=c(length(qtlgridOut),Nsims))
	wc_BR_exstay_moments <- array(0.,dim=c(Nmoments,Nsims))
	
	for( simiter in seq(1,Nsims)){	
		
		set.seed(seedint+simiter*Nsims)		
		
		wcRec <- subset(wcDF, rec==T)
		wcExp <- subset(wcDF, rec==F)
		if(std_errs == T){
			datsampR <- sample(nrow(wcRec),nrow(wcRec),replace=T) # <-,prob=wcRec$wt do not weight the sample. I will weight the regression and stuff
			datsampE <- sample(nrow(wcExp),nrow(wcExp),replace=T) #   ,prob=wcExp$wt
			wcRec <- wcRec[datsampR,]
			wcExp <- wcExp[datsampE,]
		}
		
		if(durEU==T){
			regform <- formula(paste(c("wc~factor(s)","dur","0"),collapse=" + ") )
		}else{
			regform <- formula(paste(c("wc~factor(s)","0"),collapse=" + ") )
		}
		rhere <- rq( regform ,tau= qtlgridEst, data=wcRec, weights = wt, method="fn")
		betaptsR = t(rhere$coefficients)
	
		rhere <- rq( regform ,tau= qtlgridEst, data=wcExp, weights = wt, method="fn") #
		betaptsE = t(rhere$coefficients)
		rm(rhere)
	
		if(NS==7){
			colnames(betaptsE) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay")
			colnames(betaptsR) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay")
		}else if(NS==8){
			colnames(betaptsE) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay_nosw","stay_sw")
			colnames(betaptsR) <-c("EE_sw","UE_sw","EU_sw","EE_nosw","UE_nosw","EU_nosw","stay_nosw","stay_sw")
		}else if(NS==5){
			colnames(betaptsE) <-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay")
			colnames(betaptsR) <-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay")
		}else if(NS==6){
			colnames(betaptsE) <-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay_nosw","stay_sw")
			colnames(betaptsR) <-c("EE_sw","EU_sw","EE_nosw","EU_nosw","stay_nosw","stay_sw")
		}else if(NS==4){
			colnames(betaptsE) <-c("EE","UE","EU","stay")
			colnames(betaptsR) <-c("EE","UE","EU","stay")
		}else if(NS==3){
			colnames(betaptsE) <-c("EE","EU","stay")
			colnames(betaptsR) <-c("EE","EU","stay")
		}
		
		
		betaE <- array(0.,dim=c(NS,length(qtlgridSamp)) )
		betaR <- array(0.,dim=c(NS,length(qtlgridSamp)) )
		for(si in seq(1,NS)){
			gpE <- c(T,betaptsE[2:length(qtlgridEst),si]>=betaptsE[1:length(qtlgridEst)-1,si]) #ensure monotonicity
			gpR <- c(T,betaptsR[2:length(qtlgridEst),si]>=betaptsR[1:length(qtlgridEst)-1,si]) #ensure monotonicity
			if(min(qtlgridSamp)<min(qtlgridEst[gpR]) | max(qtlgridSamp)>max(qtlgridEst[gpR])){
				betaE[si,] <- spline(x=c(min(qtlgridSamp) ,qtlgridEst[gpE], max(qtlgridSamp)), 
									 y=c(min(betaptsE[gpE,si]), betaptsE[gpE,si],max(betaptsE[gpE,si])), method="hyman", xout=qtlgridSamp)$y
			}
			if(min(qtlgridSamp)<min(qtlgridEst[gpR]) | max(qtlgridSamp)>max(qtlgridEst[gpR])){
				betaR[si,] <- spline(x=c(min(qtlgridSamp) ,qtlgridEst[gpR], max(qtlgridSamp)), 
									 y=c(min(betaptsR[gpR,si]),betaptsR[gpR,si],max(betaptsR[gpR,si])), method="hyman", xout=qtlgridSamp)$y	
			}
			if( min(qtlgridSamp)>=min(qtlgridEst) & max(qtlgridSamp)<=max(qtlgridEst[gpR]) & max(qtlgridSamp)<=max(qtlgridEst[gpE])){
				betaE[si,] <- spline(x=qtlgridEst[gpE], y=betaptsE[gpE,si], method="hyman", xout=qtlgridSamp)$y
				betaR[si,] <- spline(x=qtlgridEst[gpR], y=betaptsR[gpR,si], method="hyman", xout=qtlgridSamp)$y
			}
			if( any(gpE ==F) | any(gpR==F) ){
				warning(c('non-monotone regression coefficients, ',as.character(sum(gpE==F) + sum(gpR==F)), ' time(s)') )
			}
		}

		wc_IR <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution
		wc_BR <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution
		
		qi = 1
		for(q in qtlgridSamp){
			sumChar <- ""
			for(si in seq(1,NS) ){
				sumBetaE <- paste(sumChar,paste0("betaE[",as.character(si),",qi]*s",as.character(si) ) , collapse=" + ")
				sumBetaR <- paste(sumChar,paste0("betaR[",as.character(si),",qi]*s",as.character(si) ) , collapse=" + ")
			}
			if( durEU == T){
				sumBetaE <- paste(sumChar,paste0("betaE[",as.character(si),",qi]*s",as.character(si) ) , collapse=" + ")
				sumBetaR <- paste(sumChar,paste0("betaR[",as.character(si),",qi]*s",as.character(si) ) , collapse=" + ")
			}
			
			wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[ sampR[,qi] , eval(parse(sumBetaE))]
			wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[ sampE[,qi] , eval(parse(sumBetaR))]

			# if(NS == 6){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
			# 		  	betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
			# 		  	betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
			# 			betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6] 
			# }else if(NS==7){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
			# 			betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
			# 			betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 + 
			# 			betaR[7,qi]*s7]
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
			# 			betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 + 
			# 			betaE[7,qi]*s7]
			# }else if(NS==8){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
			# 			betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
			# 			betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 + 
			# 			betaR[7,qi]*s7 + betaR[8,qi]*s8]
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
			# 			betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 + 
			# 			betaE[7,qi]*s7 + betaE[8,qi]*s8]
			# }else if(NS==5){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
			# 			betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
			# 			betaE[3,qi]*s3 + betaE[4,qi]*s4 +
			# 			betaR[5,qi]*s5	]
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi], 
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
			# 			betaR[3,qi]*s3 + betaR[4,qi]*s4 +
			# 			betaE[5,qi]*s5	]
			# }else if(NS==4){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
			# 			betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
			# 			betaE[3,qi]*s3 + 
			# 			betaR[4,qi]*s4]
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi], 
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
			# 			betaR[3,qi]*s3 + 
			# 			betaE[4,qi]*s4]
			# }else if(NS==3){
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
			# 			betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
			# 			betaE[3,qi]*s3 ]
			# 	wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi], 
			# 			betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
			# 			betaR[3,qi]*s3 ]
			# }else{
			# 	wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ]  = 0.
			# 	for(si in seq(1,NS)){
			# 		wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],betaE[si,qi]*eval(as.name(paste0("s",si))) ] + wc_IR[ ((qi-1)*nsampR+1):(qi*nsampR) ] 
			# 		wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],betaR[si,qi]*eval(as.name(paste0("s",si))) ] + wc_BR[ ((qi-1)*nsampE+1):(qi*nsampE) ] 
			# 	}
			# }
			qi = qi+1
		}
		#clean up the space:
		wc_IR_pctile[,simiter] <- quantile(wc_IR,probs=qtlgridOut,na.rm=T)
		wc_BR_pctile[,simiter] <- quantile(wc_BR,probs=qtlgridOut,na.rm=T)
		wc_IR_moments[1,simiter] <- mean(wc_IR,na.rm = T)
		wc_IR_moments[2:5,simiter] <- wtd.4qtlmoments(xt=wc_IR,wt=array(1.,dim = dim(wc_IR)))
		wc_BR_moments[1,simiter] <- mean(wc_BR)
		wc_BR_moments[2:5,simiter] <- wtd.4qtlmoments(xt=wc_BR,wt=array(1.,dim = dim(wc_BR)))
		
		rm(wc_IR)
		rm(wc_BR)

		qi=1
		wc_rec <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
		for(q in qtlgridSamp){
			if(NS == 6){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
						betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6] 
			}else if(NS == 7){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
						betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 +
						betaR[7,qi]*s7]
			}else if(NS == 8){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
						betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 +
						betaR[7,qi]*s7 + betaR[8,qi]*s8]
			}else if(NS==4){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
						betaR[3,qi]*s3 + betaR[4,qi]*s4]
			}else if(NS==5){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
						betaR[3,qi]*s3 + betaR[4,qi]*s4 +
						betaR[5,qi]*s5]
			}else if(NS==3){
				wc_rec[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
						betaR[1,qi]*s1 + betaR[2,qi]*s2 + 
						betaR[3,qi]*s3 ]
		}
			qi = qi+1
		}
		#clean up the space:
		wc_rec_pctile[,simiter] <- quantile(wc_rec,probs=qtlgridOut,na.rm=T)
		wc_rec_moments[1,simiter] <- mean( wc_rec,na.rm = T )
		wc_rec_moments[2:5,simiter] <- wtd.4qtlmoments( xt=wc_rec,wt=array(1,dim=dim(wc_rec)) )
		rm(wc_rec)
		
		qi=1
		wc_exp <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
		for(q in qtlgridSamp){
			if(NS == 6){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
						betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
			}else if(NS == 7){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
						betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 +
						betaE[7,qi]*s7]
			}else if(NS == 8){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
						betaE[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 +
						betaE[7,qi]*s7 + betaE[8,qi]*s8]
			}else if(NS==4){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
						betaE[3,qi]*s3 + betaE[4,qi]*s4]
			}else if(NS==5){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi], 
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
						betaE[3,qi]*s3 + betaE[4,qi]*s4 +
						betaE[5,qi]*s5]
			}else if(NS==3){
				wc_exp[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi], 
						betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
						betaE[3,qi]*s3]
			}
			qi = qi+1
		}
		#clean up the space:
		wc_exp_pctile[,simiter] <- quantile(wc_exp,probs=qtlgridOut,na.rm=T)
		wc_exp_moments[1,simiter] <- mean( wc_exp,na.rm = T )
		wc_exp_moments[2:5,simiter] <- wtd.4qtlmoments( xt=wc_exp,wt=array(1,dim=dim(wc_exp)) )
		rm(wc_exp)
				
		#some additional counter-factuals
		if(NS==6|NS==8){
			wc_BR_exstay <- matrix(NA, nrow=nsampE*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
			qi=1
			for(q in qtlgridSamp){
				wc_BR_exstay[ ((qi-1)*nsampE+1):(qi*nsampE) ] <- wcExp[sampE[,qi],
								betaR[1,qi]*s1 + betaR[2,qi]*s2 + betaR[3,qi]*s3 + 
								betaR[4,qi]*s4 + betaE[5,qi]*s5 + betaR[6,qi]*s6] 
				qi = qi+1
			}
			wc_BR_exstay_pctile[,simiter] <- quantile(wc_BR_exstay,probs=qtlgridOut,na.rm=T)
			wc_BR_exstay_moments[1,simiter] <- mean( wc_BR_exstay,na.rm = T )
			wc_BR_exstay_moments[2:5,simiter] <- wtd.4qtlmoments( xt=wc_BR_exstay,wt=array(1,dim=dim(wc_BR_exstay)) )
			
		}
		
		if(no_occ ==F){
			qi=1
			wc_IR_sw <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only switch
			wc_IR_sw_moments <- matrix(NA,Nmoments,Nsims)
			for(q in qtlgridSamp){
				if(NS == 6){
					wc_IR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
								betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
								betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6] 
				}else if(NS == 7){
					wc_IR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
								betaE[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
								betaR[4,qi]*s4 + betaR[5,qi]*s5 + betaR[6,qi]*s6 + 
								betaR[7,qi]*s7]
				}else if(NS==4){
					wc_IR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
								betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
								betaR[3,qi]*s3 + betaR[4,qi]*s4]
				}else if(NS==5){
					wc_IR_sw[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
								betaE[1,qi]*s1 + betaE[2,qi]*s2 + 
								betaR[3,qi]*s3 + betaR[4,qi]*s4 + 
								betaR[5,qi]*s5]
				}
				qi = qi+1
			}
			#clean up the space:
			wc_IR_sw_pctile[,simiter] <- quantile(wc_IR_sw,probs=qtlgridOut,na.rm=T)
			wc_IR_sw_moments[1,simiter] <- mean( wc_IR_sw,na.rm = T )
			wc_IR_sw_moments[2:5,simiter] <- wtd.4qtlmoments( xt=wc_IR_sw,wt=array(1,dim=dim(wc_IR_sw)) )
			rm(wc_IR_sw)
			
			
			qi=1
			wc_IR_un <- matrix(NA, nrow=nsampR*length(qtlgridSamp),ncol=1) #storing the counter-factual distribution - only unemployment
			for(q in qtlgridSamp){
				if(NS == 6){
					wc_IR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
								betaR[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
								betaR[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6] 
				}else if(NS == 7){
					wc_IR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi],
								betaR[1,qi]*s1 + betaE[2,qi]*s2 + betaE[3,qi]*s3 + 
								betaR[4,qi]*s4 + betaE[5,qi]*s5 + betaE[6,qi]*s6 +
								betaR[7,qi]*s7]
				}else if(NS==4){
					wc_IR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
								betaR[1,qi]*s1 + betaE[2,qi]*s2 + 
								betaR[3,qi]*s3 + betaE[4,qi]*s4]
				}else if(NS==5){
					wc_IR_un[ ((qi-1)*nsampR+1):(qi*nsampR) ] <- wcRec[sampR[,qi], 
								betaR[1,qi]*s1 + betaE[2,qi]*s2 + 
								betaR[3,qi]*s3 + betaE[4,qi]*s4 +
								betaR[5,qi]*s5]
				}
				qi = qi+1
			}
			#clean up the space:
			wc_IR_un_pctile[,simiter] <- quantile(wc_IR_un,probs=qtlgridOut,na.rm=T)
			rm(wc_IR_un)
		}
	}#simiter	
	
	if(no_occ==F){	
		return(list(betaptsE = betaptsE,betaptsR=betaptsR,wc_IR = wc_IR_pctile, wc_IR_sw= wc_IR_sw_pctile, wc_IR_un= wc_IR_un_pctile,
				wc_rec = wc_rec_pctile,wc_exp = wc_exp_pctile, wc_BR = wc_BR_pctile,wc_BR_mmts = wc_BR_moments,wc_IR_mmts = wc_IR_moments, 
				wc_exp_mmts = wc_exp_moments, wc_rec_mmts = wc_rec_moments, wc_IR_sw_mmts = wc_IR_sw_moments))
	}else{
		return(list(betaptsE = betaptsE,betaptsR=betaptsR,wc_IR = wc_IR_pctile, wc_rec = wc_rec_pctile,wc_exp = wc_exp_pctile, wc_BR = wc_BR_pctile,
					wc_BR_mmts = wc_BR_moments,wc_IR_mmts = wc_IR_moments,wc_exp_mmts = wc_exp_moments, wc_rec_mmts = wc_rec_moments))
	}
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
	return( list(mean_hr,median_hr,mad_hr,GroenMeed_hr,Moors_hr))
}

########################################################################
########################################################################
## Seams version -------------------------------

DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))

DTseam <- subset(DTseam, changermo==T|stayer==T)

toKeep <- c("truncweight","cycweight","wpfinwgt","EU_wave","UE_wave","EE_wave","switchedOcc_wave","switched_wave","wagechange_wave","wagechangeEUE_wave",
			"max.unempdur_wave","wagechange_wvan","wagechange_anan","recIndic_wave","recIndic2_wave","changer","changermo","stayer","date")


# select toKeep columns only
#DTseam <- DTseam[, toKeep, with = FALSE]
DTseam <- subset(DTseam, is.finite(wagechange_wave) & is.finite(EU_wave) & is.finite(UE_wave)& is.finite(EE_wave))
DTseam[ , EE := EE_wave]
DTseam[ , EU := EU_wave]
DTseam[ , UE := UE_wave]
DTseam[ stayer ==T, switch := switchedOcc_wave]
DTseam[ changer==T, switch := switched_wave]
DTseam[ EU==T, dur := max.unempdur_wave]
DTseam[ !EU==T | !is.finite(dur), dur := 0.]

MM_waveall_betaE_betaR_IR <- MMdecomp(DTseam,7,"recIndic2_wave","wagechange_wave","truncweight",std_errs = MMstd_errs, no_occ = F)
MM_waveallEUE_betaE_betaR_IR <- MMdecomp(DTseam,6,"recIndic2_wave",wcname=wc,wtname=wt,std_errs = MMstd_errs, no_occ = T)

saveRDS(MM_waveall_betaE_betaR_IR,paste0(outputdir,"/MM_waveall.RData"))
saveRDS(MM_waveallEUE_betaE_betaR_IR,paste0(outputdir,"/MM_waveallEUE.RData"))

MM_wavenooc_betaE_betaR_IR <- MMdecomp(DTseam,4,"recIndic_wave","wagechange_wave","truncweight",std_errs = MMstd_errs, no_occ = T)
MM_wavenoocEUE_betaE_betaR_IR <- MMdecomp(DTseam,3,"recIndic_wave","wagechangeEUE_wave","truncweight",std_errs = MMstd_errs, no_occ = T)

saveRDS(MM_wavenooc_betaE_betaR_IR,paste0(outputdir,"/MM_wavenooc.RData"))
saveRDS(MM_wavenoocEUE_betaE_betaR_IR,paste0(outputdir,"/MM_wavenoocEUE.RData"))


#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EU,UE all observations (7 categories) -----------------------------------
wcExp <- subset(DTseam,recIndic_wave==F)
wcRec <- subset(DTseam,recIndic_wave==T)
dist_exp      <- wcExp[ , wtd.quantile(wagechange_wave,probs=qtlgridOut,weights=truncweight, na.rm=T)]
dist_rec      <- wcRec[ , wtd.quantile(wagechange_wave,probs=qtlgridOut,weights=truncweight, na.rm=T)]
dist_IR       <- rowMeans(MM_waveall_betaE_betaR_IR$wc_IR   )
if(MMstd_errs ==T){
	CI_IR         <- apply(MM_waveall_betaE_betaR_IR$wc_IR,1,quantile,0.95) - apply(MM_waveall_betaE_betaR_IR$wc_IR,1,quantile,0.05)
}
dist_BR       <- rowMeans(MM_waveall_betaE_betaR_IR$wc_BR   )
dist_IR_sw    <- rowMeans(MM_waveall_betaE_betaR_IR$wc_IR_un)
dist_IR_un    <- rowMeans(MM_waveall_betaE_betaR_IR$wc_IR_sw)
dist_sim_rec  <- rowMeans(MM_waveall_betaE_betaR_IR$wc_rec  )
dist_sim_exp  <- rowMeans(MM_waveall_betaE_betaR_IR$wc_exp  )

#moving window to summarize the quantiles
Nwin = 87
win_IR      <- array(0.,dim=Nwin)
win_BR      <- array(0.,dim=Nwin)
win_IR_un   <- array(0.,dim=Nwin)
win_IR_sw   <- array(0.,dim=Nwin)
win_exp     <- array(0.,dim=Nwin)
win_rec     <- array(0.,dim=Nwin)
win_rec_dat <- array(0.,dim=Nwin)
win_exp_dat <- array(0.,dim=Nwin)
di = 1
for(lb in seq(1,Nwin,1)){
	ub = lb+10
	win_IR[di]      <- mean(dist_IR     [ lb:ub ])
	win_BR[di]      <- mean(dist_BR     [ lb:ub ])
	win_rec[di]     <- mean(dist_sim_rec[ lb:ub ])
	win_exp[di]     <- mean(dist_sim_exp[ lb:ub ])
	win_IR_un[di]   <- mean(dist_IR_un  [ lb:ub ])
	win_IR_sw[di]   <- mean(dist_IR_sw  [ lb:ub ])	
	win_rec_dat[di] <- mean(dist_rec    [ lb:ub ])
	win_exp_dat[di] <- mean(dist_exp    [ lb:ub ])	
	di= di+1
}
win_dif_rec_exp_dat <- win_rec_dat - win_exp_dat
win_dif_rec_exp_sim <- win_rec     - win_exp
win_dif_IR_exp      <- win_IR      - win_exp
win_dif_IR_un_exp   <- win_IR_un   - win_exp
win_dif_IR_sw_exp   <- win_IR_sw   - win_exp
win_dif_BR_exp      <- win_BR      - win_exp

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EUE wage changes (5 categories) -----------------------
distEUE_exp      <- wcExp[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distEUE_rec      <- wcRec[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distEUE_IR       <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_IR   )
distEUE_BR       <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_BR   )
distEUE_IR_sw    <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_IR_un)
distEUE_IR_un    <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_IR_sw)
distEUE_sim_rec  <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_rec  )
distEUE_sim_exp  <- rowMeans(MM_waveallEUE_betaE_betaR_IR$wc_exp  )

#moving window to summarize the quantiles
winEUE_IR      <- array(0.,dim=Nwin)
winEUE_BR      <- array(0.,dim=Nwin)
winEUE_IR_un   <- array(0.,dim=Nwin)
winEUE_IR_sw   <- array(0.,dim=Nwin)
winEUE_exp     <- array(0.,dim=Nwin)
winEUE_rec     <- array(0.,dim=Nwin)
winEUE_rec_dat <- array(0.,dim=Nwin)
winEUE_exp_dat <- array(0.,dim=Nwin)
di = 1
for(lb in seq(1,Nwin,1)){
	ub = lb+10
	winEUE_IR[di]      <- mean(distEUE_IR     [ lb:ub ])
	winEUE_BR[di]      <- mean(distEUE_BR     [ lb:ub ])
	winEUE_rec[di]     <- mean(distEUE_sim_rec[ lb:ub ])
	winEUE_exp[di]     <- mean(distEUE_sim_exp[ lb:ub ])
	winEUE_IR_un[di]   <- mean(distEUE_IR_un  [ lb:ub ])
	winEUE_IR_sw[di]   <- mean(distEUE_IR_sw  [ lb:ub ])	
	winEUE_rec_dat[di] <- mean(distEUE_rec    [ lb:ub ])
	winEUE_exp_dat[di] <- mean(distEUE_exp    [ lb:ub ])	
	di= di+1
}
winEUE_dif_rec_exp_dat <- winEUE_rec_dat - winEUE_exp_dat
winEUE_dif_rec_exp_sim <- winEUE_rec     - winEUE_exp
winEUE_dif_IR_exp      <- winEUE_IR      - winEUE_exp
winEUE_dif_BR_exp      <- winEUE_BR      - winEUE_exp

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EU,UE wage changes - no occupations (4 categories) --------------------
distnooc_exp      <- wcExp[ , wtd.quantile(wagechange_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distnooc_rec      <- wcRec[ , wtd.quantile(wagechange_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distnooc_IR       <- rowMeans(MM_wavenooc_betaE_betaR_IR$wc_IR   )
distnooc_BR       <- rowMeans(MM_wavenooc_betaE_betaR_IR$wc_BR   )
distnooc_sim_rec  <- rowMeans(MM_wavenooc_betaE_betaR_IR$wc_rec  )
distnooc_sim_exp  <- rowMeans(MM_wavenooc_betaE_betaR_IR$wc_exp  )

winnooc_IR      <- array(0.,dim=Nwin)
winnooc_BR      <- array(0.,dim=Nwin)
winnooc_exp     <- array(0.,dim=Nwin)
winnooc_rec     <- array(0.,dim=Nwin)
winnooc_rec_dat <- array(0.,dim=Nwin)
winnooc_exp_dat <- array(0.,dim=Nwin)
di = 1
for(lb in seq(1,Nwin,1)){
	ub = lb+10
	winnooc_IR[di]      <- mean(distnooc_IR     [ lb:ub ])
	winnooc_BR[di]      <- mean(distnooc_BR     [ lb:ub ])
	winnooc_rec[di]     <- mean(distnooc_sim_rec[ lb:ub ])
	winnooc_exp[di]     <- mean(distnooc_sim_exp[ lb:ub ])
	winnooc_rec_dat[di] <- mean(distnooc_rec    [ lb:ub ])
	winnooc_exp_dat[di] <- mean(distnooc_exp    [ lb:ub ])	
	di= di+1
}
winnooc_dif_rec_exp_dat <- winnooc_rec_dat - winnooc_exp_dat
winnooc_dif_rec_exp_sim <- winnooc_rec     - winnooc_exp
winnooc_dif_IR_exp      <- winnooc_IR      - winnooc_exp
winnooc_dif_BR_exp      <- winnooc_BR      - winnooc_exp

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EUE wage changes - no occupations (3 categories) --------------
distnoocEUE_exp      <- wcExp[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distnoocEUE_rec      <- wcRec[ , wtd.quantile(wagechangeEUE_wave,probs=seq(0.02,.98,0.01),weights=truncweight, na.rm=T)]
distnoocEUE_IR       <- rowMeans(MM_wavenoocEUE_betaE_betaR_IR$wc_IR   )
distnoocEUE_BR       <- rowMeans(MM_wavenoocEUE_betaE_betaR_IR$wc_BR   )
distnoocEUE_sim_rec  <- rowMeans(MM_wavenoocEUE_betaE_betaR_IR$wc_rec  )
distnoocEUE_sim_exp  <- rowMeans(MM_wavenoocEUE_betaE_betaR_IR$wc_exp  )

winnoocEUE_IR      <- array(0.,dim=Nwin)
winnoocEUE_BR      <- array(0.,dim=Nwin)
winnoocEUE_exp     <- array(0.,dim=Nwin)
winnoocEUE_rec     <- array(0.,dim=Nwin)
winnoocEUE_rec_dat <- array(0.,dim=Nwin)
winnoocEUE_exp_dat <- array(0.,dim=Nwin)
di = 1
for(lb in seq(1,Nwin,1)){
	ub = lb+10
	winnoocEUE_IR[di]      <- mean(distnoocEUE_IR     [ lb:ub ])
	winnoocEUE_BR[di]      <- mean(distnoocEUE_BR     [ lb:ub ])
	winnoocEUE_rec[di]     <- mean(distnoocEUE_sim_rec[ lb:ub ])
	winnoocEUE_exp[di]     <- mean(distnoocEUE_sim_exp[ lb:ub ])
	winnoocEUE_rec_dat[di] <- mean(distnoocEUE_rec    [ lb:ub ])
	winnoocEUE_exp_dat[di] <- mean(distnoocEUE_exp    [ lb:ub ])	
	di= di+1
}
winnoocEUE_dif_rec_exp_dat <- winnoocEUE_rec_dat - winnoocEUE_exp_dat
winnoocEUE_dif_rec_exp_sim <- winnoocEUE_rec     - winnoocEUE_exp
winnoocEUE_dif_IR_exp      <- winnoocEUE_IR      - winnoocEUE_exp
winnoocEUE_dif_BR_exp      <- winnoocEUE_BR      - winnoocEUE_exp



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



#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# only month changes data set
#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
# wagechanges <- readRDS(paste0(datadir,"/balancedwagechanges.RData"))
# wagechanges <- subset(wagechanges, select=keep)
# 
# 
# # create vector of recession dates : above 7%
# recDates3 <- as.Date(c("1991-10-01", "1993-07-01","2008-12-01", "2013-11-01"))
# wagechanges[, recIndic3 := (date > recDates3[1] & date < recDates3[2]) | 
# 				(date > recDates3[3] & date < recDates3[4])]
# 
# 
# shift_share <- DHLdecomp(wagechanges,6,"recIndic","wagechange","balanceweight")
# pct_share <- shift_share[[3]]/shift_share[[1]]
# pct_shift <- shift_share[[2]]/matrix(shift_share[[1]],nrow=9,ncol=6 )
# 
# shift_share_EUE <- DHLdecomp(wagechanges,4,"recIndic","wagechange_EUE","balanceweightEUE")
# pct_share_EUE <- shift_share_EUE[[3]]/shift_share_EUE[[1]]
# pct_shift_EUE <- shift_share_EUE[[2]]/matrix(shift_share_EUE[[1]],nrow=9,ncol=4 )
# 
# MM_betaE_betaR_IR <- MMdecomp(wagechanges,6,"recIndic","wagechange","balanceweight")
# 
# MMEUE_betaE_betaR_IR <- MMdecomp(subset(wagechanges,EE==T|EU==T),4,"recIndic","wagechange_EUE","balanceweightEUE")
# 
# mmtabqtls <- seq(0.1,0.9,0.1)
# wcExp <- subset(wagechanges,recIndic==F)
# wcRec <- subset(wagechanges,recIndic==T)
# distEUE_exp <- wcExp[EE==T | EU==T, wtd.quantile(wagechange_EUE,probs=mmtabqtls,weights=balanceweightEUE, na.rm=T)]
# distEUE_rec <- wcRec[EE==T | EU==T, wtd.quantile(wagechange_EUE,probs=mmtabqtls,weights=balanceweightEUE, na.rm=T)]
# distEUE_IR  <- MMEUE_betaE_betaR_IR$wc_IR[mmtabqtls*100] 
# distEUE_IR_sw  <- MMEUE_betaE_betaR_IR$wc_IR_un[mmtabqtls*100]
# distEUE_IR_un  <- MMEUE_betaE_betaR_IR$wc_IR_sw[mmtabqtls*100]
# distEUE_pct <- (distEUE_IR - distEUE_exp)/(distEUE_rec - distEUE_exp)
# distEUE_pct_sw <- (distEUE_IR_sw - distEUE_exp)/(distEUE_rec - distEUE_exp)
# distEUE_pct_un <- (distEUE_IR_un - distEUE_exp)/(distEUE_rec - distEUE_exp)
# 
# dist_exp <- wtd.quantile(wcExp$wagechange,probs=mmtabqtls,weights=wcExp$balanceweight, na.rm=T)
# dist_rec <- wtd.quantile(wcRec$wagechange,probs=mmtabqtls,weights=wcRec$balanceweight, na.rm=T)
# dist_IR  <- MM_betaE_betaR_IR$wc_IR[mmtabqtls*100]
# dist_IR_sw  <- MM_betaE_betaR_IR$wc_IR_sw[mmtabqtls*100]
# dist_IR_un  <- MM_betaE_betaR_IR$wc_IR_un[mmtabqtls*100] 
# dist_pct <- (dist_IR - dist_exp)/(dist_rec - dist_exp)
# dist_pct_sw <- (dist_IR_sw - dist_exp)/(dist_rec - dist_exp)
# dist_pct_un <- (dist_IR_un - dist_exp)/(dist_rec - dist_exp)
# 
# 
# #-- spit these out into tables
# 
# #DHW: EE, EU, UE
# pct_share_shift <- data.table(cbind(seq(0.1,0.9,0.1),pct_share,pct_shift))
# names(pct_share_shift) <- c("Decile","Share","EE","UE","EU","EE","UE","EU")
# addswitched <-list(pos = list(-1),command="\\hline\\hline& & \\multicolumn{3}{c}{Switched Occ} & \\multicolumn{3}{c}{Not Switched Occ} \\\\ ")
# pct_share_shift <- xtable(pct_share_shift, label="tab:pct_share_shift", digits=2, 
# 						  align="ll|l|lll|lll", caption="Shift-Share (DHL) decomposition, including unemployment")
# print(pct_share_shift,add.to.row=addswitched, include.rownames=F, hline.after= c(0,nrow(pct_shift)), file="./Figures/pct_share_shift.tex")
# 
# #DHW: EE, EUE
# pct_share_shift_EUE <- data.table(cbind(seq(0.1,0.9,0.1),pct_share_EUE,pct_shift_EUE))
# names(pct_share_shift_EUE) <- c("Decile","Share","EE","EUE","EE","EUE")
# addswitched <-list(pos = list(-1),command="\\hline\\hline& & \\multicolumn{2}{c}{Switched Occ} & \\multicolumn{2}{c}{Not Switched Occ} \\\\ ")
# pct_share_shift_EUE <- xtable(pct_share_shift_EUE, label="tab:pct_share_shift_EUE", digits=2, 
# 							  align="ll|l|ll|ll", caption="Shift-Share (DHL) decomposition, connecting across unemployment")
# print(pct_share_shift_EUE,add.to.row=addswitched, include.rownames=F, hline.after= c(0,nrow(pct_shift)), file="./Figures/pct_share_shift_EUE.tex")
# 
# #MM : EE, EU, UE
# MM_tab <- data.table(cbind( mmtabqtls,(dist_IR),dist_rec,dist_exp,dist_rec- dist_exp, 
# 							dist_pct,dist_pct_sw,dist_pct_un))
# names(MM_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
# rownames(MM_tab) <- mmtabqtls
# MM_tab <- xtable(MM_tab, digits=2, 
# 				 align="ll|lll|l|lll", caption="Machado-Mata, including unemployment \\label{tab:MM_tab}")
# print(MM_tab,include.rownames=T, hline.after= c(0,nrow(MM_tab)), file="./Figures/MM.tex")
# ks.test(wcRec$wagechange,wcExp$wagechange,alternative = "greater")
# # plot the coefficients
# MMcoef <- data.table(cbind( mmtabqtls,MM_betaE_betaR_IR$betaptsE,MM_betaE_betaR_IR$betaptsE))
# names(MMcoef) <- c("Quantiles","ExpSwEE","ExpSwUE","ExpSwEU","ExpNSwEE","ExpNSwUE","ExpNSwEU","RecSwEE","RecSwUE","RecSwEU","RecNSwEE","RecNSwUE","RecNSwEU")
# EUfrac <- wagechanges[recIndic==F, wtd.mean(EU,na.rm=T,weights=balanceweight)]
# UEfrac <- wagechanges[recIndic==F, wtd.mean(UE,na.rm=T,weights=balanceweight)]
# EEfrac <- wagechanges[recIndic==F, wtd.mean(EE,na.rm=T,weights=balanceweight)]
# Swfrac <- wagechanges[recIndic==F, wtd.mean(switchedOcc,na.rm=T,weights=balanceweight)]
# 
# 
# #MM : EE,EUE
# MMEUE_tab <- data.table(cbind( mmtabqtls,(distEUE_IR),distEUE_rec,distEUE_exp,distEUE_rec- distEUE_exp, 
# 							   distEUE_pct,distEUE_pct_sw,distEUE_pct_un ))
# names(MMEUE_tab) <- c("Quantile","CF\ Rec","Rec","Exp","Rec-Exp","Pct\ CF","Pct\ CF\ OccSw","Pct\ CF\ Unemp")
# rownames(MMEUE_tab) <- mmtabqtls
# MMEUE_tab <- xtable(MMEUE_tab, digits=2, 
# 					align="ll|lll|l|lll", caption="Machado-Mata, connecting across unemployment \\label{tab:MMEUE_tab}")
# print(MMEUE_tab,include.rownames=F, hline.after= c(0,nrow(MMEUE_tab)), file="./Figures/MMEUE.tex")
# 
# ks.test(wcRec$wagechange_EUE,wcExp$wagechange_EUE,alternative = "greater")
# 
# rm(list = c("wcRec","wcExp"))
# 
# 
# rm(list=ls())