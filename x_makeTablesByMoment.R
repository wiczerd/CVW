# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)
library(xtable)
library(quantreg)
library(texreg)
library(ggplot2)

wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Figures"
setwd(wd0)


recDef <- "recIndic_wave"
wt <- "truncweight"
wclab <- "res" #raw or res

demolbl <- 0 #or choose number from categories in demotxt
demotxt <- c("Young", "Prime","Old","HS","Col","Male","Female")

bootse <- T #compute bootstrapped standard errors or no?
seedint = 941987
Nsim = 40 #make this bigger: just for diagnostic!!!!!!!!!!!!!!!!!!!

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
	Moors  <- wtd.Moors(xt,wt,samp = T)
	return(c(median,mad,GrnMd,Moors))
}


CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
CPSunempRt$unrt <- CPSunempRt$unrt/100

#*********************************************************************

toKeep_wave <- c("switchedOcc_wave",
            "ageGrp","HSCol",
            "recIndic","recIndic_wave","recIndic2_wave","recIndic_stint",
            "wagechange_month","wagechange_wave","wagechangeEUE_wave","rawwgchange_wave","rawwgchangeEUE_wave",
            "wagechange_wave_bad2","wagechange_wave_low","wagechange_wave_high","wagechange_wave_jcbad",
            "EE_wave","EU_wave","UE_wave","changer","stayer",
            "wpfinwgt","perwt","cycweight","truncweight","cleaningtruncweight",
			"lfstat_wave","next.lfstat_wave","wave","id","date","panel")
DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))
DTseam <- merge(DTseam, CPSunempRt, by = "date", all.x = TRUE)

# select toKeep columns only
DTseam <- DTseam[, toKeep_wave, with = FALSE]
DTseam <- subset(DTseam, (stayer|changer))

DTseam[ , idx := as.integer(as.factor(id))]

# select deomgraphic
if(demolbl==1){
	DTseam[(stayer==T|changer==T), demo:= (ageGrp==1)]
}else if(demolbl==2){
	DTseam[(stayer==T|changer==T), demo:= (ageGrp==2)]
}else if(demolbl==3){
	DTseam[(stayer==T|changer==T), demo:= (ageGrp==3)]
}else if(demolbl==4){
	DTseam[(stayer==T|changer==T), demo:= (HSCol==1)]
}else if(demolbl==5){
	DTseam[(stayer==T|changer==T), demo:= (HSCol==2)]
}else if(demolbl==6){
	DTseam[(stayer==T|changer==T), demo:= (female==F)]
}else if(demolbl==7){
	DTseam[(stayer==T|changer==T), demo:= (female==T)]
}else{
	DTseam[(stayer==T|changer==T), demo := T]
}

if(recDef == "recIndic_wave"){
	reclab = "NBER"
}else if( recDef == "recIndic2_wave"){
	reclab = "urt"
}else if(wc == "recIndic_stint"){
	reclab = "recstint"
}


#Wage moments --------------------------------------------------------
Nmoments <- 5 # the number of moemnts to spit out for each subset
Nsubsamp <- 6 # the number of population subsamples: All, stayer, EE, EUE, All-EUUE, EUUE
Nt <- 3 # time periods, all, exp, rec

# how to weights EUE's: 2x for an EUUE
DTseam[ ,wtEUE:= eval(as.name(wt))]
DTseam[ UE_wave==T,wtEUE:= 0.]
DTseam[ EU_wave==T,wtEUE:= 2.*eval(as.name(wt))]
origwt = wt

if(bootse == T){
	set.seed(seedint)
	#draw the sample
	nsamp  = DTseam[ , max(idx)]
	se_wavemoments <- array(0.,dim = c(Nsubsamp,Nt,Nmoments,Nsim))
	sampidx <-array(1L, dim = c(nsamp,Nsim))
	for(si in seq(1,Nsim)){ sampidx[ , si] = sample(nsamp,nsamp,replace = T) }
}
tab_wavemoments <- array(0., dim=c(Nsubsamp,Nt,Nmoments))

for( si in seq(0,bootse*Nsim) ){
	if(si>0){
		set.seed(seedint+si)
		DThr <- DTseam[idx %in% sampidx[ , si], ] 
	}else{
		DThr <- DTseam
	}

	for(bi in seq(1,Nsubsamp)){
		#subsamples: All, stayer, EE, EUE, All-EUUE, EUUE, 
		if(bi==1|bi==5){
			DThr[(changer==T|stayer==T), subIndic := T]
		}else if(bi==2){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[stayer==T, subIndic := T]
		}else if(bi==3){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[(changer==T & EE_wave==T)&!is.na(switchedOcc_wave), subIndic := T]
		}else if(bi==4| bi==6){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[changer==T & (EU_wave==T|UE_wave==T)&!is.na(switchedOcc_wave), subIndic := T]
		}
		if(bi==5|bi==6){
			wt = origwt
			if(wclab == "res"){
				wc = "wagechange_wave"
			}else{
				wc = "rawwgchange_wave"
			}
		}else{
			wt = "wtEUE"
			if(wclab=="res"){
					wc = "wagechangeEUE_wave"
			}else{
				wc = "rawwgchangeEUE_wave"
			}
		}
		
		tab_wavemoments[bi,1,1]      <- DThr[subIndic==T&demo==T,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
		tab_wavemoments[bi,1,2:5]    <- DThr[subIndic==T&demo==T,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		#expansion/recession
		for(rI in c(F,T)){
			rix = as.integer(rI)+2
			tab_wavemoments[bi,rix,1]      <- DThr[eval(as.name(recDef)) == rI & subIndic==T&demo==T,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
			tab_wavemoments[bi,rix,2:5]    <- DThr[eval(as.name(recDef)) == rI & subIndic==T&demo==T,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		}
	}
	if(si>0){
		se_wavemoments[,,,si] = tab_wavemoments
	}else{
		dat_wavemoments <- tab_wavemoments	
	}
} #si, simulation iteration

#output it to tables
outputtable <-array(NA, dim=c(Nt*Nmoments,Nsubsamp))
for(bi in seq(1,Nsubsamp)){
	for(mi in seq(0,Nmoments-1)){
		for(ti in seq(1,Nt)){
			outputtable[ti+ mi*Nt,bi] <- dat_wavemoments[bi,ti,mi+1]
		}
	}
}

outputtable <- data.table(outputtable)
names(outputtable) <- c("All","Stayers","EE", "EUE", "All-EU,UE","EU,UE")
rnames0 <- c("All\ Periods", "Expansion","Recession")
rnames1 <- rnames0
rnames <- rnames0
for(mi in seq(2,Nmoments)){
	for(ti in seq(1,Nt)){
		rnames1[ti] <- paste0(rnames1[ti],"\ ")
	}
	rnames<- c(rnames,rnames1)
}
rownames(outputtable) <- rnames
momentsnames <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
momentcommand <- paste0("\\hline  \\color{Maroon}{",momentsnames[1],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n")
for(mi in seq(2,Nmoments)){
	momentcommand <- c(momentcommand,paste0("\\hline  \\color{Maroon}{",momentsnames[mi],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n"))
}

rowtitles <- list( pos=as.list(seq(0,Nmoments-1)*Nt), command=momentcommand )

nametab <- "moments"

outputtable <- xtable(outputtable, digits=2, 
					   align="l|llll|ll", caption=paste0("Moments of earnings change distribution \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
if(demolbl>=1 & demolbl<=7){
	print(outputtable,include.rownames=T, hline.after= c(nrow(outputtable)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
}else{
	print(outputtable,include.rownames=T, hline.after= c(nrow(outputtable)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
}

tab_wavemomentsci <- array(NA,dim=c(2*nrow(outputtable),2*ncol(outputtable)))
tab_wavemomentsse <- array(0.,dim=c(2*nrow(outputtable),ncol(outputtable)))
if(bootse == T){

	for(bi in seq(1,Nsubsamp)){
		for(mi in seq(0,Nmoments-1)){
			for(ti in seq(0,Nt-1)){
				tab_wavemomentsse[(ti+ mi*Nt)*2+1,bi] <- dat_wavemoments[bi,ti+1,mi+1]
				tab_wavemomentsse[(ti+ mi*Nt)*2+2,bi] <- var(se_wavemoments[bi,ti+1,mi+1, ],na.rm=T)^.5
				tab_wavemomentsci[(ti+ mi*Nt)*2+1,(bi-1)*2+1] <- dat_wavemoments[bi,ti+1,mi+1]
				ci <- quantile(se_wavemoments[bi,ti+1,mi+1, ],na.rm=T,probs=c(0.025,0.975))
				tab_wavemomentsci[(ti+ mi*Nt)*2+2,(bi-1)*2+1]<-ci[1]
				tab_wavemomentsci[(ti+ mi*Nt)*2+2,bi*2]<-ci[2]
			}
		}
	}

	for(ri in seq(1,nrow(tab_wavemomentsci))){
		for(ci in seq(1,ncol(tab_wavemomentsse))){
			if(ci %%2==1){
				tab_wavemomentsse[ri,ci]<- round(tab_wavemomentsse[ri,ci],2)	
			}else{
				tab_wavemomentsse[ri,ci]<- round(tab_wavemomentsse[ri,ci],3)				
			}
		}
		for(ci in seq(1,ncol(tab_wavemomentsci))){
			if(ci %%2==1){
				tab_wavemomentsci[ri,ci]<- round(tab_wavemomentsci[ri,ci],2)	
			}else{
				tab_wavemomentsci[ri,ci]<- round(tab_wavemomentsci[ri,ci],3)
			}
		}
	}
	
	
	parens <- function(x, i, j){
		x[i,j] <- sprintf("(%s)", x[i,j])
		x
	}
	parensLR <- function(x, i, j){
		x[i,j] <- sprintf("(%s,%s)", x[i,j],x[i,j+1])
		x
	}
	for (ri in seq(2,nrow(tab_wavemomentsse),2)){
		for(ci in seq(1,ncol(tab_wavemomentsse))){
			tab_wavemomentsse <-parens(tab_wavemomentsse,ri,ci)
			tab_wavemomentsci <-parensLR(tab_wavemomentsci,ri,(ci-1)*2+1)
		}
	}
	tab_wavemomentsse <-data.table(tab_wavemomentsse)
	tab_wavemomentsci <-data.table(tab_wavemomentsci[,seq(1,ncol(tab_wavemomentsci),2)])
	
	names(tab_wavemomentsse) <- c("All","Stayers","EE", "EUE", "All-EU,UE","EU,UE")
	names(tab_wavemomentsci) <- c("All","Stayers","EE", "EUE", "All-EU,UE","EU,UE")
	
	rnames0 <- c("All\ Periods","% ","Expansion","% %","Recession","% % %")
	rnames1 <- rnames0
	rnames <- rnames0
	for(mi in seq(2,Nmoments)){
		for(ti in seq(1,Nt*2)){
			rnames1[ti] <- paste0(rnames1[ti],"\ ")
		}
		rnames<- c(rnames,rnames1)
	}
	rownames(tab_wavemomentsse) <- rnames
	rownames(tab_wavemomentsci) <- rnames
	momentsnames <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
	momentcommand <- paste0("\\hline  \\color{Maroon}{",momentsnames[1],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n")
	for(mi in seq(2,Nmoments)){
		momentcommand <- c(momentcommand,paste0("\\hline  \\color{Maroon}{",momentsnames[mi],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n"))
	}
	
	rowtitles <- list( pos=as.list(seq(0,Nmoments-1)*Nt*2), command=momentcommand )
	
	nametabse <- "momentsse"
	nametabci <- "momentsci"
	tab_wavemomentsse <- xtable(tab_wavemomentsse, digits=2, 
							align="l|llll|ll", caption=paste0("Moments of earnings change distribution \\label{tab:",nametabse,"_",wclab,"_",reclab,"}"))
	tab_wavemomentsci <- xtable(tab_wavemomentsci, digits=2, 
							align="l|cccc|cc", caption=paste0("Moments of earnings change distribution \\label{tab:",nametabci,"_",wclab,"_",reclab,"}"))
	
	if(demolbl>=1 & demolbl<=7){
		print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		print(tab_wavemomentsci,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		
	}else{
		print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,"_",wclab,"_",reclab,".tex"))
		print(tab_wavemomentsci,include.rownames=T, hline.after= c(nrow(tab_wavemomentsci)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabci,"_",wclab,"_",reclab,".tex"))
	}
	
}

#**********************************************************************************
#**********************************************************************************
# moments among job changers------------------------------------------------------
#**********************************************************************************
#**********************************************************************************

Nmoments <- 5 # the number of moemnts to spit out for each subset
Nsubsamp <- 6 # the number of population subsamples: EE-sw,EE-nosw, EUE-sw, EUE-nosw, EUUE-sw, EUUE-nosw
Nt <- 3 # time periods, all, exp, rec


if(bootse == T){
	set.seed(seedint)
	#draw the sample
	nsamp  = DTseam[ , max(idx)]
	se_chngmoments <- array(0.,dim = c(Nsubsamp,Nt,Nmoments,Nsim))
	sampidx <-array(1L, dim = c(nsamp,Nsim))
	for(si in seq(1,Nsim)){ sampidx[ , si] = sample(nsamp,nsamp,replace = T) }
}
tab_chngmoments <- array(0., dim=c(Nsubsamp,Nt,Nmoments))

for( si in seq(0,bootse*Nsim) ){
	if(si>0){
		DThr <- DTseam[idx %in% sampidx[ , si], ] 
	}else{
		DThr <- DTseam
	}
	
	for(bi in seq(1,Nsubsamp)){
		#subsamples: EE-sw,EE-nosw, EUE-sw, EUE-nosw, EUUE-sw, EUUE-nosw
		if(bi==1){
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[changer==T & EE_wave==T & switchedOcc_wave==T, subIndic := T]
		}else if(bi==2){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[changer==T & EE_wave==T & switchedOcc_wave==F, subIndic := T]
		}else if(bi==3| bi==5){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[changer==T &(EU_wave==T | UE_wave==T)& switchedOcc_wave==T, subIndic := T]
		}else if(bi==4| bi==6){
			DThr[, subIndic := NULL]
			DThr[(changer==T|stayer==T), subIndic := F]
			DThr[changer==T &(EU_wave==T | UE_wave==T)& switchedOcc_wave==F, subIndic := T]
		}
		if(bi==5|bi==6){
			wt = origwt
			if(wclab == "res"){
				wc = "wagechange_wave"
			}else{
				wc = "rawwgchange_wave"
			}
		}else{
			wt = "wtEUE"
			if(wclab=="res"){
				wc = "wagechangeEUE_wave"
			}else{
				wc = "rawwgchangeEUE_wave"
			}
		}
		
		tab_chngmoments[bi,1,1]      <- DThr[subIndic==T&demo==T,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
		tab_chngmoments[bi,1,2:5]    <- DThr[subIndic==T&demo==T,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		#expansion/recession
		for(rI in c(F,T)){
			rix = as.integer(rI)+2
			tab_chngmoments[bi,rix,1]      <- DThr[eval(as.name(recDef)) == rI & subIndic==T&demo==T,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
			tab_chngmoments[bi,rix,2:5]    <- DThr[eval(as.name(recDef)) == rI & subIndic==T&demo==T,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		}
	}
	if(si>0){
		se_chngmoments[,,,si] = tab_chngmoments
	}else{
		dat_chngmoments <- tab_chngmoments	
	}
} #si, simulation iteration


#output it to tables
outputtable <-array(NA, dim=c(Nt*Nmoments,Nsubsamp))
for(bi in seq(1,Nsubsamp)){
	for(mi in seq(0,Nmoments-1)){
		for(ti in seq(1,Nt)){
			outputtable[ti+ mi*Nt,bi] <- dat_chngmoments[bi,ti,mi+1]
		}
	}
}

outputtable <- data.table(outputtable)
names(outputtable) <- c("EE","EE\ ","EUE","EUE\ ","EU,UE","EU,UE\ ")
title0 <- c( "& \\multicolumn{2}{|c}{EE} & \\multicolumn{2}{|c|}{EUE} & \\multicolumn{2}{|c|}{EU,UE}  \\\\ \n", 
			"&   Switch Occ & No Switch & Switch Occ & No Switch & Switch Occ & No Switch  \\\\ \\hline \n")
rnames0 <- c("All\ Periods", "Expansion","Recession")
rnames1 <- rnames0
rnames <- rnames0
for(mi in seq(2,Nmoments)){
	for(ti in seq(1,Nt)){
		rnames1[ti] <- paste0(rnames1[ti],"\ ")
	}
	rnames<- c(rnames,rnames1)
}
rownames(outputtable) <- rnames
momentsnames <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
momentcommand <- title0 #paste0("\\hline  \\color{Maroon}{",momentsnames[1],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n")
for(mi in seq(1,Nmoments)){
	momentcommand <- c(momentcommand,paste0("\\hline  \\color{Maroon}{",momentsnames[mi],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n"))
}

rowtitles <- list( pos=as.list(c(0,0,seq(0,Nmoments-1)*Nt)), command=momentcommand )

nametab <- "chngmoments"

outputtable <- xtable(outputtable, digits=2, 
					  align="l|llll|ll", caption=paste0("Moments of earnings change distribution \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
if(demolbl>=1 & demolbl<=7){
	print(outputtable,include.rownames=T, include.colnames = F, hline.after= c(nrow(outputtable)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
}else{
	print(outputtable,include.rownames=T, include.colnames = F, hline.after= c(nrow(outputtable)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
}

tab_wavemomentsci <- array(NA,dim=c(2*nrow(outputtable),2*ncol(outputtable)))
tab_wavemomentsse <- array(0.,dim=c(2*nrow(outputtable),ncol(outputtable)))
if(bootse == T){
	
	for(bi in seq(1,Nsubsamp)){
		for(mi in seq(0,Nmoments-1)){
			for(ti in seq(0,Nt-1)){
				tab_wavemomentsse[(ti+ mi*Nt)*2+1,bi] <- dat_chngmoments[bi,ti+1,mi+1]
				tab_wavemomentsse[(ti+ mi*Nt)*2+2,bi] <- var(se_chngmoments[bi,ti+1,mi+1, ],na.rm=T)^.5
				tab_wavemomentsci[(ti+ mi*Nt)*2+1,(bi-1)*2+1] <- dat_chngmoments[bi,ti+1,mi+1]
				ci <- quantile(se_chngmoments[bi,ti+1,mi+1, ],na.rm=T,probs=c(0.025,0.975))
				tab_wavemomentsci[(ti+ mi*Nt)*2+2,(bi-1)*2+1]<-ci[1]
				tab_wavemomentsci[(ti+ mi*Nt)*2+2,bi*2]<-ci[2]
			}
		}
	}

	#round everything:
	for(ri in seq(1,nrow(tab_wavemomentsse))){
		for(ci in seq(1,ncol(tab_wavemomentsse))){
			if(ci %% 2 ==1){
				tab_wavemomentsse[ri,ci] = round(tab_wavemomentsse[ri,ci],2)
			}else{
				tab_wavemomentsse[ri,ci] = round(tab_wavemomentsse[ri,ci],3)	
			}
		}
		for(ci in seq(1,ncol(tab_wavemomentsci))){
			if(ci %% 2 ==1){
				tab_wavemomentsci[ri,ci] = round(tab_wavemomentsci[ri,ci],2)
			}else{
				tab_wavemomentsci[ri,ci] = round(tab_wavemomentsci[ri,ci],3)	
			}
			
		}
	}
	parens <- function(x, i, j){
		x[i,j] <- sprintf("(%s)", x[i,j])
		x
	}
	parensLR <- function(x, i, j){
		x[i,j] <- sprintf("(%s,%s)", x[i,j],x[i,j+1])
		x
	}
	for (ri in seq(2,nrow(tab_wavemomentsse),2)){
		for(ci in seq(1,ncol(tab_wavemomentsse))){
			tab_wavemomentsse <-parens(tab_wavemomentsse,ri,ci)
			tab_wavemomentsci <-parensLR(tab_wavemomentsci,ri,(ci-1)*2+1)
		}
	}
	tab_wavemomentsse <-data.table(tab_wavemomentsse)
	tab_wavemomentsci <-data.table(tab_wavemomentsci[,seq(1,ncol(tab_wavemomentsci),2)])
	
	
	names(tab_wavemomentsse) <- c("EE","EE\ ","EUE","EUE\ ","EU,UE","EU,UE\ ")
	names(tab_wavemomentsci) <- c("EE","EE\ ","EUE","EUE\ ","EU,UE","EU,UE\ ")
	title0 <- c( "& \\multicolumn{2}{|c|}{EE} & \\multicolumn{2}{|c|}{EUE} & \\multicolumn{2}{|c|}{EU,UE}  \\\\ \n", 
				 "&   Swtich Occ & No Switch & Swtich Occ & No Switch & Swtich Occ & No Switch  \\\\ \\hline \n")
	rnames0 <- c("All\ Periods","~","Expansion","~ ~","Recession","~ ~ ~")
	rnames1 <- rnames0
	rnames <- rnames0
	for(mi in seq(2,Nmoments)){
		for(ti in seq(1,Nt*2)){
			rnames1[ti] <- paste0(rnames1[ti],"\ ")
		}
		rnames<- c(rnames,rnames1)
	}
	rownames(tab_wavemomentsse) <- rnames
	rownames(tab_wavemomentsci) <- rnames
	momentsnames <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
	momentcommand <- title0 #paste0("\\hline  \\color{Maroon}{",momentsnames[1],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n")
	for(mi in seq(1,Nmoments)){
		momentcommand <- c(momentcommand,paste0("\\hline  \\color{Maroon}{",momentsnames[mi],"} &  \\multicolumn{",as.character(Nsubsamp),"}{|c|}{} \\\\ \n"))
	}
	
	rowtitles <- list( pos=as.list(c(0,0,seq(0,Nmoments-1)*Nt*2)), command=momentcommand )
	
	nametabse <- "chngmomentsse"
	nametabci <- "chngmomentsci"
	tab_wavemomentsse <- xtable(tab_wavemomentsse, digits=2, 
								align="l|llll|ll", caption=paste0("Moments of earnings change distribution \\label{tab:",nametabse,"_",wclab,"_",reclab,"}"))
	tab_wavemomentsci <- xtable(tab_wavemomentsci, digits=2, 
								align="l|cccc|cc", caption=paste0("Moments of earnings change distribution \\label{tab:",nametabci,"_",wclab,"_",reclab,"}"))
	
	if(demolbl>=1 & demolbl<=7){
		print(tab_wavemomentsse,include.rownames=T, include.colnames=F, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		print(tab_wavemomentsci,include.rownames=T, include.colnames=F, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		
	}else{
		print(tab_wavemomentsse,include.rownames=T, include.colnames=F, hline.after= c(nrow(tab_wavemomentsse)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabse,"_",wclab,"_",reclab,".tex"))
		print(tab_wavemomentsci,include.rownames=T, include.colnames=F, hline.after= c(nrow(tab_wavemomentsci)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametabci,"_",wclab,"_",reclab,".tex"))
	}
	
}


#Variance decomp -----------------------------------------------

tab_wavevardec <- array(NA_real_,dim=c(2,3))

totmean <- DTseam[(stayer|changer)&demo==T, wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt))) ]
totvar  <- DTseam[(stayer|changer)&demo==T, sum(eval(as.name(wt))*(eval(as.name(wc))- totmean)^2,na.rm=T) ]
tab_wavevardec[1,1] <- DTseam[stayer == T                        &demo==T, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,2] <- DTseam[changer ==T& EE_wave ==T           &demo==T, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,3] <- DTseam[changer ==T&(EU_wave==T|UE_wave==T)&demo==T, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
totwt <- DTseam[(stayer|changer)&demo==T, sum(eval(as.name(wt)),na.rm=T) ]
tab_wavevardec[2,1] <- DTseam[stayer == T                        &demo==T, sum(eval(as.name(wt)),na.rm=T) ]/totwt
tab_wavevardec[2,2] <- DTseam[changer ==T& EE_wave ==T           &demo==T, sum(eval(as.name(wt)),na.rm=T) ]/totwt
tab_wavevardec[2,3] <- DTseam[changer ==T&(EU_wave==T|UE_wave==T)&demo==T, sum(eval(as.name(wt)),na.rm=T) ]/totwt


# Full sample quantile-diff decomposition --------------------------
tot51025qtl <- DTseam[(stayer==T|changer==T)&demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=c(0.05,0.1,.25,.75,.9,.95)) ]

Nqtls <-length(tot51025qtl)
Ndifs <- Nqtls/2

tab_waveqtldec <- array(NA_real_,dim=c(Ndifs+1,3))

qtlshere <- tot51025qtl

for(rri in seq(1,Ndifs)){
  tab_waveqtldec[rri,1] <- DTseam[stayer ==T                          &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  tab_waveqtldec[rri,2] <- DTseam[EE_wave==T & changer==T             &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  tab_waveqtldec[rri,3] <- DTseam[(EU_wave==T|UE_wave==T)& changer==T &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
}
tab_waveqtldec[Ndifs+1,1] <- DTseam[stayer == T                          &demo==T, sum(eval(as.name(wt)),na.rm=T) ]
tab_waveqtldec[Ndifs+1,2] <- DTseam[EE_wave==T & changer==T              &demo==T, sum(eval(as.name(wt)),na.rm=T) ]
tab_waveqtldec[Ndifs+1,3] <- DTseam[(EU_wave==T|UE_wave==T)& changer==T  &demo==T, sum(eval(as.name(wt)),na.rm=T) ]


rsum <- rowSums(tab_waveqtldec, na.rm=T)
for(ri in seq(1,nrow(tab_waveqtldec))){
	tab_waveqtldec[ri,] <-tab_waveqtldec[ri,]/rsum[ri]
}

#output it to tables
tab_chngvarqtldec <- data.table(rbind(tab_wavevardec[1,],tab_waveqtldec[1:Ndifs,],tab_wavevardec[2,]) )
if(wc=="wagechange_wave"|wc=="rawwgchange_wave"){
	names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EU,UE")	
}else{
	names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EUE")	
}


#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngvarqtldec) <- c("Variance","0.95-0.05","0.9-0.1","0.75-0.25","Pop")# "Variance\ ","0.95-0.05\ ","0.9-0.1\ ","0.75-0.25\ ","Pct Sample")

tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
                            align="l|l|ll", caption="Decomposition of earnings change dispersion \\label{tab:wavechngvarqtldec}")

nametab <- "wavechngvarqtldec"

if(demolbl>=1 & demolbl<=7){
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/",nametab, demotxt[demolbl],"_",wclab,"_",reclab,".tex")) 
}else{
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex") ) 
}

# recession and expansion
for (rI in c(F,T)){
  totmean <- DTseam[recIndic_wave==rI, wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) ) ]
  totvar  <- DTseam[recIndic_wave==rI, sum(eval(as.name(wt))*(eval(as.name(wc))- totmean)^2,na.rm=T) ]
  tab_wavevardec[1,1] <- DTseam[(EE_wave==F&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
  tab_wavevardec[1,2] <- DTseam[(EE_wave==T&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
  tab_wavevardec[1,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T))&recIndic_wave==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
  totwt <- DTseam[recIndic_wave==rI, sum(eval(as.name(wt)),na.rm=T) ]
  tab_wavevardec[2,1] <- DTseam[(EE_wave==F&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt
  tab_wavevardec[2,2] <- DTseam[(EE_wave==T&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt
  tab_wavevardec[2,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T))&recIndic_wave==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt

  # Seam sample quantile-diff decomposition --------------------------
  tot51025qtl <- DTseam[ recIndic_wave==rI, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=c(0.05,0.1,.25,.75,.9,.95)) ]
  
  Nqtls <-length(tot51025qtl)
  Ndifs <- Nqtls/2
  
  tab_waveqtldec <- array(NA_real_,dim=c(Ndifs+1,3))
  
  qtlshere <- tot51025qtl
  for(rri in seq(1,Ndifs)){
    tab_waveqtldec[rri,1] <- DTseam[recIndic_wave==rI & (EE_wave==F&EU_wave==F&UE_wave==F)   &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
    tab_waveqtldec[rri,2] <- DTseam[recIndic_wave==rI & (EE_wave==T&EU_wave==F&UE_wave==F)   &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
    tab_waveqtldec[rri,3] <- DTseam[recIndic_wave==rI & (EE_wave==F&(EU_wave==T|UE_wave==T)) &demo==T& (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  }
  tab_waveqtldec[Ndifs+1,1] <- DTseam[recIndic_wave==rI &  (stayer==T)                        &demo==T, sum(eval(as.name(wt)),na.rm=T) ]
  tab_waveqtldec[Ndifs+1,2] <- DTseam[recIndic_wave==rI &  (EE_wave==T&changer==T)            &demo==T, sum(eval(as.name(wt)),na.rm=T) ]
  tab_waveqtldec[Ndifs+1,3] <- DTseam[recIndic_wave==rI & ((EU_wave==T|UE_wave==T)&changer==T)&demo==T, sum(eval(as.name(wt)),na.rm=T) ]

  
  rsum <- rowSums(tab_waveqtldec, na.rm=T)
  for(ri in seq(1,nrow(tab_waveqtldec))){
    tab_waveqtldec[ri,] <-tab_waveqtldec[ri,]/rsum[ri]
  }
  
  #output it to tables
  tab_chngvarqtldec <- data.table(rbind(tab_wavevardec[1,],tab_waveqtldec[1:Ndifs,]) )
  names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EU,UE")
  #rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
  rownames(tab_chngvarqtldec) <- c("Variance","0.95-0.05","0.9-0.1","0.75-0.25")
  
  if(rI){
    tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
                                align="l|l|ll", caption="Decomposition of earnings change dispersion, recession \\label{tab:wavechngvarqtldec_rec}")
  }else{
    tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
                                align="l|l|ll", caption="Decomposition of earnings change dispersion, expansion \\label{tab:wavechngvarqtldec_exp}")
  }
  if(demolbl>=1L & demolbl<=7L){ 
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/wavechngvarqtldec_rec",rI,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
  }else{
  	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/wavechngvarqtldec_rec",rI,"_",wclab,"_",reclab,".tex"))
  }
}

#*******************************************************************
#*******************************************************************
# Moments among job changers  -----------------------
#tabqtls <- c(.1,.25,.5,.75,.9)
#tN <- (length(tabqtls)+1)


if(bootse==T){
	set.seed(seedint)
	#draw the sample
	Nsim = 50
	nsampE = nrow(DTseam[eval(as.name(recDef)) == F ])
	nsampR = nrow(DTseam[eval(as.name(recDef)) == T ])
	nsamp  = nsampR+nsampE
	se_wavechngdist <- array(0.,dim = c(9,length(tabqtls)+1,Nsim))
}else{
	Nsim = 0
}

for(si in seq(1,bootse*Nsim+1)){
	if(si>1){
		DThr <- DTseam[ sample(nsamp,nsamp,replace=T)]#,prob=eval(as.name(wt))
	}else{
		DThr <- DTseam
	}
	tab_wavechngdist[1,1]    <- DThr[changer==T&demo==T ,                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[1,2:tN] <- DThr[changer==T&demo==T ,                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_wavechngdist[2,1]    <- DThr[changer==T&demo==T & switchedOcc_wave==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[2,2:tN] <- DThr[changer==T&demo==T & switchedOcc_wave==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavechngdist[3,1]    <- DThr[changer==T&demo==T & switchedOcc_wave==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[3,2:tN] <- DThr[changer==T&demo==T & switchedOcc_wave==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	for(EUhere in c(F,T)){ #loop over EE or EU,UE
		EEhere = (EUhere==F)
		eidx = as.integer(EUhere)*3
		tab_wavechngdist[4+eidx,1]    <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere)),                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist[4+eidx,2:tN] <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere)),                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
		tab_wavechngdist[5+eidx,1]    <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist[5+eidx,2:tN] <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_wavechngdist[6+eidx,1]    <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist[6+eidx,2:tN] <- DThr[changer==T&demo==T & (EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}
	if(si>1){
		se_wavechngdist[,,si-1] = tab_wavechngdist
	}else{
		dat_wavechngdist <- tab_wavechngdist
	}
	
}
#output it to tables
tab_wavechngdist <- data.table(dat_wavechngdist)
names(tab_wavechngdist) <- c("Mean",tabqtls)
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
if(wc =="wagechange_wave" | wc=="rawwgchange_wave"){
rownames(tab_wavechngdist) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job, Switch\ Occ",
									"EE,\ All", "EE,\ Same\ Occ", "EE,\ Switch\ Occ",
									"EU,UE,\ All", "EU,UE,\ Same\ Occ", "EU,UE,\ Switch\ Occ")
}else if(wc == "wagechangeEUE_wave"| wc=="rawwgchangeEUE_wave"){
rownames(tab_wavechngdist) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job, Switch\ Occ",
									"EE,\ All", "EE,\ Same\ Occ", "EE,\ Switch\ Occ",
									"EUE,\ All", "EUE,\ Same\ Occ", "EUE,\ Switch\ Occ")
}

labtxt <- "chngdist"

tab_wavechngdist <- xtable(tab_wavechngdist, label=paste0("tab:",labtxt), digits=2, 
                           align="l|l|lllll", caption="Distribution of earnings changes among job changers")
if(demolbl>=1 & demolbl<=7){
	print(tab_wavechngdist,include.rownames=T, hline.after= c(0,3,6,nrow(tab_wavechngdist)), file=paste0(outputdir,"/",labtxt,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
}else{
	print(tab_wavechngdist,include.rownames=T, hline.after= c(0,3,6,nrow(tab_wavechngdist)), file=paste0(outputdir,"/",labtxt,"_",wclab,"_",reclab,".tex"))
}
# standard errors
tab_wavechngdistci <- array(0.,dim=c(2*nrow(tab_wavechngdist),2*ncol(tab_wavechngdist)))
tab_wavechngdistse <- array(0.,dim=c(2*nrow(tab_wavechngdist),ncol(tab_wavechngdist)))
if(bootse == T){
	for( ri in seq(0,nrow(tab_wavechngdist)-1) ){
		for(ci in seq( 0,ncol(tab_wavechngdist)-1 )){
			tab_wavechngdistci[ ri*2+1,ci*2+1 ] <- tab_wavechngdist[ri+1,ci+1]
			tab_wavechngdistci[ ri*2+2,(ci*2+1):(ci*2+2) ] <- quantile(se_wavechngdist[ri+1,ci+1, ], probs=c(0.05,0.95))
			tab_wavechngdistse[ ri*2+1,ci+1 ] <- tab_wavechngdist[ri+1,ci+1]
			tab_wavechngdistse[ ri*2+2,ci+1 ] <- var(se_wavechngdist[ri+1,ci+1, ])^0.5
		}
	}
	
	# output to tables with standard errors
	#pre-round it:
	tab_wavechngdistse <- round(100*tab_wavechngdistse,digits=2)
	tab_wavechngdistci <- round(100*tab_wavechngdistci,digits=2)
	
	parens <- function(x, i, j){
		x[i,j] <- sprintf("(%s)", x[i,j])
		x
	}
	for (ri in seq(2,nrow(tab_wavechngdistse),2)){
		for(ci in seq(1,ncol(tab_wavechngdistse))){
			tab_wavechngdistse <-parens(tab_wavechngdistse,ri,ci)
		}
	}
	tab_wavechngdistse <-data.table(tab_wavechngdistse)
	
	names(tab_wavechngdistse) <- c("Mean",as.character( tabqtls))
	if(wc == "wagechangeEUE_wave"){
		rnames <- c("Chng\ Job, All"," "     , "Chng\ Job,\ Same\ Occ",","   , "Chng\ Job,\ Switch\ Occ",",,",
					"EE,\ All"      ,",,,"   , "EE,\ Same\ Occ"       ,",,,,", "EE,\ Switch\ Occ"       ,",,,,,",
					"EUE,\ All"     ,",,,,,,", "EUE,\ Same\ Occ"   ,",,,,,,,", "EUE\ Switch\ Occ"       ,",,,,,,,,")
	}else{
		rnames <- c("Chng\ Job, All"," "     , "Chng\ Job,\ Same\ Occ",","   , "Chng\ Job,\ Switch\ Occ",",,",
					"EE,\ All"      ,",,,"   , "EE,\ Same\ Occ"       ,",,,,", "EE,\ Switch\ Occ"       ,",,,,,",
					"EU,UE,\ All"   ,",,,,,,", "EU,UE,\ Same\ Occ" ,",,,,,,,", "EU,UE\ Switch\ Occ"     ,",,,,,,,,")
	}
	rownames(tab_wavechngdistse) <- rnames
	
	nametab <- "chngdistse"
	
	tab_wavechngdistse <- xtable(tab_wavechngdistse, digits=2, 
						   align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
	
	if(demolbl>=1 & demolbl<=7){
		print(tab_wavechngdistse,include.rownames=T, hline.after= c(0,6,12,nrow(tab_wavechngdistse)), 
			  file=paste0(outputdir,"/",nametab,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
	}else{
		print(tab_wavechngdistse,include.rownames=T, hline.after= c(0,6,12,nrow(tab_wavechngdistse)), 
			  file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
	}
}
#*****************************************************************
# rec and expansion, job changers


tab_wavechngdist_rec    <- array(NA,dim=c(9*3,tN))
tab_wavechngmoments_rec     <- array(NA, dim=c(9*3,5))
if(bootse == T){
	set.seed(seedint)
	#draw the sample
	Nsim = 50
	nsampE = nrow(DTseam[eval(as.name(recDef)) == F ])
	nsampR = nrow(DTseam[eval(as.name(recDef)) == T ])
	nsamp  = nsampR+nsampE
	se_wavechngdist_rec <- array(0.,dim = c(9*3,length(tabqtls)+1,Nsim))
	se_wavechngmoments_rec  <- array(0.,dim = c(9,5,Nsim))
}else{
	Nsim = 0
}

for( si in seq(1,bootse*Nsim+1) ){
	if(si>1){
		DThr <- DTseam[ sample(nsamp,nsamp,replace=T)] #,prob=eval(as.name(wt))
	}else{
		DThr <- DTseam
	}
	for(AllEEEU in c(1,2,3)){
		if(AllEEEU==1){
			DThr[demo==T, trX_indic := T]
		}else if(AllEEEU==2){
			DThr[ , trX_indic := NULL]
			DThr[                     , trX_indic := F]
			DThr[ EE_wave==T & demo==T, trX_indic := T]
		}else if(AllEEEU==3){
			DThr[ , trX_indic := NULL]
			DThr[                                , trX_indic := F]
			DThr[ (EU_wave==T|UE_wave==T)&demo==T, trX_indic := T]
		}
		eidx = 9*(AllEEEU-1)

		tab_wavechngdist_rec[1+eidx,1   ] <- DThr[changer==T &trX_indic==T                      ,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec[1+eidx,2:tN] <- DThr[changer==T &trX_indic==T                      ,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
		tab_wavechngdist_rec[2+eidx,1   ] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F ,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec[2+eidx,2:tN] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F ,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_wavechngdist_rec[3+eidx,1   ] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T ,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec[3+eidx,2:tN] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T ,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		for(rI in c(F,T)){
			ridx = as.integer(rI)*3+3
			tab_wavechngdist_rec[1+ridx+eidx,1   ] <- DThr[changer==T &trX_indic==T                      &eval(as.name(recDef))  == rI,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
			tab_wavechngdist_rec[1+ridx+eidx,2:tN] <- DThr[changer==T &trX_indic==T                      &eval(as.name(recDef))  == rI,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
			tab_wavechngdist_rec[2+ridx+eidx,1   ] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F &eval(as.name(recDef))  == rI,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
			tab_wavechngdist_rec[2+ridx+eidx,2:tN] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F &eval(as.name(recDef))  == rI,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
			tab_wavechngdist_rec[3+ridx+eidx,1   ] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T &eval(as.name(recDef))  == rI,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
			tab_wavechngdist_rec[3+ridx+eidx,2:tN] <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T &eval(as.name(recDef))  == rI,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		}
		#wage moments:
		tab_wavechngmoments_rec[1+eidx,1]    <- DThr[changer==T &trX_indic==T &!is.na(switchedOcc_wave),       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
		tab_wavechngmoments_rec[1+eidx,2:5]  <- DThr[changer==T &trX_indic==T &!is.na(switchedOcc_wave),wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		tab_wavechngmoments_rec[2+eidx,1]    <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F ,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngmoments_rec[2+eidx,2:5]  <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F ,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		tab_wavechngmoments_rec[3+eidx,1]    <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T ,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngmoments_rec[3+eidx,2:5]  <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T ,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
	
		#expansion/recession
		for(rI in c(F,T)){
			ridx = rI*3+3
			tab_wavechngmoments_rec[1+eidx+ridx,1]    <- DThr[changer==T &trX_indic==T &!is.na(switchedOcc_wave) &eval(as.name(recDef))  == rI,        wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
			tab_wavechngmoments_rec[1+eidx+ridx,2:5]  <- DThr[changer==T &trX_indic==T &!is.na(switchedOcc_wave) &eval(as.name(recDef))  == rI,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
			tab_wavechngmoments_rec[2+eidx+ridx,1]    <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F &eval(as.name(recDef))  == rI,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
			tab_wavechngmoments_rec[2+eidx+ridx,2:5]  <- DThr[changer==T &trX_indic==T &switchedOcc_wave==F &eval(as.name(recDef))  == rI,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
			tab_wavechngmoments_rec[3+eidx+ridx,1]    <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T &eval(as.name(recDef))  == rI,       wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
			tab_wavechngmoments_rec[3+eidx+ridx,2:5]  <- DThr[changer==T &trX_indic==T &switchedOcc_wave==T &eval(as.name(recDef))  == rI,wtd.4qtlmoments(eval(as.name(wc)),eval(as.name(wt)))]
		}
	}
	if(si>1){
		se_wavechngdist_rec[,,si-1] = tab_wavechngdist_rec
		se_wavechngmoments_rec[,,si-1] = tab_wavechngmoments_rec
	}else{
		dat_wavechngdist_rec <- tab_wavechngdist_rec
		dat_wavechngmoments_rec <- tab_wavechngmoments_rec
	}
}#si, sim
#output it to tables
tab_wavechngdist_rec <- dat_wavechngdist_rec
tab_wavechngmoments_rec<-dat_wavechngmoments_rec

for(AllEEEU in c(1,2,3)){
	eidx = 9*(AllEEEU-1)
	tab_wavedist <- data.table(tab_wavechngdist_rec[seq(1+eidx,9+eidx),])

	labtxt = "recchngdist"

	if(AllEEEU==1){
		rtxt <- "All\ changers"
		labtxt <- paste0(labtxt,"_all")
	}else if(AllEEEU==2){
		rtxt <- "EE"
		labtxt <- paste0(labtxt,"_EE")
	}else if(AllEEEU==3){
		if( wc=="wagechangeEUE_wave" | wc=="rawwgchangeEUE_wave"){
			rtxt <- "EUE"
			labtxt <- paste0(labtxt,"_EUE")
		}else{
			rtxt <- "EU,UE"
			labtxt <- paste0(labtxt,"_EUUE")
		}
	}
	cnames <- c("Mean",tabqtls)
	rnames <- c(paste0(rtxt,"    "), "Occ stayers    ","Occ movers    ",
		        paste0(rtxt,"\   "), "Occ stayers\   ","Occ movers\   ",
		        paste0(rtxt,"\ \ "), "Occ stayers\ \ ","Occ movers\ \ ")

	rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
		                                          "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
		                                          "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )

	names(tab_wavedist) <- cnames
	rownames(tab_wavedist) <- rnames
	tab_wavedist <- xtable(tab_wavedist, digits=2, 
		                       align="l|l|lllll", caption=paste0("Distribution of earnings changes among ", rtxt ," job changers in recession and expansion \\label{tab:",labtxt,"_",wclab,"_",reclab,"}"))
	print(tab_wavedist,include.rownames=T, hline.after= c(nrow(tab_wavedist)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",labtxt,"_",wclab,"_",reclab,".tex"))

	#standard errors
	tab_wavedistci <- array(0.,dim=c(2*nrow(tab_wavedist),2*ncol(tab_wavedist)))
	tab_wavedistse <- array(0.,dim=c(2*nrow(tab_wavedist),ncol(tab_wavedist)))
	if(bootse == T){
		for( ri in seq(0,nrow(tab_wavedist)-1) ){
			for(ci in seq( 0,ncol(tab_wavedist)-1 )){
				tab_wavedistci[ ri*2+1,ci*2+1 ] <- tab_wavedist[ri+1,ci+1]
				tab_wavedistci[ ri*2+2,(ci*2+1):(ci*2+2) ] <- quantile(se_wavechngdist_rec[eidx+ri+1,ci+1, ], probs=c(0.05,0.95))
				tab_wavedistse[ ri*2+1,ci+1 ] <- tab_wavedist[ri+1,ci+1]
				tab_wavedistse[ ri*2+2,ci+1 ] <- var(se_wavechngdist_rec[eidx+ri+1,ci+1, ])^0.5
			}
		}
		#output it to tables
	
		#pre-round it:
		tab_wavedistse <- round(100*tab_wavedistse,digits=2)
		tab_wavedistci <- round(100*tab_wavedistci,digits=2)
	
		parens <- function(x, i, j){
			x[i,j] <- sprintf("(%s)", x[i,j])
			x
		}
		for (ri in seq(2,nrow(tab_wavedistse),2)){
			for(ci in seq(1,ncol(tab_wavedistse))){
				tab_wavedistse <-parens(tab_wavedistse,ri,ci)
			}
		}
		tab_wavedistse <-data.table(tab_wavedistse)
	
		names(tab_wavedistse) <- c("Mean",as.character( tabqtls))
		#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
		rnames <- c(paste0(rtxt,"    ")," "     , "Occ stayers    ",","      ,"Occ movers    ",
			        paste0(rtxt,"\   "),",,,"   , "Occ stayers\   ",",,,,"   ,"Occ movers\   ",
			        paste0(rtxt,"\ \ "),",,,,,,", "Occ stayers\ \ ",",,,,,,,","Occ movers\ \ ")
	
		rownames(tab_wavedistse) <- rnames
	
		rowtitles <- list( pos=list(0,6,12), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
													  "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
													  "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
		labtxt <- "recchngdistse"
		if(AllEEEU==1){
			rtxt <- "All\ changers"
			labtxt <- paste0(labtxt,"_all")
		}else if(AllEEEU==2){
			rtxt <- "EE"
			labtxt <- paste0(labtxt,"_EE")
		}else if(AllEEEU==3){
			if(wc=="wagechangeEUE_wave"| wc=="rawwgchangeEUE_wave"){
				rtxt <- "EUE"
				labtxt <- paste0(labtxt,"_EUE")
			}else{	
				rtxt <- "EU,UE"
				labtxt <- paste0(labtxt,"_EUUE")
			}
		}
		tab_wavedistse <- xtable(tab_wavedistse, digits=2, 
							   align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",labtxt,"_",wclab,"_",reclab,"}"))
	
		if(demolbl>=1 & demolbl<=7){
			print(tab_wavedistse,include.rownames=T, hline.after= c(nrow(tab_wavedistse)), 
				  add.to.row=rowtitles, file=paste0(outputdir,"/",labtxt,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		}else{
			print(tab_wavedistse,include.rownames=T, hline.after= c(nrow(tab_wavedistse)), 
				  add.to.row=rowtitles, file=paste0(outputdir,"/",labtxt,"_",wclab,"_",reclab,".tex"))
		}
	}# se on dist stats?
	
	# print the moments tables
	tab_wavemoments <- data.table(tab_wavechngmoments_rec[seq(1+eidx,9+eidx),])
	
	#names(tab_wavemoments) <- c("Mean","Median","Std Dev", "Skew", "Kurtosis")
	names(tab_wavemoments) <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
	#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
	rnames <- c(paste0(rtxt,"    "), "Occ stayers    ","Occ movers    ",
				paste0(rtxt,"\   "), "Occ stayers\   ","Occ movers\   ",
				paste0(rtxt,"\ \ "), "Occ stayers\ \ ","Occ movers\ \ ")
	
	rownames(tab_wavemoments) <- rnames
	
	rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  \\multicolumn{5}{|c|}{} \\\\ \n",
												  "\\hline \\hline   \\color{Maroon}{Expansion} &  \\multicolumn{5}{|c|}{}  \\\\  \n", 
												  "\\hline \\hline   \\color{Maroon}{Recession} &  \\multicolumn{5}{|c|}{}  \\\\  \n")  )
	
	nametab <- "recchngmoments"
	if(AllEEEU==1){
		rtxt <- "All\ changers"
		nametab <- paste0(nametab,"_all")
	}else if(AllEEEU==2){
		rtxt <- "EE"
		nametab <- paste0(nametab,"_EE")
	}else if(AllEEEU==3){
		if(wc=="wagechangeEUE_wave"| wc=="rawwgchangeEUE_wave"){
			rtxt <- "EUE"
			nametab <- paste0(nametab,"_EUE")
		}else{	
			rtxt <- "EU,UE"
			nametab <- paste0(nametab,"_EUUE")
		}
	}
	tab_wavemoments <- xtable(tab_wavemoments, digits=2, 
							  align="l|lllll", caption=paste0("Moments of  ", rtxt ," job changers' earnings change distribution \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
	if(demolbl>=1 & demolbl<=7){
		print(tab_wavemoments,include.rownames=T, hline.after= c(nrow(tab_wavemoments)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
	}else{
		print(tab_wavemoments,include.rownames=T, hline.after= c(nrow(tab_wavemoments)), 
			  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
	}
	tab_wavemomentsci <- array(0.,dim=c(2*nrow(tab_wavemoments),2*ncol(tab_wavemoments)))
	tab_wavemomentsse <- array(0.,dim=c(2*nrow(tab_wavemoments),ncol(tab_wavemoments)))
	if(bootse == T){
		for( ri in seq(0,nrow(tab_wavemoments)-1) ){
			for(ci in seq( 0,ncol(tab_wavemoments)-1 )){
				tab_wavemomentsci[ ri*2+1,ci*2+1 ] <- tab_wavechngmoments[ri+1,ci+1]
				tab_wavemomentsci[ ri*2+2,(ci*2+1):(ci*2+2) ] <- quantile(se_wavechngmoments[ri+1,ci+1, ], probs=c(0.05,0.95))
				tab_wavemomentsse[ ri*2+1,ci+1 ] <- tab_wavechngmoments[ri+1,ci+1]
				tab_wavemomentsse[ ri*2+2,ci+1 ] <- var(se_wavechngmoments[ri+1,ci+1, ])^0.5
			}
		}
		
		#output it to tables
		tab_wavemomentsse[,1:3]<-tab_wavemomentsse[,1:3]*100
		tab_wavemomentsci[,1:3]<-tab_wavemomentsci[,1:3]*100
		
		parens <- function(x, i, j){
			x[i,j] <- sprintf("(%s)", x[i,j])
			x
		}
		for (ri in seq(2,nrow(tab_wavemomentsse),2)){
			for(ci in seq(1,ncol(tab_wavemomentsse))){
				tab_wavemomentsse <-parens(tab_wavemomentsse,ri,ci)
			}
		}
		tab_wavemomentsse <-data.table(tab_wavemomentsse)
		
		names(tab_wavemomentsse) <- c("Mean",as.character( tabqtls))
		#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
		rnames <- c("All\ Workers"    ," "     , "Occ\ Stayers"    ,","      ,"Occ\ Movers"    ,",,",
					"All\ Workers\ "  ,",,,"   , "Occ\ Stayers\ "  ,",,,,"   ,"Occ\ Movers\ "  ,",,,,,",
					"All\ Workers\ \ ",",,,,,,", "Occ\ Stayers\ \ ",",,,,,,,","Occ\ Movers\ \ ",",,,,,,,,")
		rownames(tab_wavemomentsse) <- rnames
		
		rowtitles <- list( pos=list(0,6,12), command=c("\\hline  \\color{Maroon}{1996-2012} &  \\multicolumn{5}{|c|}{} \\\\ \n",
													   "\\hline \\hline   \\color{Maroon}{Expansion} &  \\multicolumn{5}{|c|}{}  \\\\  \n", 
													   "\\hline \\hline   \\color{Maroon}{Recession} &  \\multicolumn{5}{|c|}{}  \\\\  \n")  )
		nametab <- "recchngmomentsse"
		if(AllEEEU==1){
			rtxt <- "All\ changers"
			nametab <- paste0(nametab,"_all")
		}else if(AllEEEU==2){
			rtxt <- "EE"
			nametab <- paste0(nametab,"_EE")
		}else if(AllEEEU==3){
			if(wc=="wagechangeEUE_wave"| wc=="rawwgchangeEUE_wave"){
				rtxt <- "EUE"
				nametab <- paste0(nametab,"_EUE")
			}else{	
				rtxt <- "EU,UE"
				nametab <- paste0(nametab,"_EUUE")
			}
		}
		tab_wavemomentsse <- xtable(tab_wavemomentsse, digits=2, 
									align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"_",wclab,"_",reclab,"}"))
		
		if(demolbl>=1 & demolbl<=7){
			print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
				  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],"_",wclab,"_",reclab,".tex"))
		}else{
			print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
				  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,"_",wclab,"_",reclab,".tex"))
		}
	}
}# all, EU,UE  & EE

if(wc == "wagechangeEUE_wave"|wc == "reswgchangeEUE_wave"){
	wt = origwt
}
