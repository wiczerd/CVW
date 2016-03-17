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
setwd(wd0)

wagechanges <- readRDS("./Data/balancedwagechanges.RData")
CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
CPSunempRt$unrt <- CPSunempRt$unrt/100

wagechanges <- merge(wagechanges, CPSunempRt, by = "date", all.x = TRUE)

wagechanges[EE==T, balanceweightEUE:=balanceweight]
wagechanges[EU==T, balanceweightEUE:=balanceweight*2]
wagechanges[UE==T, balanceweightEUE:=0.]

toKeep <- c("switchedOcc","switchedInd",
			"Young",
			"HSCol",
			"recIndic",
			"wagechange",
			"wagechange_EUE", 
			"wagechange_all", 
			"balanceweight", 
			"EE","EU","UE",
			"unrt",
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

# Full sample table-------------------------------------------------------------

tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)

tab_fulldist <- array(0., dim=c(9,length(tabqtls)+1))
tab_fulldist[1,1]    <- DTall[ , wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[1,2:tN] <- DTall[ , wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[2,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[2,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[3,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[3,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#expansion
tab_fulldist[4,1]     <- DTall[recIndic == F , wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[4,2:tN]  <- DTall[recIndic == F , wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[5,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[5,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[6,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[6,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#recession
tab_fulldist[7,1]     <- DTall[recIndic == T, wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[7,2:tN]  <- DTall[recIndic == T, wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[8,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[8,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[9,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[9,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]

# EUE wage changes
tab_fulldistEUE <- array(0., dim=c(9,length(tabqtls)+1))
tab_fulldistEUE[1,1]    <- DTall[, wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[1,2:tN] <- DTall[, wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[2,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[2,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[3,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[3,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
#expansion
tab_fulldistEUE[4,1]     <- DTall[recIndic == F , wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[4,2:tN]  <- DTall[recIndic == F , wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[5,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[5,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[6,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[6,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
#recession
tab_fulldistEUE[7,1]     <- DTall[recIndic == T , wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[7,2:tN]  <- DTall[recIndic == T , wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[8,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[8,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[9,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[9,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]


#output it to tables
tab_fulldist <- data.table(tab_fulldist)
names(tab_fulldist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_fulldist) <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
							"All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
							"All\ Workers\ \ ", "Same\ Job,\ \ ","Chng\ Job\ \ ")
tab_fulldistEUE <- data.table(tab_fulldistEUE)
names(tab_fulldistEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_fulldistEUE) <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
							   "All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
							   "All\ Workers\ \ ", "Same\ Job\ \ ","Chng\ Job\ \ ")

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
					"\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
					"\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
tab_fulldist <- xtable(tab_fulldist, label="tab:fulldist", digits=2, 
					align="l|l|lllll", caption="Distribution of earnings changes")
print(tab_fulldist,include.rownames=T, hline.after= c(nrow(tab_fulldist)), 
	  add.to.row=rowtitles, file="./Figures/fulldist.tex")


tab_fulldistEUE <- xtable(tab_fulldistEUE, label="tab:fulldistEUE", digits=2,
					   align="l|l|lllll", caption="Distribution of earnings changes, connecting unemployment spells")
print(tab_fulldistEUE,include.rownames=T, hline.after= c(nrow(tab_fulldistEUE)),
	  add.to.row= rowtitles,file="./Figures/fulldistEUE.tex")


# Full sample var decomp -------------------------------------------------------------

# take out seam effect by regression:
DTall[ , seam:= wave != shift(wave,1,type="lead"), by=id ]
seamcorr<-lm( wagechange_all ~ I(seam & wagechange_all>0) + I(seam & wagechange_all<0) , data=subset(DTall, EE==F&EU==F&UE==F), na.action=na.omit)
seamcorrEUE<-lm( wagechange_allEUE ~ I(seam & wagechange_allEUE>0) + I(seam & wagechange_allEUE<0) , data=subset(DTall, EE==F&EU==F&UE==F), na.action=na.omit)
DTall[ (EE==F&EU==F&UE==F)& wagechange_all<0 & seam==T, wagechange_noseam := wagechange_all + seamcorr$coefficients[3] ]
DTall[ (EE==F&EU==F&UE==F)& wagechange_all>0 & seam==T, wagechange_noseam := wagechange_all + seamcorr$coefficients[2] ]
DTall[ (EE==F&EU==F&UE==F)& wagechange_allEUE<0 & seam==T, wagechange_noseamEUE := wagechange_allEUE + seamcorrEUE$coefficients[3] ]
DTall[ (EE==F&EU==F&UE==F)& wagechange_allEUE>0 & seam==T, wagechange_noseamEUE := wagechange_allEUE + seamcorrEUE$coefficients[2] ]
seamcorr<-lm( wagechange_all ~ I(seam & wagechange_all>0) + I(seam & wagechange_all<0) , data=subset(DTall, EE==T|EU==T|UE==T), na.action=na.omit)
seamcorrEUE<-lm( wagechange_allEUE ~ I(seam & wagechange_allEUE>0) + I(seam & wagechange_allEUE<0) , data=subset(DTall, EE==T|EU==T|UE==T), na.action=na.omit)
DTall[ (EE==T|EU==T|UE==T)& wagechange_all<0 & seam==T, wagechange_noseam := wagechange_all + seamcorr$coefficients[3] ]
DTall[ (EE==T|EU==T|UE==T)& wagechange_all>0 & seam==T, wagechange_noseam := wagechange_all + seamcorr$coefficients[2] ]
DTall[ (EE==T|EU==T|UE==T)& wagechange_allEUE<0 & seam==T, wagechange_noseamEUE := wagechange_allEUE + seamcorrEUE$coefficients[3] ]
DTall[ (EE==T|EU==T|UE==T)& wagechange_allEUE>0 & seam==T, wagechange_noseamEUE := wagechange_allEUE + seamcorrEUE$coefficients[2] ]


totmean <- DTall[, wtd.mean(wagechange_noseam,na.rm=T,weights=allwt) ]
totvar  <- DTall[, sum(allwt*(wagechange_noseam- totmean)^2,na.rm=T) ]

tab_vardec <- array(NA_real_,dim=c(4,4))

tab_vardec[1,1] <- DTall[(EE==F&EU==F&UE==F), sum(allwt*(wagechange_noseam - totmean)^2,na.rm=T) ]/totvar
tab_vardec[1,2] <- DTall[(EE==T&EU==F&UE==F), sum(allwt*(wagechange_noseam - totmean)^2,na.rm=T) ]/totvar
tab_vardec[1,3] <- DTall[(EE==F&(EU==T|UE==T)),sum(allwt*(wagechange_noseam - totmean)^2,na.rm=T) ]/totvar
totwt <- DTall[, sum(allwt,na.rm=T) ]
tab_vardec[2,1] <- DTall[(EE==F&EU==F&UE==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[2,2] <- DTall[(EE==T&EU==F&UE==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[2,3] <- DTall[(EE==F&(EU==T|UE==T)),sum(allwt,na.rm=T) ]/totwt


totmean <- DTall[UE==F, wtd.mean(wagechange_noseamEUE,na.rm=T,weights=allwt) ]
totvar  <- DTall[UE==F, sum(allwt*(wagechange_noseamEUE- totmean)^2,na.rm=T) ]
tab_vardec[3,1] <- DTall[(EE==F&EU==F), sum(allwtEUE*(wagechange_noseamEUE - totmean)^2,na.rm=T) ]/totvar
tab_vardec[3,2] <- DTall[(EE==T&EU==F), sum(allwtEUE*(wagechange_noseamEUE - totmean)^2,na.rm=T) ]/totvar
tab_vardec[3,4] <- DTall[(EE==F&(EU==T)),sum(allwtEUE*(wagechange_noseamEUE - totmean)^2,na.rm=T) ]/totvar
totwt <- DTall[UE==F, sum(allwt,na.rm=T) ]
tab_vardec[4,1] <- DTall[(EE==F&EU==F), sum(allwtEUE,na.rm=T) ]/totwt
tab_vardec[4,2] <- DTall[(EE==T&EU==F), sum(allwtEUE,na.rm=T) ]/totwt
tab_vardec[4,4] <- DTall[(EE==F&(EU==T)),sum(allwtEUE,na.rm=T) ]/totwt

# Full sample quantile-diff decomposition --------------------------
tot51025qtl <- DTall[, wtd.quantile(wagechange_noseam,na.rm=T,weights=allwt, probs=c(0.05,0.1,.25,.75,.9,.95)) ]
tot51025qtlEUE <- DTall[, wtd.quantile(wagechange_noseamEUE,na.rm=T,weights=allwtEUE, probs=c(0.05,0.1,.25,.75,.9,.95)) ]

Nqtls <-length(tot51025qtl)
Ndifs <- Nqtls/2


tab_qtldec <- array(NA_real_,dim=c((Ndifs+1)*2,4))

for( EUEindic in c(F,T)){
	wc <- ifelse(EUEindic,"wagechange_noseamEUE","wagechange_noseam")
	wt <- ifelse(EUEindic,"allwtEUE","allwt")
	qtlshere <- ifelse(rep(EUEindic,Nqtls),tot51025qtlEUE,tot51025qtl)
	ri <- (Ndifs+1)*EUEindic
	ci <- 3+EUEindic
	for(rri in seq(1,Ndifs)){
		tab_qtldec[rri+ri,1] <- DTall[(EE==F&EU==F&UE==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
		tab_qtldec[rri+ri,2] <- DTall[(EE==T&EU==F&UE==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
		tab_qtldec[rri+ri,ci]<- DTall[(EE==F&(EU==T|UE==T)) & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
	}
	tab_qtldec[Ndifs+1+ri,1] <- DTall[(EE==F&EU==F&UE==F)   , sum(eval(as.name(wt)),na.rm=T) ]
	tab_qtldec[Ndifs+1+ri,2] <- DTall[(EE==T&EU==F&UE==F)   , sum(eval(as.name(wt)),na.rm=T) ]
	tab_qtldec[Ndifs+1+ri,ci]<- DTall[(EE==F&(EU==T|UE==T)) , sum(eval(as.name(wt)),na.rm=T) ]
}

rsum <- rowSums(tab_qtldec, na.rm=T)
for(ri in seq(1,nrow(tab_qtldec))){
	tab_qtldec[ri,] <-tab_qtldec[ri,]/rsum[ri]
}

#output it to tables
tab_chngvarqtldec <- data.table(rbind(tab_vardec[1,],tab_qtldec[1:Ndifs,],tab_vardec[3,],tab_qtldec[(2+Ndifs):nrow(tab_qtldec),]) )
names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EU,UE","EUE")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngvarqtldec) <- c("Variance","0.95-0.05","0.9-0.1","0.75-0.25", "Variance\ ","0.95-0.05\ ","0.9-0.1\ ","0.75-0.25\ ","Pct Sample")

tab_chngvarqtldec <- xtable(tab_chngvarqtldec, label="tab:chngvarqtldec", digits=2, 
					   align="l|l|lll", caption="Decomposition of earnings changes")
print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)), file="./Figures/chngvarqtldec.tex")




##########################################################################################

# Only job-changers -----------------------
wagechanges<- wagechanges[ is.finite(switchedOcc), ]
wagechanges[EE==T, balanceweightEUE:=balanceweight]
wagechanges[EU==T, balanceweightEUE:=balanceweight*2]
wagechanges[UE==T, balanceweightEUE:=0.]


tab_chngdist    <- array(0.,dim=c(3,tN))
tab_chngdistEUE <- array(0.,dim=c(3,tN))
tab_chngdistEUE_EEEUE <- array(0.,dim=c(6,tN))
tab_chngdist[1,1]    <- wagechanges[(EE|EU|UE),wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
tab_chngdist[1,2:tN] <- wagechanges[ EE|EU|UE,wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs= tabqtls)]
tab_chngdist[2,1]    <- wagechanges[(EE|EU|UE)& switchedOcc==F,wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
tab_chngdist[2,2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==F,wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]
tab_chngdist[3,1]    <- wagechanges[(EE|EU|UE)& switchedOcc==T,wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
tab_chngdist[3,2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==T,wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]

tab_chngdistEUE[1,1]    <- wagechanges[(EE|EU|UE),wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
tab_chngdistEUE[1,2:tN] <- wagechanges[ EE|EU|UE,wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs= tabqtls)]
tab_chngdistEUE[2,1]    <- wagechanges[(EE|EU|UE)& switchedOcc==F,wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
tab_chngdistEUE[2,2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==F,wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
tab_chngdistEUE[3,1]    <- wagechanges[(EE|EU|UE)& switchedOcc==T,wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
tab_chngdistEUE[3,2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==T,wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]

for(EUhere in c(F,T)){
	EEhere = ifelse(EUhere==F, T,F)
	
	tab_chngdistEUE_EEEUE[(1+EUhere*3),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)),wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_EEEUE[(1+EUhere*3),2:tN] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)),wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs= tabqtls)]
	tab_chngdistEUE_EEEUE[(2+EUhere*3),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==F,wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_EEEUE[(2+EUhere*3),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==F,wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
	tab_chngdistEUE_EEEUE[(3+EUhere*3),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==T,wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_EEEUE[(3+EUhere*3),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==T,wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
}

#output it to tables
tab_chngdist <- data.table(tab_chngdist)
names(tab_chngdist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngdist) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")
tab_chngdistEUE <- data.table(tab_chngdistEUE)
names(tab_chngdistEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_chngdistEUE) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")

tab_chngdist <- xtable(tab_chngdist, label="tab:fulldist", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job changers")
print(tab_chngdist,include.rownames=T, hline.after= c(0,nrow(tab_chngdist)), file="./Figures/chngdist.tex")

tab_chngdistEUE <- xtable(tab_chngdistEUE, label="tab:chngdistEUE", digits=2, 
						  align="l|l|lllll", caption="Distribution of earnings changes among job changers, connecting unemployment spells")
print(tab_chngdistEUE,include.rownames=T, hline.after= c(0,nrow(tab_chngdistEUE)), file="./Figures/chngdistEUE.tex")


tab_chngdist_rec    <- array(0.,dim=c(6,tN))
tab_chngdistEUE_rec <- array(0.,dim=c(6,tN))

tab_chngdist_rec_EEEUE    <- array(0.,dim=c(12,tN))
tab_chngdistEUE_rec_EEEUE <- array(0.,dim=c(12,tN))

ri = 0
for(recHere in c(F,T)){
	tab_chngdist_rec[(1+ri*3),1]    <- wagechanges[(EE|EU|UE) & recIndic == recHere,
											   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
	tab_chngdist_rec[(1+ri*3),2:tN] <- wagechanges[ (EE|EU|UE) & recIndic == recHere,
												wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs= tabqtls)]
	tab_chngdist_rec[(2+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
											   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
	tab_chngdist_rec[(2+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
												  wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]
	tab_chngdist_rec[(3+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
											   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
	tab_chngdist_rec[(3+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
												  wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]
	
	tab_chngdistEUE_rec[(1+ri*3),1]    <- wagechanges[(EE|EU|UE) & recIndic == recHere,
												   wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_rec[(1+ri*3),2:tN] <- wagechanges[ (EE|EU|UE) & recIndic == recHere,
													wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs= tabqtls)]
	tab_chngdistEUE_rec[(2+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
												   wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_rec[(2+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
													  wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
	tab_chngdistEUE_rec[(3+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
												   wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_rec[(3+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
													  wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
	for(EUhere in c(F,T)){
		EEhere = ifelse(EUhere==F, T,F)
	
		tab_chngdist_rec_EEEUE[(1+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		recIndic == recHere,
													   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
		tab_chngdist_rec_EEEUE[(1+ri*3 + 6*EUhere),2:tN] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) 
																		& recIndic == recHere,
														wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs= tabqtls)]
		tab_chngdist_rec_EEEUE[(2+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		switchedOcc==F & recIndic == recHere,
													   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
		tab_chngdist_rec_EEEUE[(2+ri*3 + 6*EUhere),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   switchedOcc==F & recIndic == recHere,
														  wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]
		tab_chngdist_rec_EEEUE[(3+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		switchedOcc==T & recIndic == recHere,
													   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
		tab_chngdist_rec_EEEUE[(3+ri*3 + 6*EUhere),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   switchedOcc==T & recIndic == recHere,
														  wtd.quantile(wagechange,na.rm=T,weights=balanceweight, probs=tabqtls)]
		
		tab_chngdistEUE_rec_EEEUE[(1+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   recIndic == recHere,
														  wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
		tab_chngdistEUE_rec_EEEUE[(1+ri*3 + 6*EUhere),2:tN] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   recIndic == recHere,
														   wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs= tabqtls)]
		tab_chngdistEUE_rec_EEEUE[(2+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   switchedOcc==F & recIndic == recHere,
														  wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
		tab_chngdistEUE_rec_EEEUE[(2+ri*3 + 6*EUhere),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																			  switchedOcc==F & recIndic == recHere,
															 wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
		tab_chngdistEUE_rec_EEEUE[(3+ri*3 + 6*EUhere),1]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																		   switchedOcc==T & recIndic == recHere,
														  wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
		tab_chngdistEUE_rec_EEEUE[(3+ri*3 + 6*EUhere),2:tN]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)) &
																			  switchedOcc==T & recIndic == recHere,
															 wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
		
	}
	ri = ri+1
}

#output it to tables
tab_chngdist_rec <- data.table(tab_chngdist_rec)
names(tab_chngdist_rec) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_chngdist_rec) <- c("Chng\ Job, All -- Exp", "Chng\ Job, Same\ Occ -- Exp", "Chng\ Job & Switch\ Occ -- Exp",
								"Chng\ Job, All -- Rec", "Chng\ Job, Same\ Occ -- Rec", "Chng\ Job & Switch\ Occ -- Rec")
tab_chngdistEUE_rec <- data.table(tab_chngdistEUE_rec)
names(tab_chngdistEUE_rec) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_chngdistEUE_rec) <- c("Chng\ Job, All -- Exp", "Chng\ Job, Same\ Occ -- Exp", "Chng\ Job & Switch\ Occ -- Exp",
								   "Chng\ Job, All -- Rec", "Chng\ Job, Same\ Occ -- Rec", "Chng\ Job & Switch\ Occ -- Rec")

tab_chngdist_rec <- xtable(tab_chngdist_rec, label="tab:chngdist_rec", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job changers in recession and expansion")
print(tab_chngdist_rec,include.rownames=T, 
	  hline.after= c(0,nrow(tab_chngdist_rec)/2,nrow(tab_chngdist_rec)), file="./Figures/chngdist_rec.tex")

tab_chngdistEUE_rec <- xtable(tab_chngdistEUE_rec, label="tab:chngdistEUE", digits=2, 
						  align="l|l|lllll", caption="Distribution of earnings changes among job changers in recession and expansion, connecting unemployment spells")
print(tab_chngdistEUE_rec,include.rownames=T, 
	  hline.after= c(0,nrow(tab_chngdistEUE_rec)/2,nrow(tab_chngdistEUE_rec)), file="./Figures/chngdistEUE_rec.tex")

tab_chngdistEUE_recEE <- data.table(rbind(tab_chngdistEUE_EEEUE[1:3,],tab_chngdistEUE_rec_EEEUE[1:6,]))
tab_chngdistEUE_recEUE <- data.table(rbind(tab_chngdistEUE_EEEUE[4:6,],tab_chngdistEUE_rec_EEEUE[7:12,]))
names(tab_chngdistEUE_recEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
names(tab_chngdistEUE_recEE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rnvec <- c("All job changers", "Occ stayers", "Occ movers",
		   "All job changers\ ", "Occ stayers\ ", "Occ movers\ ",
		   "All job changers\ \ ", "Occ stayers\ \ ", "Occ movers\ \ ")
rownames(tab_chngdistEUE_recEE) <-rnvec
rownames(tab_chngdistEUE_recEUE) <-rnvec
rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
											  "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
											  "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
tab_chngdistEUE_recEE <- xtable(tab_chngdistEUE_recEE, label="tab:chngdistEUE_recEE", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job-job changers")
print(tab_chngdistEUE_recEE,include.rownames=T, hline.after= c(nrow(tab_chngdistEUE_recEE)), 
	  add.to.row=rowtitles, file="./Figures/chngdistEUE_recEE.tex")
tab_chngdistEUE_recEUE <- xtable(tab_chngdistEUE_recEUE, label="tab:chngdistEUE_recEE", digits=2, 
								align="l|l|lllll", caption="Distribution of earnings changes among those transitioning through unemployment")
print(tab_chngdistEUE_recEUE,include.rownames=T, hline.after= c(nrow(tab_chngdistEUE_recEUE)), 
	  add.to.row=rowtitles, file="./Figures/chngdistEUE_recEUE.tex")


# Quantile Regressions ---------------

#qr_EUEunrt <- rq( wagechange_EUE ~  factor(EU) + factor(switchedOcc) + unrt,tau= tabqtls, data=wagechanges, weights = balanceweight)
lm_EUEunrt <- lm( wagechange_EUE ~  factor(EU) + factor(switchedOcc) + unrt, data=wagechanges, weights = balanceweightEUE)
lmqr_EUEunrt <- list(lm_EUEunrt)
for(ti in seq(1,length(tabqtls))){
	lmqr_EUEunrt[[ti+1]] <- rq( wagechange_EUE ~  factor(EU) + factor(switchedOcc) + unrt,tau= tabqtls[ti], data=wagechanges, weights = balanceweightEUE)
}

texreg(lmqr_EUEunrt,custom.model.names=c("OLS","0.1","0.25","0.5","0.75","0.9"), reorder.coef=c(3,2,4,1) ,
	   custom.coef.names=c("Const","Unemp Indic","Switched Occ", "Unemp Rt"), file="./Figures/lmqr_EUEunrt.tex")

qr_swEUE_rec <- rq( wagechange_EUE ~  factor(switchedOcc)+ factor(EU),tau= tabqtls, data=subset(wagechanges,recIndic==T), weights = balanceweightEUE)
qr_swEUE_exp <- rq( wagechange_EUE ~  factor(switchedOcc)+ factor(EU),tau= tabqtls, data=subset(wagechanges,recIndic==F), weights = balanceweightEUE)
lm_swEUE_rec <- lm( wagechange_EUE ~  factor(switchedOcc)+ factor(EU),data=subset(wagechanges,recIndic==T), weights = balanceweightEUE)
lm_swEUE_exp <- lm( wagechange_EUE ~  factor(switchedOcc)+ factor(EU),data=subset(wagechanges,recIndic==F), weights = balanceweightEUE)

tab_lmqr_swEUE_exprec <- cbind(t(qr_swEUE_exp$coefficients), t(qr_swEUE_rec$coefficients) )
tab_lmqr_swEUE_exprec <- rbind(cbind(t(lm_swEUE_exp$coefficients),t(lm_swEUE_rec$coefficients)),
							   tab_lmqr_swEUE_exprec)
rownames(tab_lmqr_swEUE_exprec) <- c("OLS","0.1","0.25","0.5","0.75","0.9")
colnames(tab_lmqr_swEUE_exprec) <- c("Const","Switched\ Occ", "Unemp\ Indic","Const","Switched\ Occ", "Unemp\ Indic")
tab_lmqr_swEUE_exprec <- tab_lmqr_swEUE_exprec[, c(2,3,1,5,6,4)]
tab_lmqr_swEUE_exprec <- xtable(tab_lmqr_swEUE_exprec, label="tab:lmqr_swEUE_exprec", digits=2, align="l|lll|lll",
								caption="Earnings change regressions in expansions and recessions")
addcyc <-list(pos = list(-1),command="\\hline\\hline & \\multicolumn{3}{c}{\\color{Maroon}{Expansion} } &
			  \\multicolumn{3}{c}{\\color{Maroon}{Recession} \\\\ \n")
print(tab_lmqr_swEUE_exprec,include.rownames=T, add.to.row=addcyc,
	  hline.after= c(nrow(tab_lmqr_swEUE_exprec)), file="./Figures/lmqr_swEUE_exprec.tex")

qtlregsco <- data.table(cbind( rbind(array(0,dim=c(length(tabqtls),1)), array(1,dim=c(length(tabqtls),1))),
						rbind(t(qr_swEUE_exp$coefficients), t(qr_swEUE_rec$coefficients) )) )
qtlregsco <- cbind(rep(tabqtls,2),qtlregsco)
names(qtlregsco) <- c("Quantile","recIndic","Const","Switched","EUE")
qtlregsco[ , recIndic := as.factor(recIndic)]
ggplot(qtlregsco, aes(x=Quantile,y=Switched, color= recIndic)) +
	geom_line() +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Return to occupation switch") +
	xlab("Quantile") 
ggsave(filename = "Figures/lmqr_swEUE_sw.png",height= 5,width=10)
ggsave(filename = "Figures/lmqr_swEUE_sw.eps",height= 5,width=10)

ggplot(qtlregsco, aes(x=Quantile,y=EUE, color= recIndic)) +
	geom_line() +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Cost of unemployment") +
	xlab("Quantile") 
ggsave(filename = "Figures/lmqr_swEUE_EUE.png",height= 5,width=10)
ggsave(filename = "Figures/lmqr_swEUE_EUE.eps",height= 5,width=10)


#lmqr_OccIndunrt0 <-rq( wagechange ~  EU + UE + switchedOcc+ switchedInd  + unrt,tau= tabqtls, data=wagechanges, weights = balanceweight)
#lmqr_Occunrt0 <-rq( wagechange ~  EU + UE + switchedOcc  + unrt,tau= tabqtls, data=wagechanges, weights = balanceweight)
#lmqr_Indunrt0 <-rq( wagechange ~  EU + UE + switchedInd  + unrt,tau= tabqtls, data=wagechanges, weights = balanceweight)

lm_EUEOccIndunrt <- lm( wagechange_EUE ~  EU + switchedOcc + switchedInd + unrt, data=wagechanges, weights = balanceweightEUE)
lmqr_EUEOccIndnrt0 <-rq( wagechange_EUE ~  EU + switchedOcc+ switchedInd  + unrt,tau= tabqtls, data=wagechanges, weights = balanceweightEUE)
lmqr_EUEOccIndnrt <- list(lm_EUEOccIndunrt)
for(ti in seq(1,length(tabqtls))){
	lmqr_EUEOccIndnrt[[ti+1]] <- rq( wagechange_EUE ~  factor(EU) + factor(switchedOcc)+ factor(switchedInd)  + unrt,tau= tabqtls[ti], data=wagechanges, weights = balanceweightEUE)
}