# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)
library(xtable)
library(quantreg)


wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

wagechanges <- readRDS("./Data/balancedwagechanges.RData")
CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
wagechanges <- merge(wagechanges, CPSunempRt, by = "date", all.x = TRUE)


toKeep <- c("switchedOcc",
			"Young",
			"HSCol",
			"recIndic",
			"wagechange",
			"wagechange_EUE", 
			"wagechange_all", 
			"balanceweight", 
			"EE","EU","UE",
			"unrt")

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

tab_fulldist <- array(0., dim=c(6,length(tabqtls)+1))
tab_fulldist[1,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[1,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[2,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[2,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#expansion
tab_fulldist[3,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[3,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[4,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[4,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#recession
tab_fulldist[5,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[5,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[6,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[6,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]

# EUE wage changes
tab_fulldistEUE <- tab_fulldist
tab_fulldistEUE[1,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[1,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[2,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[2,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
#expansion
tab_fulldistEUE[3,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[3,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[4,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[4,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
#recession
tab_fulldistEUE[5,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[5,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]
tab_fulldistEUE[6,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwtEUE)]
tab_fulldistEUE[6,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwtEUE, probs=tabqtls)]


#output it to tables
tab_fulldist <- data.table(tab_fulldist)
names(tab_fulldist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_fulldist) <- c("Same\ Job","Chng\ Job","Same\ Job,\ Exp","Chng\ Job,\ Exp","Same\ Job,\ Rec","Chng\ Job,\ Rec")
tab_fulldistEUE <- data.table(tab_fulldistEUE)
names(tab_fulldistEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_fulldistEUE) <- c("Same\ Job","Chng\ Job","Same\ Job,\ Exp","Chng\ Job,\ Exp","Same\ Job,\ Rec","Chng\ Job,\ Rec")

tab_fulldist <- xtable(tab_fulldist, label="tab:fulldist", digits=2, 
					align="l|l|lllll", caption="Distribution of earnings changes")
print(tab_fulldist,include.rownames=T, hline.after= c(0,nrow(tab_fulldist)), file="./Figures/fulldist.tex")

tab_fulldistEUE <- xtable(tab_fulldistEUE, label="tab:fulldistEUE", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes, connecting unemployment spells")
print(tab_fulldistEUE,include.rownames=T, hline.after= c(0,nrow(tab_fulldist)), file="./Figures/fulldistEUE.tex")


# Full sample var decomp -------------------------------------------------------------

totmean <- DTall[, wtd.mean(wagechange_all,na.rm=T,weights=allwt) ]
totvar  <- DTall[, sum(allwt*(wagechange_all- totmean)^2,na.rm=T) ]

tab_vardec <- array(0.,dim=c(4,4))

tab_vardec[1,1] <- DTall[(EE==F&EU==F&UE==F), sum(allwt*(wagechange_all - totmean)^2,na.rm=T) ]/totvar
tab_vardec[1,2] <- DTall[(EE==T&EU==F&UE==F), sum(allwt*(wagechange_all - totmean)^2,na.rm=T) ]/totvar
tab_vardec[1,3] <- DTall[(EE==F&(EU==T|UE==T)),sum(allwt*(wagechange_all - totmean)^2,na.rm=T) ]/totvar
totwt <- DTall[, sum(allwt,na.rm=T) ]
tab_vardec[2,1] <- DTall[(EE==F&EU==F&UE==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[2,2] <- DTall[(EE==T&EU==F&UE==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[2,3] <- DTall[(EE==F&(EU==T|UE==T)),sum(allwt,na.rm=T) ]/totwt


totmean <- DTall[UE==F, wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt) ]
totvar  <- DTall[UE==F, sum(allwt*(wagechange_allEUE- totmean)^2,na.rm=T) ]
tab_vardec[3,1] <- DTall[(EE==F&EU==F), sum(allwt*(wagechange_allEUE - totmean)^2,na.rm=T) ]/totvar
tab_vardec[3,2] <- DTall[(EE==T&EU==F), sum(allwt*(wagechange_allEUE - totmean)^2,na.rm=T) ]/totvar
tab_vardec[3,4] <- DTall[(EE==F&(EU==T)),sum(allwt*(wagechange_allEUE - totmean)^2,na.rm=T) ]/totvar
totwt <- DTall[UE==F, sum(allwt,na.rm=T) ]
tab_vardec[4,1] <- DTall[(EE==F&EU==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[4,2] <- DTall[(EE==T&EU==F), sum(allwt,na.rm=T) ]/totwt
tab_vardec[4,4] <- DTall[(EE==F&(EU==T)),sum(allwt,na.rm=T) ]/totwt

# Full sample quantile-diff decomposition --------------------------
tot1025qtl <- DTall[, wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=c(.1,.25,.75,.9)) ]
tot1025qtlEUE <- DTall[, wtd.quantile(wagechange_all,na.rm=T,weights=allwtEUE, probs=c(.1,.25,.75,.9)) ]



##########################################################################################

# Only job-changers -----------------------
wagechanges<- wagechanges[ is.finite(switchedOcc), ]

tab_chngdist    <- array(0.,dim=c(3,tN))
tab_chngdistEUE <- array(0.,dim=c(3,tN))
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
ri = 0
for(recHere in c(F,T)){
	tab_chngdist_rec[(1+ri*3),1]    <- wagechanges[(EE|EU|UE) & recIndic == recHere,
											   wtd.mean(wagechange,na.rm=T,weights=balanceweight)]
	tab_chngdist_rec[(1+ri*3),2:tN] <- wagechanges[ (EE|EU|UE & recIndic == recHere) & recIndic == recHere,
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
	tab_chngdistEUE_rec[(1+ri*3),2:tN] <- wagechanges[ (EE|EU|UE & recIndic == recHere) & recIndic == recHere,
													wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs= tabqtls)]
	tab_chngdistEUE_rec[(2+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
												   wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_rec[(2+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==F & recIndic == recHere,
													  wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
	tab_chngdistEUE_rec[(3+ri*3),1]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
												   wtd.mean(wagechange_EUE,na.rm=T,weights=balanceweight)]
	tab_chngdistEUE_rec[(3+ri*3),2:tN]    <- wagechanges[(EE|EU|UE)& switchedOcc==T & recIndic == recHere,
													  wtd.quantile(wagechange_EUE,na.rm=T,weights=balanceweight, probs=tabqtls)]
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


# Quantile Regressions ---------------

	rhere <- rq( wagechange_EUE ~  factor(EU) + factor(switchedOcc) + unrt,tau= tabqtls, data=wagechanges, weights = balanceweight)
	betaptsR[,ri] = rhere$coefficients



