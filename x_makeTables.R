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

CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
CPSunempRt$unrt <- CPSunempRt$unrt/100


##########################################################################################
# By wave -----------------------
toKeep_wave <- c("switchedOcc_wave","switchedInd",
            "Young",
            "HSCol",
            "recIndic","recIndic_wave",
            "wagechange",
            "wagechange_wave", 
            "wagechange_seam",
            "EE_wave","EU_wave","UE_wave",
            "unrt","wpfinwgt","waveweight",
            "wave","id")
DTseam <- readRDS("./Data/DTseam.RData")
DTseam <- merge(DTseam, CPSunempRt, by = "date", all.x = TRUE)

# select toKeep columns only
DTseam <- DTseam[, toKeep_wave, with = FALSE]
DTseam <- subset(DTseam, is.finite(wpfinwgt) & is.finite(wagechange_wave))

DTseam<-DTseam[ is.finite(EE_wave)&is.finite(EU_wave)&is.finite(UE_wave),]
tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)

tab_wavedist <- array(0., dim=c(9,length(tabqtls)+1))

tab_wavedist[1,1]    <- DTseam[                                   ,     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
tab_wavedist[1,2:tN] <- DTseam[                                   , wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
tab_wavedist[2,1]    <- DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
tab_wavedist[2,2:tN] <- DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
tab_wavedist[3,1]    <- DTseam[  EU_wave==T|UE_wave==T|EE_wave==T ,     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
tab_wavedist[3,2:tN] <- DTseam[  EU_wave==T|UE_wave==T|EE_wave==T , wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
#expansion/recession
for(rI in c(T,F)){
  rix = rI*3+3
  tab_wavedist[1+rix,1]   <- DTseam[recIndic_wave == rI,     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
  tab_wavedist[1+rix,2:tN]<- DTseam[recIndic_wave == rI, wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
  tab_wavedist[2+rix,1]   <- DTseam[recIndic_wave == rI & !(EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
  tab_wavedist[2+rix,2:tN]<- DTseam[recIndic_wave == rI & !(EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
  tab_wavedist[3+rix,1]   <- DTseam[recIndic_wave == rI &  (EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=waveweight)]
  tab_wavedist[3+rix,2:tN]<- DTseam[recIndic_wave == rI &  (EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=tabqtls)]
}



#output it to tables
dat_wavedist <- tab_wavedist
tab_wavedist <- data.table(tab_wavedist)
names(tab_wavedist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rnames <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
            "All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
            "All\ Workers\ \ ", "Same\ Job,\ \ ","Chng\ Job\ \ ")
rownames(tab_wavedist) <- rnames

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
                                              "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
                                              "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
tab_wavedist <- xtable(tab_wavedist, digits=2, 
                       align="l|l|lllll", caption="Distribution of earnings changes \\label{tab:wavedist}")
print(tab_wavedist,include.rownames=T, hline.after= c(nrow(tab_wavedist)), 
      add.to.row=rowtitles, file="./Figures/wavedist.tex")

tab_wavevardec <- array(NA_real_,dim=c(2,3))

totmean <- DTseam[, wtd.mean(wagechange_wave,na.rm=T,weights=waveweight) ]
totvar  <- DTseam[, sum(waveweight*(wagechange_wave- totmean)^2,na.rm=T) ]
tab_wavevardec[1,1] <- DTseam[(EE_wave==F& EU_wave==F&UE_wave==F) , sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,2] <- DTseam[(EE_wave==T& EU_wave==F&UE_wave==F) , sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T)), sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
totwt <- DTseam[, sum(waveweight,na.rm=T) ]
tab_wavevardec[2,1] <- DTseam[(EE_wave==F& EU_wave==F&UE_wave==F) , sum(waveweight,na.rm=T) ]/totwt
tab_wavevardec[2,2] <- DTseam[(EE_wave==T& EU_wave==F&UE_wave==F) , sum(waveweight,na.rm=T) ]/totwt
tab_wavevardec[2,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T)), sum(waveweight,na.rm=T) ]/totwt


# Full sample quantile-diff decomposition --------------------------
tot51025qtl <- DTseam[, wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=c(0.05,0.1,.25,.75,.9,.95)) ]

Nqtls <-length(tot51025qtl)
Ndifs <- Nqtls/2

tab_waveqtldec <- array(NA_real_,dim=c(Ndifs+1,3))

qtlshere <- tot51025qtl

for(rri in seq(1,Ndifs)){
  tab_waveqtldec[rri,1] <- DTseam[(EE_wave==F& EU_wave==F&UE_wave==F)  & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
  tab_waveqtldec[rri,2] <- DTseam[(EE_wave==T& EU_wave==F&UE_wave==F)  & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
  tab_waveqtldec[rri,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T)) & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
}
tab_waveqtldec[Ndifs+1,1] <- DTseam[(EE_wave==F& EU_wave==F&UE_wave==F)  , sum(waveweight,na.rm=T) ]
tab_waveqtldec[Ndifs+1,2] <- DTseam[(EE_wave==T& EU_wave==F&UE_wave==F)  , sum(waveweight,na.rm=T) ]
tab_waveqtldec[Ndifs+1,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T)) , sum(waveweight,na.rm=T) ]


rsum <- rowSums(tab_waveqtldec, na.rm=T)
for(ri in seq(1,nrow(tab_waveqtldec))){
  tab_waveqtldec[ri,] <-tab_waveqtldec[ri,]/rsum[ri]
}

#output it to tables
tab_chngvarqtldec <- data.table(rbind(tab_wavevardec[1,],tab_waveqtldec[1:Ndifs,]) )
names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EU,UE")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngvarqtldec) <- c("Variance","0.95-0.05","0.9-0.1","0.75-0.25", "Variance\ ","0.95-0.05\ ","0.9-0.1\ ","0.75-0.25\ ","Pct Sample")

tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
                            align="l|l|ll", caption="Decomposition of earnings change dispersion \\label{tab:wavechngvarqtldec}")
print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file="./Figures/wavechngvarqtldec.tex")


# recession and expansion
for (rI in c(F,T)){
   
  totmean <- DTseam[recIndic_wave==rI, wtd.mean(wagechange_wave,na.rm=T,weights=waveweight) ]
  totvar  <- DTseam[recIndic_wave==rI, sum(waveweight*(wagechange_wave- totmean)^2,na.rm=T) ]
  tab_wavevardec[1,1] <- DTseam[(EE_wave==F&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
  tab_wavevardec[1,2] <- DTseam[(EE_wave==T&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
  tab_wavevardec[1,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T))&recIndic_wave==rI, sum(waveweight*(wagechange_wave - totmean)^2,na.rm=T) ]/totvar
  totwt <- DTseam[recIndic_wave==rI, sum(waveweight,na.rm=T) ]
  tab_wavevardec[2,1] <- DTseam[(EE_wave==F&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(waveweight,na.rm=T) ]/totwt
  tab_wavevardec[2,2] <- DTseam[(EE_wave==T&EU_wave==F&UE_wave==F) & recIndic_wave==rI, sum(waveweight,na.rm=T) ]/totwt
  tab_wavevardec[2,3] <- DTseam[(EE_wave==F&(EU_wave==T|UE_wave==T))&recIndic_wave==rI, sum(waveweight,na.rm=T) ]/totwt

  # Seam sample quantile-diff decomposition --------------------------
  tot51025qtl <- DTseam[ recIndic_wave==rI, wtd.quantile(wagechange_wave,na.rm=T,weights=waveweight, probs=c(0.05,0.1,.25,.75,.9,.95)) ]
  
  Nqtls <-length(tot51025qtl)
  Ndifs <- Nqtls/2
  
  tab_waveqtldec <- array(NA_real_,dim=c(Ndifs+1,3))
  
  qtlshere <- tot51025qtl
  for(rri in seq(1,Ndifs)){
    tab_waveqtldec[rri,1] <- DTseam[recIndic_wave==rI & (EE_wave==F&EU_wave==F&UE_wave==F)   & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
    tab_waveqtldec[rri,2] <- DTseam[recIndic_wave==rI & (EE_wave==T&EU_wave==F&UE_wave==F)   & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
    tab_waveqtldec[rri,3] <- DTseam[recIndic_wave==rI & (EE_wave==F&(EU_wave==T|UE_wave==T)) & (wagechange_wave > qtlshere[Nqtls-rri+1] | wagechange_wave < qtlshere[rri]), sum(waveweight ,na.rm=T)]
  }
  tab_waveqtldec[Ndifs+1,1] <- DTseam[recIndic_wave==rI & (EE_wave==F&EU_wave==F&UE_wave==F)  , sum(waveweight,na.rm=T) ]
  tab_waveqtldec[Ndifs+1,2] <- DTseam[recIndic_wave==rI & (EE_wave==T&EU_wave==F&UE_wave==F)  , sum(waveweight,na.rm=T) ]
  tab_waveqtldec[Ndifs+1,3] <- DTseam[recIndic_wave==rI & (EE_wave==F&(EU_wave==T|UE_wave==T)), sum(waveweight,na.rm=T) ]

  
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
  print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0("./Figures/wavechngvarqtldec_rec",rI,".tex"))
  
}


# Quantiles by job change  -----------------------
tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)

tab_wavechngdist    <- array(0.,dim=c(9,tN))


wc <- "wagechange_wave"
wt <- "waveweight"

tab_wavechngdist[1,1]    <- DTseam[(EE_wave|EU_wave|UE_wave),                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[1,2:tN] <- DTseam[ EE_wave|EU_wave|UE_wave ,                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
tab_wavechngdist[2,1]    <- DTseam[(EE_wave|EU_wave|UE_wave)& switchedOcc_wave==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[2,2:tN] <- DTseam[(EE_wave|EU_wave|UE_wave)& switchedOcc_wave==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
tab_wavechngdist[3,1]    <- DTseam[(EE_wave|EU_wave|UE_wave)& switchedOcc_wave==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[3,2:tN] <- DTseam[(EE_wave|EU_wave|UE_wave)& switchedOcc_wave==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
for(EUhere in c(F,T)){
	EEhere = ifelse(EUhere==F, T,F)
	eidx = as.integer(EUhere)*3
	tab_wavechngdist[4+eidx,1]    <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere)),                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[4+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere)),                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_wavechngdist[5+eidx,1]    <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[5+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavechngdist[6+eidx,1]    <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist[6+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))& switchedOcc_wave==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
}

dat_wavechngdist <- tab_wavechngdist

#output it to tables
tab_wavechngdist_out <- data.table(tab_wavechngdist[1:3, ])
names(tab_wavechngdist_out) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_wavechngdist_out) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")

tab_wavechngdist_out <- xtable(tab_wavechngdist_out, label="tab:wavechngdist", digits=2, 
                           align="l|l|lllll", caption="Distribution of earnings changes among job changers")
print(tab_wavechngdist_out,include.rownames=T, hline.after= c(0,nrow(tab_wavechngdist_out)), file="./Figures/wavechngdist.tex")

# rec and expansion
tab_wavechngdist_rec    <- array(NA,dim=c(6,tN))
tab_wavechngdist_rec_EEEU    <- array(NA,dim=c(12,tN))
for(recHere in c(F,T)){
	ridx = as.integer(recHere)*3
	wc <- "wagechange_wave"
	wt <- "waveweight"
	
	tab_wavechngdist_rec[1+ridx,1   ] <- DTseam[(EE_wave|EU_wave|UE_wave)                & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[1+ridx,2:tN] <- DTseam[(EE_wave|EU_wave|UE_wave)                & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_wavechngdist_rec[2+ridx,1   ] <- DTseam[(EE_wave|EU_wave|UE_wave)&switchedOcc_wave==F & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[2+ridx,2:tN] <- DTseam[(EE_wave|EU_wave|UE_wave)&switchedOcc_wave==F & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavechngdist_rec[3+ridx,1   ] <- DTseam[(EE_wave|EU_wave|UE_wave)&switchedOcc_wave==T & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[3+ridx,2:tN] <- DTseam[(EE_wave|EU_wave|UE_wave)&switchedOcc_wave==T & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	
	for(EUhere in c(F,T)){
		EEhere = ifelse(EUhere==F, T,F)
		eidx = 6*as.integer(EUhere)
		tab_wavechngdist_rec_EEEU[1+ridx+eidx,1   ] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))                &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec_EEEU[1+ridx+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))                &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
		tab_wavechngdist_rec_EEEU[2+ridx+eidx,1   ] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))&switchedOcc_wave==F &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec_EEEU[2+ridx+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))&switchedOcc_wave==F &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_wavechngdist_rec_EEEU[3+ridx+eidx,1   ] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))&switchedOcc_wave==T &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavechngdist_rec_EEEU[3+ridx+eidx,2:tN] <- DTseam[(EE_wave==EEhere&(EU_wave==EUhere|UE_wave==EUhere))&switchedOcc_wave==T &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}
}
#output it to tables
dat_wavechngdist_rec <-tab_wavechngdist_rec
dat_wavechngdist_rec_EEEU <-tab_wavechngdist_rec_EEEU
cnames <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rnames <- c("All job changers    ", "Occ stayers    ","Occ movers    ",
            "All job changers\   ", "Occ stayers\   ","Occ movers\   ",
            "All job changers\ \ ", "Occ stayers\ \ ","Occ movers\ \ ")

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
                                              "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
                                              "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )


tab_wavechngdist_rec <- data.table(rbind(dat_wavechngdist[1:3, ],dat_wavechngdist_rec[1:6, ]))
names(tab_wavechngdist_rec) <- cnames
rownames(tab_wavechngdist_rec) <- rnames
tab_wavechngdist_rec <- xtable(tab_wavechngdist_rec, digits=2, 
                           align="l|l|lllll", caption="Distribution of earnings changes among job changers in recession and expansion \\label{tab:wavechngdist_rec}")
print(tab_wavechngdist_rec,include.rownames=T, hline.after= c(nrow(tab_wavechngdist_rec)), 
      add.to.row=rowtitles, file="./Figures/wavechngdist_rec.tex")



tab_wavechngdist_recEE     <- data.table(rbind(dat_wavechngdist[4:6, ],dat_wavechngdist_rec_EEEU[1:6 , ]))
tab_wavechngdist_recEUUE   <- data.table(rbind(dat_wavechngdist[7:9, ],dat_wavechngdist_rec_EEEU[7:12, ]))

names(tab_wavechngdist_recEUUE) <- cnames
names(tab_wavechngdist_recEE) <- cnames

rownames(tab_wavechngdist_recEE) <-rnames
rownames(tab_wavechngdist_recEUUE) <-rnames


tab_wavechngdist_recEE <- xtable(tab_wavechngdist_recEE, label="tab:wavechngdist_recEE", digits=2, 
                             align="l|l|lllll", caption="Distribution of earnings changes among job-job changers")
print(tab_wavechngdist_recEE,include.rownames=T, hline.after= c(nrow(tab_wavechngdist_recEE)), 
      add.to.row=rowtitles, file="./Figures/wavechngdist_recEE.tex")

tab_wavechngdist_recEUUE <- xtable(tab_wavechngdist_recEUUE, label="tab:wavechngdistEUE_recEE", digits=2, 
                               align="l|l|lllll", caption="Distribution of earnings changes among those transitioning through unemployment")
print(tab_wavechngdist_recEUUE,include.rownames=T, hline.after= c(nrow(tab_wavechngdist_recEUUE)), 
      add.to.row=rowtitles, file="./Figures/wavechngdist_recEUUE.tex")


# Are these the same distributions ? -----------------------------------



# first do KS test with each subsample 
DTseam[ , s1:= T]
DTseam[ , s2:= ifelse(switchedOcc==F, T,F) ]
DTseam[ , s3:= ifelse(switchedOcc==T, T,F)]
#expansion
DTseam[ , s4:= ifelse(recIndic == F,T,F)]
DTseam[ , s5:= ifelse(recIndic == F & switchedOcc==F,T,F)]
DTseam[ , s6:= ifelse(recIndic == F & switchedOcc==T,T,F)]
#recession
DTseam[ , s7 := ifelse(recIndic == T, T,F)]
DTseam[ , s8 := ifelse(recIndic == T & switchedOcc==F,T,F)]
DTseam[ , s9 := ifelse(recIndic == T & switchedOcc==T,T,F)]

NS = 9

tab_chngdist_ks <- array(NA_real_,dim=c(NS-1,NS-1,2))
for(EUEindic in c(F,T)){
	idx = as.integer(EUEindic)+1
	wc <- ifelse(EUEindic,"wagechange_EUE","wagechange")
	wt <- ifelse(EUEindic,"balanceweightEUE","balanceweight")
	
	for( si in seq(1,NS-1)){
		for(ki in seq(si+1,NS)){
			kshere <- ks.test( DTseam[get(paste0("s",si))==T,eval(as.name(wc)) ] , DTseam[get(paste0("s",ki))==T,eval(as.name(wc))])
			tab_chngdist_ks[si,ki-1] = kshere$p.value
		}
	}
}

tab_chngdist_moodqtl     <- array(0.,dim=c(NS,NS,length(tabqtls),2))
tab_chngdistEE_moodqtl   <- array(0.,dim=c(NS,NS,length(tabqtls)))
tab_chngdistEUUE_moodqtl <- array(0.,dim=c(NS,NS,length(tabqtls),2))
for(EUEindic in c(F,T)){
	idx = as.integer(EUEindic)+1
	wc <- ifelse(EUEindic,"wagechange_EUE","wagechange")
	wt <- ifelse(EUEindic,"balanceweightEUE","balanceweight")
	#full dist				 
	for( qi in seq(1,length(tabqtls)) ){
		for( si in seq(1,NS)){
			for(ki in seq(1,NS)){
				if(si != ki){
					Finvq <- tab_chngdist_rec[si,qi+1]
					Nbase <- DTseam[ get(paste0("s",si))==T, sum(eval(as.name(wc))>0, na.rm=T) ]
					prob_other <- DTseam[ get(paste0("s",ki))==T, wtd.mean(eval(as.name(wc)) < Finvq, weights=eval(as.name(wt)), na.rm=T) ]
					test_tab <- cbind( c( tabqtls[qi], 1.-tabqtls[qi] ), c(prob_other,1.-prob_other) )*Nbase
					tab_chngdist_moodqtl[si,ki,qi,idx] <- chisq.test(test_tab)$p.value 	
				}
			}
		}
	}
	# EUE and EUUE
	for( qi in seq(1,length(tabqtls)) ){
		for( si in seq(1,NS)){
			for(ki in seq(1,NS)){
				if(si != ki){
					if(EUEindic){
						Finvq <- tab_chngdistEUE_recEUE[si,qi+1]
					}else{
						Finvq <- tab_chngdist_recEUUE[si,qi+1]
					}
					Nbase <- DTseam[ get(paste0("s",si))==T, sum(eval(as.name(wc))>0, na.rm=T) ]
					prob_other <- DTseam[ get(paste0("s",ki))==T, wtd.mean(eval(as.name(wc)) < Finvq, weights=eval(as.name(wt)), na.rm=T) ]
					test_tab <- cbind( c( tabqtls[qi], 1.-tabqtls[qi] ), c(prob_other,1.-prob_other) )*Nbase
					tab_chngdistEUUE_moodqtl[si,ki,qi,idx] <- chisq.test(test_tab)$p.value 	
				}
			}
		}
	}
}
# EE
for( qi in seq(1,length(tabqtls)) ){
	for( si in seq(1,NS)){
		for(ki in seq(1,NS)){
			if(si != ki){
				Finvq <- tab_chngdist_recEE[si,qi+1]
				Nbase <- DTseam[ get(paste0("s",si))==T, sum(wagechange>0, na.rm=T) ]
				prob_other <- DTseam[ get(paste0("s",ki))==T, wtd.mean(wagechange < Finvq, weights=balanceweight, na.rm=T) ]
				test_tab <- cbind( c( tabqtls[qi], 1.-tabqtls[qi] ), c(prob_other,1.-prob_other) )*Nbase
				tab_chngdistEE_moodqtl[si,ki,qi] <- chisq.test(test_tab)$p.value 	
			}
		}
	}
}

#output them to latex tables
rnames <- c("All, 1996-2012","Stayers, 1996-2012","Changers, 1996-2012",
			"All, Expansion","Stayers, Expansion","Changers, Expansion",
			"All, Recession","Stayers, Recession","Changers, Recession")
#cnames <- c("All, 96-12","Stay, 96-12","Chng, 96-12",
#			"All, Exp","Stay, Exp","Chng, Exp",
#			"All, Rec","Stay, Rec","Chng, Rec")
cnames <- c("All    ","Stay    ","Chng    ",
			"All\   ","Stay\   ","Chng\   ",
			"All\ \ ","Stay\ \ ","Chng\ \ ")
rowtitles = list( pos=list(-1), command = c("& \\multicolumn{3}{|c|}{1996-2012} & \\multicolumn{3}{|c|}{Expansion} &\\multicolumn{3}{|c}{Recession} \\\\ \n  "))
for(qi in seq(1,length(tabqtls))) {
	for(EUEindic in c(F,T)){
		EUEtxt <- ifelse(EUEindic," connecting across unemployment spells ","")
		EUEfnam <- ifelse(EUEindic,"EUE","")
		tab_chngdist_mood <- data.table(tab_chngdist_moodqtl[,,qi,EUEindic+1])
		names(tab_chngdist_mood) <- cnames
		rownames(tab_chngdist_mood) <- rnames
		tab_chngdist_mood <- xtable(tab_chngdist_mood, digits=4, 
									align="l|lll|lll|lll", caption=paste0("Difference in the ", tabqtls[qi]*100 ," pctile of earnings changes among job changers",EUEtxt,", p-values \\label{tab:chng_moodqtl",qi,"}"))
		print(tab_chngdist_mood,include.rownames=T, hline.after= c(0,nrow(tab_chngdist_mood)), add.to.row = rowtitles, 
			  file=paste0("./Figures/chng",EUEfnam,"_moodqtl",qi,".tex"))
		#EUE
		tab_chngdist_mood <- data.table(tab_chngdistEUUE_moodqtl[,,qi,EUEindic+1])
		names(tab_chngdist_mood) <- cnames
		rownames(tab_chngdist_mood) <- rnames
		tab_chngdist_mood <- xtable(tab_chngdist_mood, digits=4, 
									align="l|lll|lll|lll", caption=paste0("Difference in the ", tabqtls[qi]*100 ," pctile of earnings changes among job changers through unemployment",EUEtxt,", p-values \\label{tab:chng_moodqtl",qi,"}"))
		print(tab_chngdist_mood,include.rownames=T, hline.after= c(0,nrow(tab_chngdist_mood)), add.to.row = rowtitles, 
			  file=paste0("./Figures/chng",EUEfnam,"_EUUE_moodqtl",qi,".tex"))
	}
	#EE
	tab_chngdist_mood <- data.table(tab_chngdistEE_moodqtl[,,qi])
	names(tab_chngdist_mood) <- cnames
	rownames(tab_chngdist_mood) <- rnames
	tab_chngdist_mood <- xtable(tab_chngdist_mood, digits=4, 
								align="l|lll|lll|lll", caption=paste0("Difference in the ", tabqtls[qi]*100 ," pctile of earnings changes among job-to-job changers, p-values \\label{tab:chng_moodqtl",qi,"}"))
	print(tab_chngdist_mood,include.rownames=T, hline.after= c(0,nrow(tab_chngdist_mood)), add.to.row = rowtitles, 
		  file=paste0("./Figures/chng_EE_moodqtl",qi,".tex"))
}


DTseam[ , c("s1","s2","s3","s4","s5","s6","s7","s8","s9"):=NULL]

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