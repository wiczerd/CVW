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


toKeep <- c("switchedOcc","switchedInd",
			"Young",
			"HSCol",
			"recIndic",
			"wagechange",
			"wagechange_EUE", 
			"balanceweight","balanceweightEUE",
			"wagechange_all",
			"EE","EU","UE",
			"unrt",
			"wave","id")

# select toKeep columns only
wagechanges <- wagechanges[, toKeep, with = FALSE]
wagechanges<- wagechanges[ is.finite(switchedOcc), ]


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
            "unrt","wpfinwgt",
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

tab_wavedist[1,1]    <- DTseam[                                   ,     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
tab_wavedist[1,2:tN] <- DTseam[                                   , wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
tab_wavedist[2,1]    <- DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
tab_wavedist[2,2:tN] <- DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
tab_wavedist[3,1]    <- DTseam[  EU_wave==T|UE_wave==T|EE_wave==T ,     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
tab_wavedist[3,2:tN] <- DTseam[  EU_wave==T|UE_wave==T|EE_wave==T , wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
#expansion/recession
for(rI in c(T,F)){
  rix = rI*3+3
  tab_wavedist[1+rix,1]   <- DTseam[recIndic_wave == rI                       ,     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
  tab_wavedist[1+rix,2:tN]<- DTseam[recIndic_wave == rI                       , wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
  tab_wavedist[2+rix,1]   <- DTseam[recIndic_wave == rI & !(EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
  tab_wavedist[2+rix,2:tN]<- DTseam[recIndic_wave == rI & !(EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
  tab_wavedist[3+rix,1]   <- DTseam[recIndic_wave == rI &  (EU_wave==T|UE_wave==T|EE_wave==T),     wtd.mean(wagechange_wave,na.rm=T,weights=wpfinwgt)]
  tab_wavedist[3+rix,2:tN]<- DTseam[recIndic_wave == rI &  (EU_wave==T|UE_wave==T|EE_wave==T), wtd.quantile(wagechange_wave,na.rm=T,weights=wpfinwgt, probs=tabqtls)]
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


##########################################################################################
# Full sample table-------------------------------------------------------------
DTall <- readRDS("./Data/DTall_6.RData")
DTall <- merge(DTall, CPSunempRt, by = "date", all.x = TRUE)

toKeep <- c(toKeep,"wpfinwgt","switchedJob","allwt","allwtEUE","balanceweight","balanceweightEUE",
            "wagechange_seam","wagechange_allEUE")


# select toKeep columns only
DTall <- DTall[, toKeep, with = FALSE]
DTall <- subset(DTall, is.finite(wpfinwgt) & is.finite(wagechange_all))

DTall<-DTall[ is.finite(EE)&is.finite(EU)&is.finite(UE),]
tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)

tab_fulldist <- array(0., dim=c(9,length(tabqtls)+1,2))
for(EUEindic in seq(T,F)){
	wc <- ifelse(EUEindic,"wagechange_allEUE","wagechange_all")
	wt <- ifelse(EUEindic,"allwtEUE","allwt")
	idx = as.integer(EUEindic)+1
	tab_fulldist[1,1,idx]    <- DTall[                    ,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_fulldist[1,2:tN,idx] <- DTall[                    , wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_fulldist[2,1,idx]    <- DTall[!(EU==T|UE==T|EE==T),     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_fulldist[2,2:tN,idx] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_fulldist[3,1,idx]    <- DTall[  EU==T|UE==T|EE==T ,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_fulldist[3,2:tN,idx] <- DTall[  EU==T|UE==T|EE==T , wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	#expansion/recession
	for(rI in c(T,F)){
		rix = rI*3+3
		tab_fulldist[1+rix,1,idx]   <- DTall[recIndic == rI                       ,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_fulldist[1+rix,2:tN,idx]<- DTall[recIndic == rI                       , wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_fulldist[2+rix,1,idx]   <- DTall[recIndic == rI & !(EU==T|UE==T|EE==T),     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_fulldist[2+rix,2:tN,idx]<- DTall[recIndic == rI & !(EU==T|UE==T|EE==T), wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_fulldist[3+rix,1,idx]   <- DTall[recIndic == rI &  (EU==T|UE==T|EE==T),     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_fulldist[3+rix,2:tN,idx]<- DTall[recIndic == rI &  (EU==T|UE==T|EE==T), wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}
}


#output it to tables
dat_fulldist <- tab_fulldist[,,1]
tab_fulldist <- data.table(tab_fulldist)
names(tab_fulldist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rnames <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
			"All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
			"All\ Workers\ \ ", "Same\ Job,\ \ ","Chng\ Job\ \ ")
rownames(tab_fulldist) <- rnames
dat_fulldistEUE <- tab_fulldist[,,2]
tab_fulldistEUE <- data.table(tab_fulldistEUE)
names(tab_fulldistEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_fulldistEUE) <- rnames

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
					"\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
					"\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
tab_fulldist <- xtable(tab_fulldist, digits=2, 
					align="l|l|lllll", caption="Distribution of earnings changes \\label{tab:fulldist}")
print(tab_fulldist,include.rownames=T, hline.after= c(nrow(tab_fulldist)), 
	  add.to.row=rowtitles, file="./Figures/fulldist.tex")


tab_fulldistEUE <- xtable(tab_fulldistEUE, digits=2,
					   align="l|l|lllll", caption="Distribution of earnings changes, connecting unemployment spells \\label{tab:fulldistEUE}")
print(tab_fulldistEUE,include.rownames=T, hline.after= c(nrow(tab_fulldistEUE)),
	  add.to.row= rowtitles,file="./Figures/fulldistEUE.tex")

# are these different? --------------------------------

# first do KS test with each subsample 
DTall[ , s1:= T]
DTall[ , s2:= ifelse(!(EU==T|UE==T|EE==T), T,F) ]
DTall[ , s3:= ifelse(EU==T|UE==T|EE==T, T,F)]
#expansion
DTall[ , s4:= ifelse(recIndic == F,T,F)]
DTall[ , s5:= ifelse(recIndic == F & !(EU==T|UE==T|EE==T),T,F)]
DTall[ , s6:= ifelse(recIndic == F &  (EU==T|UE==T|EE==T),T,F)]
#recession
DTall[ , s7 := ifelse(recIndic == T, T,F)]
DTall[ , s8 := ifelse(recIndic == T & !(EU==T|UE==T|EE==T),T,F)]
DTall[ , s9 := ifelse(recIndic == T &  (EU==T|UE==T|EE==T),T,F)]

NS = 9

tab_fulldist_ks <- array(NA_real_,dim=c(NS-1,NS-1))
for( si in seq(1,NS-1)){
	for(ki in seq(si+1,NS)){
		kshere <- ks.test( DTall[get(paste0("s",si))==T,wagechange_all ] , DTall[get(paste0("s",ki))==T,wagechange_all])
		tab_fulldist_ks[si,ki-1] = kshere$p.value
	}
}

tab_fulldist_moodqtl <- array(0.,dim=c(NS,NS,length(tabqtls)))
for( qi in seq(1,length(tabqtls)) ){
	for( si in seq(1,NS)){
		for(ki in seq(1,NS)){
			if(si != ki){
				Finvq <- dat_fulldist[si,qi+1]
				Nbase <- DTall[ get(paste0("s",si))==T, sum(allwt>0, na.rm=T) ]
				prob_other <- DTall[ get(paste0("s",ki))==T, wtd.mean(wagechange_all < Finvq, weights=allwt, na.rm=T) ]
				test_tab <- cbind( c( tabqtls[qi], 1.-tabqtls[qi] ), c(prob_other,1.-prob_other) )*Nbase
				tab_fulldist_moodqtl[si,ki,qi] <- chisq.test(test_tab)$p.value 	
			}
		}
	}
}

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
	tab_fulldist_mood <- data.table(tab_fulldist_moodqtl[,,qi])
	names(tab_fulldist_mood) <- cnames
	rownames(tab_fulldist_mood) <- rnames
	tab_fulldist_mood <- xtable(tab_fulldist_mood, digits=4, 
						   align="l|lll|lll|lll", caption=paste0("Difference in the ", tabqtls[qi]*100 ," pctile of earnings changes, p-values \\label{tab:moodqtl",qi,"}"))
	print(tab_fulldist_mood,include.rownames=T, hline.after= c(0,nrow(tab_fulldist_mood)), add.to.row = rowtitles, 
		  file=paste0("./Figures/tab_moodqtl",qi,".tex"))
}
tab_fulldistEUE_moodqtl <- array(0.,dim=c(NS,NS,length(tabqtls)))
for( qi in seq(1,length(tabqtls)) ){
	for( si in seq(1,NS)){
		for(ki in seq(1,NS)){
			if(si != ki){
				Finvq <- dat_fulldistEUE[si,qi+1]
				Nbase <- DTall[ get(paste0("s",si))==T, sum(allwt>0, na.rm=T) ]
				prob_other <- DTall[ get(paste0("s",ki))==T, wtd.mean(wagechange_allEUE < Finvq, weights=allwt, na.rm=T) ]
				test_tab <- cbind( c( tabqtls[qi], 1.-tabqtls[qi] ), c(prob_other,1.-prob_other) )*Nbase
				tab_fulldistEUE_moodqtl[si,ki,qi] <- chisq.test(test_tab)$p.value 	
			}
		}
	}
}

for(qi in seq(1,length(tabqtls))) {
	tab_fulldistEUE_mood <- data.table(tab_fulldistEUE_moodqtl[,,qi])
	names(tab_fulldistEUE_mood) <- cnames
	rownames(tab_fulldistEUE_mood) <- rnames
	tab_fulldistEUE_mood <- xtable(tab_fulldistEUE_mood, digits=4, 
								align="l|lll|lll|lll", caption=paste0("Difference in the ", tabqtls[qi]*100 ," pctile of earnings changes connecting unemployment spells, p-values \\label{tab:moodqtl",qi,"}"))
	print(tab_fulldistEUE_mood,include.rownames=T, hline.after= c(0,nrow(tab_fulldistEUE_mood)), add.to.row = rowtitles, 
		  file=paste0("./Figures/tab_EUE_moodqtl",qi,".tex"))
}


DTall[ , c("s1","s2","s3","s4","s5","s6","s7","s8","s9") := NULL]

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



tab_vardec <- array(NA_real_,dim=c(4,4))

totmean <- DTall[, wtd.mean(wagechange_noseam,na.rm=T,weights=allwt) ]
totvar  <- DTall[, sum(allwt*(wagechange_noseam- totmean)^2,na.rm=T) ]
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

tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
					   align="l|l|lll", caption="Decomposition of earnings change dispersion \\label{tab:chngvarqtldec}")
print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file="./Figures/chngvarqtldec.tex")


# recession and expansion
for (rI in c(F,T)){
	for(EUEindic in c(F,T)){
		wc <- ifelse(EUEindic,"wagechange_noseamEUE","wagechange_noseam")
		wt <- ifelse(EUEindic,"allwtEUE","allwt")
		
		totmean <- DTall[recIndic==rI, wtd.mean(wagechange_noseam,na.rm=T,weights=allwt) ]
		totvar  <- DTall[recIndic==rI, sum(allwt*(wagechange_noseam- totmean)^2,na.rm=T) ]
		idx = 2*as.integer(EUEindic)
		tab_vardec[1+idx,1] <- DTall[(EE==F&EU==F&UE==F) & recIndic==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
		tab_vardec[1+idx,2] <- DTall[(EE==T&EU==F&UE==F) & recIndic==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
		tab_vardec[1+idx,3] <- DTall[(EE==F&(EU==T|UE==T))&recIndic==rI, sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
		totwt <- DTall[recIndic==rI, sum(eval(as.name(wt)),na.rm=T) ]
		tab_vardec[2+idx,1] <- DTall[(EE==F&EU==F&UE==F) & recIndic==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt
		tab_vardec[2+idx,2] <- DTall[(EE==T&EU==F&UE==F) & recIndic==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt
		tab_vardec[2+idx,3] <- DTall[(EE==F&(EU==T|UE==T))&recIndic==rI, sum(eval(as.name(wt)),na.rm=T) ]/totwt
		
	}
	
	# Full sample quantile-diff decomposition --------------------------
	tot51025qtl <- DTall[ recIndic==rI, wtd.quantile(wagechange_noseam,na.rm=T,weights=allwt, probs=c(0.05,0.1,.25,.75,.9,.95)) ]
	tot51025qtlEUE <- DTall[ recIndic==rI, wtd.quantile(wagechange_noseamEUE,na.rm=T,weights=allwtEUE, probs=c(0.05,0.1,.25,.75,.9,.95)) ]
	
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
			tab_qtldec[rri+ri,1] <- DTall[recIndic==rI & (EE==F&EU==F&UE==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
			tab_qtldec[rri+ri,2] <- DTall[recIndic==rI & (EE==T&EU==F&UE==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
			tab_qtldec[rri+ri,ci]<- DTall[recIndic==rI & (EE==F&(EU==T|UE==T)) & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
		}
		tab_qtldec[Ndifs+1+ri,1] <- DTall[recIndic==rI & (EE==F&EU==F&UE==F)   , sum(eval(as.name(wt)),na.rm=T) ]
		tab_qtldec[Ndifs+1+ri,2] <- DTall[recIndic==rI & (EE==T&EU==F&UE==F)   , sum(eval(as.name(wt)),na.rm=T) ]
		tab_qtldec[Ndifs+1+ri,ci]<- DTall[recIndic==rI & (EE==F&(EU==T|UE==T)) , sum(eval(as.name(wt)),na.rm=T) ]
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
	
	if(rI){
		tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
								align="l|l|lll", caption="Decomposition of earnings change dispersion, recession \\label{tab:chngvarqtldec_rec}")
	}else{
		tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
									align="l|l|lll", caption="Decomposition of earnings change dispersion, expansion \\label{tab:chngvarqtldec_exp}")
	}
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,Ndifs+1, nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0("./Figures/chngvarqtldec_rec",rI,".tex"))
	
}


##########################################################################################

# Only job-changers -----------------------
tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)

tab_chngdist    <- array(0.,dim=c(9,tN,2))

for( EUEindic in c(F,T)){
	wc <- ifelse(EUEindic,"wagechange_EUE","wagechange")
	wt <- ifelse(EUEindic,"balanceweightEUE","balanceweight")
	idx = 1+as.integer(EUEindic)
	tab_chngdist[1,1,idx]    <- wagechanges[(EE|EU|UE),                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist[1,2:tN,idx] <- wagechanges[ EE|EU|UE ,                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_chngdist[2,1,idx]    <- wagechanges[(EE|EU|UE)& switchedOcc==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist[2,2:tN,idx] <- wagechanges[(EE|EU|UE)& switchedOcc==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_chngdist[3,1,idx]    <- wagechanges[(EE|EU|UE)& switchedOcc==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist[3,2:tN,idx] <- wagechanges[(EE|EU|UE)& switchedOcc==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	for(EUhere in c(F,T)){
		EEhere = ifelse(EUhere==F, T,F)
		eidx = as.integer(EUhere)*3
		tab_chngdist[4+eidx,1,idx]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)),                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist[4+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere)),                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
		tab_chngdist[5+eidx,1,idx]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist[5+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_chngdist[6+eidx,1,idx]    <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist[6+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))& switchedOcc==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}

}
dat_chngdist <- tab_chngdist

#output it to tables
tab_chngdist_out <- data.table(tab_chngdist[1:3,,1])
names(tab_chngdist_out) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngdist_out) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")
tab_chngdistEUE_out <- data.table(tab_chngdist[1:3,,2])
names(tab_chngdistEUE_out) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_chngdistEUE_out) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")

tab_chngdist_out <- xtable(tab_chngdist_out, label="tab:chngdist", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job changers")
print(tab_chngdist_out,include.rownames=T, hline.after= c(0,nrow(tab_chngdist_out)), file="./Figures/chngdist.tex")

tab_chngdistEUE_out <- xtable(tab_chngdistEUE_out, label="tab:chngdistEUE", digits=2, 
						  align="l|l|lllll", caption="Distribution of earnings changes among job changers, connecting unemployment spells")
print(tab_chngdistEUE_out,include.rownames=T, hline.after= c(0,nrow(tab_chngdistEUE_out)), file="./Figures/chngdistEUE.tex")


tab_chngdist_rec    <- array(NA,dim=c(6,tN,2))
tab_chngdist_rec_EEEU    <- array(NA,dim=c(12,tN,2))

for(EUEindic in c(F,T)){
for(recHere in c(F,T)){
	ridx = as.integer(recHere)*3
	idx = as.integer(EUEindic) +1
	wc <- ifelse(EUEindic,"wagechange_EUE","wagechange")
	wt <- ifelse(EUEindic,"balanceweightEUE","balanceweight")
	
	tab_chngdist_rec[1+ridx,1   ,idx] <- wagechanges[(EE|EU|UE)                & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist_rec[1+ridx,2:tN,idx] <- wagechanges[(EE|EU|UE)                & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_chngdist_rec[2+ridx,1   ,idx] <- wagechanges[(EE|EU|UE)&switchedOcc==F & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist_rec[2+ridx,2:tN,idx] <- wagechanges[(EE|EU|UE)&switchedOcc==F & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_chngdist_rec[3+ridx,1   ,idx] <- wagechanges[(EE|EU|UE)&switchedOcc==T & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_chngdist_rec[3+ridx,2:tN,idx] <- wagechanges[(EE|EU|UE)&switchedOcc==T & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]

	for(EUhere in c(F,T)){
		EEhere = ifelse(EUhere==F, T,F)
		eidx = 6*as.integer(EUhere)
		tab_chngdist_rec_EEEU[1+ridx+eidx,1,   idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))                &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist_rec_EEEU[1+ridx+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))                &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
		tab_chngdist_rec_EEEU[2+ridx+eidx,1,   idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))&switchedOcc==F &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist_rec_EEEU[2+ridx+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))&switchedOcc==F &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_chngdist_rec_EEEU[3+ridx+eidx,1,   idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))&switchedOcc==T &recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_chngdist_rec_EEEU[3+ridx+eidx,2:tN,idx] <- wagechanges[(EE==EEhere&(EU==EUhere|UE==EUhere))&switchedOcc==T &recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}
}
}
#output it to tables
dat_chngdist_rec <-tab_chngdist_rec
dat_chngdist_rec_EEEU <-tab_chngdist_rec_EEEU
cnames <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rnames <- c("All job changers    ", "Occ stayers    ","Occ movers    ",
			"All job changers\   ", "Occ stayers\   ","Occ movers\   ",
			"All job changers\ \ ", "Occ stayers\ \ ","Occ movers\ \ ")

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
											  "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
											  "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )


tab_chngdist_rec <- data.table(rbind(dat_chngdist[1:3,,1],dat_chngdist_rec[1:6,,1]))
names(tab_chngdist_rec) <- cnames
rownames(tab_chngdist_rec) <- rnames
tab_chngdist_rec <- xtable(tab_chngdist_rec, digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job changers in recession and expansion \\label{tab:chngdist_rec}")
print(tab_chngdist_rec,include.rownames=T, hline.after= c(nrow(tab_chngdist_rec)), 
	  add.to.row=rowtitles, file="./Figures/chngdist_rec.tex")


tab_chngdistEUE_rec <- data.table(rbind(dat_chngdist[1:3,,2],dat_chngdist_rec[1:6,,2]))
names(tab_chngdistEUE_rec) <- cnames
rownames(tab_chngdistEUE_rec) <-rnames

tab_chngdistEUE_rec <- xtable(tab_chngdistEUE_rec, digits=2, 
						  align="l|l|lllll", caption="Distribution of earnings changes among job changers in recession and expansion, connecting unemployment spells \\label{tab:chngdistEUE_rec}")
print(tab_chngdistEUE_rec,include.rownames=T, hline.after= c(nrow(tab_chngdistEUE_rec)), 
	  add.to.row=rowtitles, file="./Figures/chngdistEUE_rec.tex")

tab_chngdistEUE_recEE  <- data.table(rbind(dat_chngdist[4:6,,2],dat_chngdist_rec_EEEU[1:6 ,,2]))
tab_chngdistEUE_recEUE <- data.table(rbind(dat_chngdist[7:9,,2],dat_chngdist_rec_EEEU[7:12,,2]))
tab_chngdist_recEE     <- data.table(rbind(dat_chngdist[4:6,,1],dat_chngdist_rec_EEEU[1:6 ,,1]))
tab_chngdist_recEUUE   <- data.table(rbind(dat_chngdist[7:9,,1],dat_chngdist_rec_EEEU[7:12,,1]))

names(tab_chngdistEUE_recEUE) <- cnames
names(tab_chngdistEUE_recEE) <- cnames
names(tab_chngdist_recEUUE) <- cnames
names(tab_chngdist_recEE) <- cnames

rownames(tab_chngdistEUE_recEE) <-rnames
rownames(tab_chngdistEUE_recEUE) <-rnames
rownames(tab_chngdist_recEE) <-rnames
rownames(tab_chngdist_recEUUE) <-rnames


tab_chngdist_recEE <- xtable(tab_chngdist_recEE, label="tab:chngdist_recEE", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes among job-job changers")
print(tab_chngdist_recEE,include.rownames=T, hline.after= c(nrow(tab_chngdist_recEE)), 
	  add.to.row=rowtitles, file="./Figures/chngdist_recEE.tex")
tab_chngdistEUE_recEUE <- xtable(tab_chngdistEUE_recEUE, label="tab:chngdistEUE_recEE", digits=2, 
								align="l|l|lllll", caption="Distribution of earnings changes among those transitioning through unemployment")
print(tab_chngdistEUE_recEUE,include.rownames=T, hline.after= c(nrow(tab_chngdistEUE_recEUE)), 
	  add.to.row=rowtitles, file="./Figures/chngdistEUE_recEUE.tex")
tab_chngdist_recEUUE <- xtable(tab_chngdist_recEUUE, label="tab:chngdistEUE_recEE", digits=2, 
								 align="l|l|lllll", caption="Distribution of earnings changes among those transitioning through unemployment")
print(tab_chngdist_recEUUE,include.rownames=T, hline.after= c(nrow(tab_chngdist_recEUUE)), 
	  add.to.row=rowtitles, file="./Figures/chngdist_recEUUE.tex")

# Are these the same distributions ? -----------------------------------



# first do KS test with each subsample 
wagechanges[ , s1:= T]
wagechanges[ , s2:= ifelse(switchedOcc==F, T,F) ]
wagechanges[ , s3:= ifelse(switchedOcc==T, T,F)]
#expansion
wagechanges[ , s4:= ifelse(recIndic == F,T,F)]
wagechanges[ , s5:= ifelse(recIndic == F & switchedOcc==F,T,F)]
wagechanges[ , s6:= ifelse(recIndic == F & switchedOcc==T,T,F)]
#recession
wagechanges[ , s7 := ifelse(recIndic == T, T,F)]
wagechanges[ , s8 := ifelse(recIndic == T & switchedOcc==F,T,F)]
wagechanges[ , s9 := ifelse(recIndic == T & switchedOcc==T,T,F)]

NS = 9

tab_chngdist_ks <- array(NA_real_,dim=c(NS-1,NS-1,2))
for(EUEindic in c(F,T)){
	idx = as.integer(EUEindic)+1
	wc <- ifelse(EUEindic,"wagechange_EUE","wagechange")
	wt <- ifelse(EUEindic,"balanceweightEUE","balanceweight")
	
	for( si in seq(1,NS-1)){
		for(ki in seq(si+1,NS)){
			kshere <- ks.test( wagechanges[get(paste0("s",si))==T,eval(as.name(wc)) ] , wagechanges[get(paste0("s",ki))==T,eval(as.name(wc))])
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
					Nbase <- wagechanges[ get(paste0("s",si))==T, sum(eval(as.name(wc))>0, na.rm=T) ]
					prob_other <- wagechanges[ get(paste0("s",ki))==T, wtd.mean(eval(as.name(wc)) < Finvq, weights=eval(as.name(wt)), na.rm=T) ]
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
					Nbase <- wagechanges[ get(paste0("s",si))==T, sum(eval(as.name(wc))>0, na.rm=T) ]
					prob_other <- wagechanges[ get(paste0("s",ki))==T, wtd.mean(eval(as.name(wc)) < Finvq, weights=eval(as.name(wt)), na.rm=T) ]
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
				Nbase <- wagechanges[ get(paste0("s",si))==T, sum(wagechange>0, na.rm=T) ]
				prob_other <- wagechanges[ get(paste0("s",ki))==T, wtd.mean(wagechange < Finvq, weights=balanceweight, na.rm=T) ]
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


wagechanges[ , c("s1","s2","s3","s4","s5","s6","s7","s8","s9"):=NULL]

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