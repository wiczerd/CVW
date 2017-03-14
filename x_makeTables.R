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
wtd.Moors <- function(xt, wt){  
	wt  <- wt[is.na(xt)==F]
	xt  <- xt[is.na(xt)==F]
	octls<-wtd.quantile(xt,weights=wt,probs = c(seq(1/8,3/8,1/8),seq(5/8,7/8,1/8)),na.rm = T) #excludes the median
	((octls[6]-octls[4]) - (octls[3]-octls[1]))/(octls[5]-octls[2]) - 1.23
}
wtd.kurtosis <- function(xt, wt){  
	wt  <- wt[is.na(xt)==F]
	xt      <- xt[is.na(xt)==F]
	(sum( wt*(xt-wtd.mean(xt,weights=wt))^4 ) / sum(wt) ) / wtd.var(xt,weights=wt)^(2) - 3.
}


CPSunempRt <- readRDS("./InputData/CPSunempRt.RData")
CPSunempRt$unrt <- CPSunempRt$unrt/100

recDef <- "recIndic_wave"
wt <- "truncweight"
wc <- "wagechangeEUE_wave"

demolbl <- 0 #or choose number from categories in demotxt
demotxt <- c("Young", "Prime","Old","HS","Col","Male","Female")

bootse <- F #compute bootstrapped standard errors or no?
seedint = 941987

##########################################################################################
# By wave -----------------------
toKeep_wave <- c("switchedOcc_wave",
            "ageGrp","HSCol",
            "recIndic","recIndic_wave","recIndic2_wave","recIndic_stint",
            "wagechange_month","wagechange_wave","wagechangeEUE_wave",
            "wagechange_wave_bad","wagechange_wave_bad2","wagechange_wave_low","wagechange_wave_high","wagechange_wave_jcbad",
            "EE_wave","EU_wave","UE_wave","changer","stayer",
            "unrt","wpfinwgt","perwt","cycweight","truncweight","cleaningtruncweight",
			"lfstat_wave","next.lfstat_wave","wave","id","date","panel")
DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))
DTseam <- merge(DTseam, CPSunempRt, by = "date", all.x = TRUE)

# select toKeep columns only
DTseam <- DTseam[, toKeep_wave, with = FALSE]
DTseam <- subset(DTseam, is.finite(wpfinwgt) & is.finite(wagechange_wave))
DTseam<-DTseam[ is.finite(EE_wave)&is.finite(EU_wave)&is.finite(UE_wave), ]


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

# how to weights EUE's? 2x for an EUUE?
#if(wc == "wagechangeEUE_wave"){
#	DTseam[ UE_wave==T,eval(as.name(wt)):= 0.]
#}

#DTseam <- DTseam[ (stayer|changer) & demo==T, ]

wc_wave <- DTseam[ recIndic_wave==T& stayer==T, .(wc_wave = weighted.mean(wagechange_wave, wpfinwgt, na.rm = TRUE)), by = date]
ggplot(wc_wave, aes(date, wc_wave)) +
	geom_point() +
	geom_line() +xlab("") + ylab("mean wage change, stayers wave-frequency")

wc_wave <- DTseam[ stayer==T, .(wc_wave = wtd.mad(wagechange_wave, wpfinwgt)), by = list(date,panel)]
ggplot(wc_wave, aes(date, wc_wave,color=panel,group=panel)) +
	geom_point() + geom_line() +
	ggtitle("Wage Growth Disperison, stayers w/ cleaning")+xlab("") + ylab("median abs dev(wage change), stayers wave-frequency")
#ggsave("mad_wagegrowth_clean.png",height=5,width=10)
#ggsave("mad_wagegrowth_clean.eps",height=5,width=10)

wc_wave <- DTseam[ EE_wave==F & lfstat_wave==1 & next.lfstat_wave==1, .(wc_wave = wtd.mad(wagechange_wave, wpfinwgt)), by = list(date,panel)]
ggplot(wc_wave, aes(date, wc_wave,color=panel,group=panel)) +
	geom_point() + ggtitle("Wage Growth Disperison, stayers no cleaning")+ ylim(c(0,.2))+
	geom_line() +xlab("") + ylab("median abs dev(wage change), stayers wave-frequency")
#ggsave("mad_wagegrowth_noclean.png",height=5,width=10)
#ggsave("mad_wagegrowth_noclean.eps",height=5,width=10)


# wage quantiles ---------------------------------------------------------------
tabqtls <- c(.1,.25,.5,.75,.9)
tN <- (length(tabqtls)+1)
tab_wavedist <- array(0., dim=c(9,length(tabqtls)+1))

if(bootse == T){
	set.seed(seedint)
	#draw the sample
	Nsim = 50
	nsampE = nrow(DTseam[eval(as.name(recDef)) == F ])
	nsampR = nrow(DTseam[eval(as.name(recDef)) == T ])
	nsamp  = nsampR+nsampE
	se_wavedist <- array(0.,dim = c(9,length(tabqtls)+1,Nsim))
}else{
	Nsim = 0
}

for( si in seq(1,bootse*Nsim+1) ){
	if(si>1){
		DThr <- DTseam[ sample(nsamp,nsamp,replace=T,prob=eval(as.name(wt)))]
	}else{
		DThr <- DTseam
	}
	tab_wavedist[1,1]    <- DThr[(stayer|changer)&demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
	tab_wavedist[1,2:tN] <- DThr[(stayer|changer)&demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavedist[2,1]    <- DThr[  stayer ==T    &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavedist[2,2:tN] <- DThr[  stayer ==T    &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavedist[3,1]    <- DThr[  changer==T    &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavedist[3,2:tN] <- DThr[  changer==T    &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	#expansion/recession
	for(rI in c(T,F)){
		rix = rI*3+3
		tab_wavedist[1+rix,1]   <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavedist[1+rix,2:tN]<- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_wavedist[2+rix,1]   <- DThr[eval(as.name(recDef)) == rI & stayer ==T      &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavedist[2+rix,2:tN]<- DThr[eval(as.name(recDef)) == rI & stayer ==T      &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
		tab_wavedist[3+rix,1]   <- DThr[eval(as.name(recDef)) == rI & changer==T      &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavedist[3+rix,2:tN]<- DThr[eval(as.name(recDef)) == rI & changer==T      &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	}
	if(si>1){
		se_wavedist[,,si-1] = tab_wavedist
	}else{
		dat_wavedist <- tab_wavedist	
	}
}
#output it to tables
tab_wavedist <- dat_wavedist
tab_wavedist <- data.table(tab_wavedist)
names(tab_wavedist) <- c("Mean",as.character( tabqtls))
#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rnames <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
            "All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
            "All\ Workers\ \ ", "Same\ Job\ \ ","Chng\ Job\ \ ")
rownames(tab_wavedist) <- rnames

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
                                              "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
                                              "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )

if(wc == "wagechangeEUE_wave"){
	nametab <- "wavedistEUE"
}else{
	nametab <- "wavedist"
}
tab_wavedist <- xtable(tab_wavedist, digits=2, 
                       align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"}"))
if(demolbl>=1 & demolbl<=7){
	print(tab_wavedist,include.rownames=T, hline.after= c(nrow(tab_wavedist)), 
		add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],".tex"))
}else{
	print(tab_wavedist,include.rownames=T, hline.after= c(nrow(tab_wavedist)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,".tex"))
}
tab_wavedistci <- array(0.,dim=c(2*nrow(tab_wavedist),2*ncol(tab_wavedist)))
tab_wavedistse <- array(0.,dim=c(2*nrow(tab_wavedist),ncol(tab_wavedist)))
if(bootse == T){
	for( ri in seq(0,nrow(tab_wavedist)-1) ){
		for(ci in seq( 0,ncol(tab_wavedist)-1 )){
			tab_wavedistci[ ri*2+1,ci*2+1 ] <- tab_wavedist[ri+1,ci+1]
			tab_wavedistci[ ri*2+2,(ci*2+1):(ci*2+2) ] <- quantile(se_wavedist[ri+1,ci+1, ], probs=c(0.05,0.95))
			tab_wavedistse[ ri*2+1,ci+1 ] <- tab_wavedist[ri+1,ci+1]
			tab_wavedistse[ ri*2+2,ci+1 ] <- var(se_wavedist[ri+1,ci+1, ])^0.5
		}
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
rnames <- c("All\ Workers"    ," "     , "Same\ Job"    ,","      ,"Chng\ Job"    ,",,",
			"All\ Workers\ "  ,",,,"   , "Same\ Job\ "  ,",,,,"   ,"Chng\ Job\ "  ,",,,,,",
			"All\ Workers\ \ ",",,,,,,", "Same\ Job\ \ ",",,,,,,,","Chng\ Job\ \ ",",,,,,,,,")
rownames(tab_wavedistse) <- rnames

rowtitles <- list( pos=list(0,6,12), command=c("\\hline  \\color{Maroon}{1996-2012} &  & & & & & \\\\ \n",
											  "\\hline \\hline   \\color{Maroon}{Expansion} &  & & & & & \\\\  \n", 
											  "\\hline \\hline   \\color{Maroon}{Recession} &  & & & & & \\\\  \n")  )
if(wc == "wagechangeEUE_wave"){
	nametab <- "wavedistseEUE"
}else{
	nametab <- "wavedistse"
}
tab_wavedistse <- xtable(tab_wavedistse, digits=2, 
					   align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"}"))

if(demolbl>=1 & demolbl<=7){
	print(tab_wavedistse,include.rownames=T, hline.after= c(nrow(tab_wavedistse)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],".tex"))
}else{
	print(tab_wavedistse,include.rownames=T, hline.after= c(nrow(tab_wavedistse)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,".tex"))
}
#Wage moments --------------------------------------------------------

if(bootse == T){
	set.seed(seedint)
	#draw the sample
	nsampE = nrow(DTseam[eval(as.name(recDef)) == F ])
	nsampR = nrow(DTseam[eval(as.name(recDef)) == T ])
	nsamp  = nsampR+nsampE
	se_wavemoments <- array(0.,dim = c(9,5,Nsim))
}
tab_wavemoments <- array(0., dim=c(9,5))
for( si in seq(1,bootse*Nsim+1) ){
	if(si>1){
		DThr <- DTseam[ sample(nsamp,nsamp,replace=T,prob=eval(as.name(wt)))]
	}else{
		DThr <- DTseam
	}
	tab_wavemoments[1,1]    <- DThr[(stayer|changer)&demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
	tab_wavemoments[1,2]    <- DThr[(stayer|changer)&demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
	tab_wavemoments[1,3]    <- DThr[(stayer|changer)&demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[1,2])]
	tab_wavemoments[1,4]    <- DThr[(stayer|changer)&demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[1,2])]
	tab_wavemoments[1,5]    <- DThr[(stayer|changer)&demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
	tab_wavemoments[2,1]    <- DThr[  stayer ==T    &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavemoments[2,2]    <- DThr[  stayer ==T    &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
	tab_wavemoments[2,3]    <- DThr[  stayer ==T    &demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[2,2])]
	tab_wavemoments[2,4]    <- DThr[  stayer ==T    &demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[2,2])]
	tab_wavemoments[2,5]    <- DThr[  stayer ==T    &demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
	tab_wavemoments[3,1]    <- DThr[  changer==T    &demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavemoments[3,2]    <- DThr[  changer==T    &demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
	tab_wavemoments[3,3]    <- DThr[  changer==T    &demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[3,2])]
	tab_wavemoments[3,4]    <- DThr[  changer==T    &demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[3,2])]
	tab_wavemoments[3,5]    <- DThr[  changer==T    &demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
	
	#expansion/recession
	for(rI in c(T,F)){
		rix = rI*3+3
		tab_wavemoments[1+rix,1]    <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T,     wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)) )]
		tab_wavemoments[1+rix,2]    <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T, wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
		tab_wavemoments[1+rix,3]    <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[1+rix,2])]
		tab_wavemoments[1+rix,4]    <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[1+rix,2])]
		tab_wavemoments[1+rix,5]    <- DThr[eval(as.name(recDef)) == rI & (stayer|changer)&demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
		tab_wavemoments[2+rix,1]    <- DThr[eval(as.name(recDef)) == rI &   stayer ==T    &demo==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavemoments[2+rix,2]    <- DThr[eval(as.name(recDef)) == rI &   stayer ==T    &demo==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
		tab_wavemoments[2+rix,3]    <- DThr[eval(as.name(recDef)) == rI &   stayer ==T    &demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[2+rix,2])]
		tab_wavemoments[2+rix,4]    <- DThr[eval(as.name(recDef)) == rI &   stayer ==T    &demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[2+rix,2])]
		tab_wavemoments[2+rix,5]    <- DThr[eval(as.name(recDef)) == rI &   stayer ==T    &demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
		tab_wavemoments[3+rix,1]    <- DThr[eval(as.name(recDef)) == rI &   changer==T    &demo==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
		tab_wavemoments[3+rix,2]    <- DThr[eval(as.name(recDef)) == rI &   changer==T    &demo==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=0.5)]
		tab_wavemoments[3+rix,3]    <- DThr[eval(as.name(recDef)) == rI &   changer==T    &demo==T,            wtd.mad(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[3+rix,2])]
		tab_wavemoments[3+rix,4]    <- DThr[eval(as.name(recDef)) == rI &   changer==T    &demo==T,wtd.GroenveldMeeden(eval(as.name(wc)),eval(as.name(wt)),tab_wavemoments[3+rix,2])]
		tab_wavemoments[3+rix,5]    <- DThr[eval(as.name(recDef)) == rI &   changer==T    &demo==T,       wtd.Moors(eval(as.name(wc)),eval(as.name(wt)))]
	}
	if(si>1){
		se_wavemoments[,,si-1] = tab_wavemoments
	}else{
		dat_wavemoments <- tab_wavemoments	
	}
}

#output it to tables
tab_wavemoments<-dat_wavemoments
tab_wavemoments <- data.table(tab_wavemoments)
#names(tab_wavemoments) <- c("Mean","Median","Std Dev", "Skew", "Kurtosis")
names(tab_wavemoments) <- c("Mean","Median","Med Abs Dev", "Groenv-Meeden", "Moors")
#rownames(tab_wavedist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rnames <- c("All\ Workers",      "Same\ Job",     "Chng\ Job",
			"All\ Workers\ ",   "Same\ Job\ ",  "Chng\ Job\ ",
			"All\ Workers\ \ ", "Same\ Job\ \ ","Chng\ Job\ \ ")
rownames(tab_wavemoments) <- rnames

rowtitles <- list( pos=list(0,3,6), command=c("\\hline  \\color{Maroon}{1996-2012} &  \\multicolumn{5}{|c|}{} \\\\ \n",
											  "\\hline \\hline   \\color{Maroon}{Expansion} &  \\multicolumn{5}{|c|}{}  \\\\  \n", 
											  "\\hline \\hline   \\color{Maroon}{Recession} &  \\multicolumn{5}{|c|}{}  \\\\  \n")  )

if(wc == "wagechangeEUE_wave"){
	nametab <- "wavemomentsEUE"
}else{
	nametab <- "wavemoments"
}
tab_wavemoments <- xtable(tab_wavemoments, digits=2, 
					   align="l|lllll", caption=paste0("Moments of earnings change distribution \\label{tab:",nametab,"}"))
if(demolbl>=1 & demolbl<=7){
	print(tab_wavemoments,include.rownames=T, hline.after= c(nrow(tab_wavemoments)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],".tex"))
}else{
	print(tab_wavemoments,include.rownames=T, hline.after= c(nrow(tab_wavemoments)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,".tex"))
}
tab_wavemomentsci <- array(0.,dim=c(2*nrow(tab_wavemoments),2*ncol(tab_wavemoments)))
tab_wavemomentsse <- array(0.,dim=c(2*nrow(tab_wavemoments),ncol(tab_wavemoments)))
if(bootse == T){
	for( ri in seq(0,nrow(tab_wavemoments)-1) ){
		for(ci in seq( 0,ncol(tab_wavemoments)-1 )){
			tab_wavemomentsci[ ri*2+1,ci*2+1 ] <- tab_wavedist[ri+1,ci+1]
			tab_wavemomentsci[ ri*2+2,(ci*2+1):(ci*2+2) ] <- quantile(se_wavemoments[ri+1,ci+1, ], probs=c(0.05,0.95))
			tab_wavemomentsse[ ri*2+1,ci+1 ] <- tab_wavemoments[ri+1,ci+1]
			tab_wavemomentsse[ ri*2+2,ci+1 ] <- var(se_wavemoments[ri+1,ci+1, ])^0.5
		}
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
rnames <- c("All\ Workers"    ," "     , "Same\ Job"    ,","      ,"Chng\ Job"    ,",,",
			"All\ Workers\ "  ,",,,"   , "Same\ Job\ "  ,",,,,"   ,"Chng\ Job\ "  ,",,,,,",
			"All\ Workers\ \ ",",,,,,,", "Same\ Job\ \ ",",,,,,,,","Chng\ Job\ \ ",",,,,,,,,")
rownames(tab_wavemomentsse) <- rnames

rowtitles <- list( pos=list(0,6,12), command=c("\\hline  \\color{Maroon}{1996-2012} &  \\multicolumn{5}{|c|}{} \\\\ \n",
											   "\\hline \\hline   \\color{Maroon}{Expansion} &  \\multicolumn{5}{|c|}{}  \\\\  \n", 
											   "\\hline \\hline   \\color{Maroon}{Recession} &  \\multicolumn{5}{|c|}{}  \\\\  \n")  )
if(wc == "wagechangeEUE_wave"){
	nametab <- "wavemomentsseEUE"
}else{
	nametab <- "wavemomentsse"
}
tab_wavemomentsse <- xtable(tab_wavemomentsse, digits=2, 
						 align="l|l|lllll", caption=paste0("Distribution of earnings changes \\label{tab:",nametab,"}"))

if(demolbl>=1 & demolbl<=7){
	print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,demotxt[demolbl],".tex"))
}else{
	print(tab_wavemomentsse,include.rownames=T, hline.after= c(nrow(tab_wavemomentsse)), 
		  add.to.row=rowtitles, file=paste0(outputdir,"/",nametab,".tex"))
}


#Variance decomp -----------------------------------------------

tab_wavevardec <- array(NA_real_,dim=c(2,3))

totmean <- DTseam[(stayer|changer), wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt))) ]
totvar  <- DTseam[(stayer|changer), sum(eval(as.name(wt))*(eval(as.name(wc))- totmean)^2,na.rm=T) ]
tab_wavevardec[1,1] <- DTseam[stayer == T                        , sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,2] <- DTseam[changer ==T& EE_wave ==T           , sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
tab_wavevardec[1,3] <- DTseam[changer ==T&(EU_wave==T|UE_wave==T), sum(eval(as.name(wt))*(eval(as.name(wc)) - totmean)^2,na.rm=T) ]/totvar
totwt <- DTseam[(stayer|changer), sum(eval(as.name(wt)),na.rm=T) ]
tab_wavevardec[2,1] <- DTseam[stayer == T                        , sum(eval(as.name(wt)),na.rm=T) ]/totwt
tab_wavevardec[2,2] <- DTseam[changer ==T& EE_wave ==T           , sum(eval(as.name(wt)),na.rm=T) ]/totwt
tab_wavevardec[2,3] <- DTseam[changer ==T&(EU_wave==T|UE_wave==T), sum(eval(as.name(wt)),na.rm=T) ]/totwt


# Full sample quantile-diff decomposition --------------------------
tot51025qtl <- DTseam[(stayer==T|changer==T), wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=c(0.05,0.1,.25,.75,.9,.95)) ]

Nqtls <-length(tot51025qtl)
Ndifs <- Nqtls/2

tab_waveqtldec <- array(NA_real_,dim=c(Ndifs+1,3))

qtlshere <- tot51025qtl

for(rri in seq(1,Ndifs)){
  tab_waveqtldec[rri,1] <- DTseam[stayer ==T                          & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  tab_waveqtldec[rri,2] <- DTseam[EE_wave==T & changer==T             & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  tab_waveqtldec[rri,3] <- DTseam[(EU_wave==T|UE_wave==T)& changer==T & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
}
tab_waveqtldec[Ndifs+1,1] <- DTseam[stayer == T                          , sum(eval(as.name(wt)),na.rm=T) ]
tab_waveqtldec[Ndifs+1,2] <- DTseam[EE_wave==T & changer==T              , sum(eval(as.name(wt)),na.rm=T) ]
tab_waveqtldec[Ndifs+1,3] <- DTseam[(EU_wave==T|UE_wave==T)& changer==T  , sum(eval(as.name(wt)),na.rm=T) ]


rsum <- rowSums(tab_waveqtldec, na.rm=T)
for(ri in seq(1,nrow(tab_waveqtldec))){
	tab_waveqtldec[ri,] <-tab_waveqtldec[ri,]/rsum[ri]
}

#output it to tables
tab_chngvarqtldec <- data.table(rbind(tab_wavevardec[1,],tab_waveqtldec[1:Ndifs,]) )
names(tab_chngvarqtldec) <- c("Job\ Stayers","EE","EU,UE")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_chngvarqtldec) <- c("Variance","0.95-0.05","0.9-0.1","0.75-0.25")# "Variance\ ","0.95-0.05\ ","0.9-0.1\ ","0.75-0.25\ ","Pct Sample")

tab_chngvarqtldec <- xtable(tab_chngvarqtldec, digits=2, 
                            align="l|l|ll", caption="Decomposition of earnings change dispersion \\label{tab:wavechngvarqtldec}")
if(demolbl>=1 & demolbl<=7){
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/wavechngvarqtldec", demotxt[demolbl],".tex")) 
}else{
	print(tab_chngvarqtldec,include.rownames=T, hline.after= c(0,nrow(tab_chngvarqtldec)-1, nrow(tab_chngvarqtldec)), file=paste0(outputdir,"/wavechngvarqtldec.tex") ) 
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
    tab_waveqtldec[rri,1] <- DTseam[recIndic_wave==rI & (EE_wave==F&EU_wave==F&UE_wave==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
    tab_waveqtldec[rri,2] <- DTseam[recIndic_wave==rI & (EE_wave==T&EU_wave==F&UE_wave==F)   & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
    tab_waveqtldec[rri,3] <- DTseam[recIndic_wave==rI & (EE_wave==F&(EU_wave==T|UE_wave==T)) & (eval(as.name(wc)) > qtlshere[Nqtls-rri+1] | eval(as.name(wc)) < qtlshere[rri]), sum(eval(as.name(wt)) ,na.rm=T)]
  }
  tab_waveqtldec[Ndifs+1,1] <- DTseam[recIndic_wave==rI &  (stayer==T)                        , sum(eval(as.name(wt)),na.rm=T) ]
  tab_waveqtldec[Ndifs+1,2] <- DTseam[recIndic_wave==rI &  (EE_wave==T&changer==T)            , sum(eval(as.name(wt)),na.rm=T) ]
  tab_waveqtldec[Ndifs+1,3] <- DTseam[recIndic_wave==rI & ((EU_wave==T|UE_wave==T)&changer==T), sum(eval(as.name(wt)),na.rm=T) ]

  
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

#*******************************************************************
#*******************************************************************
# Quantiles among job changers  -----------------------
#tabqtls <- c(.1,.25,.5,.75,.9)
#tN <- (length(tabqtls)+1)

tab_wavechngdist    <- array(0.,dim=c(9,tN))


wc <- "wagechangeEUE_wave"
wt <- "truncweight"
if(wc =="wagechange_wave"){
	labtxt = "wavechngdist"
}else if(wc == "wagechangeEUE_wave"){
	labtxt = "wavechngdistEUE"
}else{
	labtxt = "chgndist"
}


tab_wavechngdist[1,1]    <- DTseam[changer==T ,                    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[1,2:tN] <- DTseam[changer==T ,                wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
tab_wavechngdist[2,1]    <- DTseam[changer==T & switchedOcc_wave==F,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[2,2:tN] <- DTseam[changer==T & switchedOcc_wave==F,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
tab_wavechngdist[3,1]    <- DTseam[changer==T & switchedOcc_wave==T,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
tab_wavechngdist[3,2:tN] <- DTseam[changer==T & switchedOcc_wave==T,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
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
names(tab_wavechngdist_out) <- c("Mean",tabqtls)
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_wavechngdist_out) <- c("Chng\ Job, All", "Chng\ Job, Same\ Occ", "Chng\ Job & Switch\ Occ")

tab_wavechngdist_out <- xtable(tab_wavechngdist_out, label=paste0("tab:",labtxt), digits=2, 
                           align="l|l|lllll", caption="Distribution of earnings changes among job changers")
if(demolbl>=1 & demolbl<=7){
	print(tab_wavechngdist_out,include.rownames=T, hline.after= c(0,nrow(tab_wavechngdist_out)), file=paste0(outputdir,"/",labtxt,demotxt[demolbl],".tex"))
}else{
	print(tab_wavechngdist_out,include.rownames=T, hline.after= c(0,nrow(tab_wavechngdist_out)), file=paste0(outputdir,"/",labtxt,".tex"))
}
# rec and expansion
tab_wavechngdist_rec    <- array(NA,dim=c(6,tN))
tab_wavechngdist_rec_EEEU    <- array(NA,dim=c(12,tN))
for(recHere in c(F,T)){
	ridx = as.integer(recHere)*3

	tab_wavechngdist_rec[1+ridx,1   ] <- DTseam[changer==T & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[1+ridx,2:tN] <- DTseam[changer==T & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs= tabqtls)]
	tab_wavechngdist_rec[2+ridx,1   ] <- DTseam[changer==T &switchedOcc_wave==F & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[2+ridx,2:tN] <- DTseam[changer==T &switchedOcc_wave==F & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	tab_wavechngdist_rec[3+ridx,1   ] <- DTseam[changer==T &switchedOcc_wave==T & recIndic == recHere,    wtd.mean(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)))]
	tab_wavechngdist_rec[3+ridx,2:tN] <- DTseam[changer==T &switchedOcc_wave==T & recIndic == recHere,wtd.quantile(eval(as.name(wc)),na.rm=T,weights=eval(as.name(wt)), probs=tabqtls)]
	
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
cnames <- c("Mean",tabqtls)
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
                           align="l|l|lllll", caption=paste0("Distribution of earnings changes among job changers in recession and expansion \\label{tab:",labtxt,"_rec}"))
print(tab_wavechngdist_rec,include.rownames=T, hline.after= c(nrow(tab_wavechngdist_rec)), 
      add.to.row=rowtitles, file=paste0(outputdir,"/",labtxt,"_rec.tex"))



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

