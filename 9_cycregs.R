# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(xtable)
library(Hmisc)
library(quantreg)
library(ggplot2)
library(mFilter)

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

wtd.Kelley <- function(xt, wt, samp=F){  
  wt  <- wt[is.na(xt)==F]
  xt  <- xt[is.na(xt)==F]
  if(samp==T){
    sampidx <- sample(length(xt), length(xt)/4,replace = T)
    xt <- xt[sampidx]
    wt <- wt[sampidx]
  }
  qtls<-wtd.quantile(xt,weights=wt,probs = c(.10,.50,.90),na.rm = T) 
  ((qtls[3]- qtls[2]) - (qtls[2]- qtls[1]))/(qtls[3]- qtls[1])
}

wtd.9010 <- function(xt, wt, samp=F){  
  wt  <- wt[is.na(xt)==F]
  xt  <- xt[is.na(xt)==F]
  if(samp==T){
    sampidx <- sample(length(xt), length(xt)/4,replace = T)
    xt <- xt[sampidx]
    wt <- wt[sampidx]
  }
  qtls<-wtd.quantile(xt,weights=wt,probs = c(.10,.90),na.rm = T) 
  (qtls[2] - qtls[1])
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



DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))
#create some vars to handle net flows
DTseam[ , last.occ :=shift(occ), by=id]
DTseam[ !is.finite(occL), occL := occ]
DTseam[ EE_wave==T , occL := last.occ]
DTseam[ , next.occ :=shift(occ,type="lead"), by=id]
DTseam[ !is.finite(occD), occD := next.occ]



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

#GMskewregs_mvsw = vector(mode="list",length=(4))
#MadDispregs_mvsw = vector(mode="list",length=(4))
#Kskewregs_mvsw = vector(mode="list",length=(4))
#DecDispregs_mvsw = vector(mode="list",length=(4))
GMskewregs_mvsw = array(data=NA, dim=c(4,2))
MadDispregs_mvsw = array(data=NA, dim=c(4,2))
Kskewregs_mvsw = array(data=NA, dim=c(4,2))
DecDispregs_mvsw = array(data=NA, dim=c(4,2))


# Move / no move, Switch / no switch
for( chI in c(T,F)  ){
  for( swI in c(T,F)){


    #the sensitivity of GM to recession
    GMskewregs_mvsw[chI*2 + swI + 1,1]  = DTseam[ recIndic2_wave==F & ch==chI  & sw==swI & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) ]
    GMskewregs_mvsw[chI*2 + swI + 1,2]  = DTseam[ recIndic2_wave==T & ch==chI  & sw==swI & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) ]
    
    MadDispregs_mvsw[chI*2 + swI + 1,1] = DTseam[ recIndic2_wave==F & ch==chI  & sw==swI & truncweight>0, wtd.mad(wagechange_anan,wt =truncweight) ]
    MadDispregs_mvsw[chI*2 + swI + 1,2] = DTseam[ recIndic2_wave==T & ch==chI  & sw==swI & truncweight>0, wtd.mad(wagechange_anan,wt =truncweight) ]

    Kskewregs_mvsw[chI *2 + swI + 1,1] = DTseam[ recIndic2_wave==F & ch==chI  & sw==swI & truncweight>0, wtd.Kelley(wagechange_anan,wt =truncweight) ]
    Kskewregs_mvsw[chI *2 + swI + 1,2] = DTseam[ recIndic2_wave==T & ch==chI  & sw==swI & truncweight>0, wtd.Kelley(wagechange_anan,wt =truncweight) ]

    DecDispregs_mvsw[chI *2 + swI + 1,1] = DTseam[ recIndic2_wave==F & ch==chI  & sw==swI & truncweight>0, wtd.9010(wagechange_anan,wt =truncweight) ]
    DecDispregs_mvsw[chI *2 + swI + 1,2] = DTseam[ recIndic2_wave==T & ch==chI  & sw==swI & truncweight>0, wtd.9010(wagechange_anan,wt =truncweight) ]
    
  }
}
rownames(GMskewregs_mvsw)  = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(Kskewregs_mvsw)   = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(MadDispregs_mvsw) = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(DecDispregs_mvsw) = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")

(GMskewregs_mvsw[,1] - GMskewregs_mvsw[,2])/abs(GMskewregs_mvsw[,1])
(Kskewregs_mvsw[,1] - Kskewregs_mvsw[,2])/abs(Kskewregs_mvsw[,1])

(MadDispregs_mvsw[,1] - MadDispregs_mvsw[,2])/abs(MadDispregs_mvsw[,1])
(DecDispregs_mvsw[,1] - DecDispregs_mvsw[,2])/abs(DecDispregs_mvsw[,1])


#++++++++++++++++++++++++++++++++++
# Now do it with the model output:
#++++++++++++++++++++++++++++++++++

outdir <- "~/Dropbox/Carrillo_Visschers_Wiczer/Computation/SimulationFigures/"
readdir <- "~/workspace/CVW/c/exper"



simDat <- readRDS(paste0(readdir,"/simDat.RData"))

setkeyv(simDat, c("id","t"))
simDat[ , last.unemp:= shift(unemp), by=id]
simDat[ , last.occ  := shift(occ), by=id]
simDat[ , next.occ  := shift(occ,type="lead"), by=id]
simDat[ , last.trx  := shift(trx), by=id]

simDat[ wage <= - 99, wage:= NA]
#simDat[ , last.wage := shift(wage), by=id]
#simDat[ , last.wage_an:= shift(wage)+shift(wage,2)+shift(wage,3), by=id]
#simDat[ , next.wage_an:=  wage +shift(wage,1,type="lead")+shift(wage,2,type="lead"), by=id]
#simDat[ unemp==0 & last.unemp==0, wagechng := log(wage) - log(last.wage), by=id]
#simDat[ next.wage_an>0 & last.wage_an>0, wagechng_anan := log(next.wage_an) - log(last.wage_an), by=id]

simDat[ , modelyr := ((t-1) - (t-1)%%12)/12]

simDat[ , occ_chng := shift(occ,type="lead") != occ, by=id]
simDat[ , eps_chng := epsilon != shift(epsilon), by=id]
simDat[ , EE       := shift(J2J,type="lead")]
simDat[ , EU       := unemp==0 & shift(unemp,type="lead")==1]
simDat[ , UE       := unemp==1 & shift(unemp,type="lead")==0]
simDat[ unemp==0 | UE==F, occ_chng := NA ]
simDat[ unemp==0, occL:=as.integer(-1)]
simDat[ unemp==0, occD:=as.integer(-1)]
simDat[ EU==T, occL := last.occ]
simDat[ EU==T, occD := NA]
simDat[ UE==T, occD := next.occ]
simDat[ , occL := na.locf0(occL), by=id]
simDat[ , occD := na.locf(occD,fromLast = T,na.rm = F), by=id]
simDat[ UE==T, occ_chngUE:= occL!=occD]
simDat[ EU==T, occ_chngUE:= occL!=occD]


simDat[ , smax := pmax(s0,s1,s2,s3)]
simDat[ smax>0, sdirected:= (smax-.33333)/ (1-.3333)]
simDat[ s0>0|s1>0|s2>0 | s3>0, smin := pmin(s0,s1,s2,s3)]
simDat[ is.finite(smin) & smax>0, smindirected:= (smin-.33333)/ (1-.3333)]

#HP filter unemployment in the model:
unempTS <- simDat[ , mean(unemp), by=t]

#make residual wages
simDat[ wage>0 & unemp==0, rwage:= lm(log(wage)~  as.factor(occ))$residuals]
simDat[ unemp==1 | wage==0, rwage:= -100]
simDat[ is.finite(wchng), lwave_rwage:= log((exp(rwage)+exp(shift(rwage))+ exp(shift(rwage,2))+ exp(shift(rwage,3)))/4), by =id]
simDat[ is.finite(wchng), lwave_wage:= log((wage+(shift(wage))+ (shift(wage,2))+ (shift(wage,3)))/4), by =id]
simDat[ rwage<=-99 & shift(rwage)<=-99 & shift(rwage,2)<=-99 & shift(rwage,3)<=-99, lwave_rwage:=NA]
simDat[ rwage<=-99 & shift(rwage)<=-99 & shift(rwage,2)<=-99 & shift(rwage,3)<=-99, lwave_wage:=NA]
simDat[ , nextw.occ_chng := trx==5|trx==6|trx==4|
          shift(trx,type="lead")==5|shift(trx,type="lead")==6|shift(trx,type="lead")==4|
          shift(trx,type="lead",2)==5|shift(trx,type="lead",2)==6|shift(trx,type="lead",2)==4 , by=id]

simDat[ , next.wage_an:=
          wage                      + shift(wage,type="lead"  )  + shift(wage,type="lead",2) + 
          shift(wage,type="lead",3) + shift(wage,type="lead",4)  + shift(wage,type="lead",5) + 
          shift(wage,type="lead",6) + shift(wage,type="lead",7)  + shift(wage,type="lead",8) + 
          shift(wage,type="lead",9) + shift(wage,type="lead",10) + shift(wage,type="lead",11) , by=id]

simDat[ , last.wage_an:=      
          shift(wage)               + shift(wage,type="lag",2)  + shift(wage,type="lag",3) + 
          shift(wage,type="lag",4)  + shift(wage,type="lag",5)  + shift(wage,type="lag",6) + 
          shift(wage,type="lag",7)  + shift(wage,type="lag",8)  + shift(wage,type="lag",9) + 
          shift(wage,type="lag",10) + shift(wage,type="lag",11) + shift(wage,type="lag",12) , by=id]

simDat[ , wagechange_anan := next.wage_an - last.wage_an]

simDat[ , wave:= floor(t/4)]
simDat[ , seam:= wave!=shift(wave,type="lag"),by=id]
simDat[ , EE_wave := max(trx==0|trx==4), by=list(wave,id)]
simDat[ , EU_wave := max(trx==1|trx==5), by=list(wave,id)]
simDat[ , UE_wave := max(trx==2|trx==6), by=list(wave,id)]
simDat[ , sw_wave := max(trx>=4&trx<=7), by=list(wave,id)]
simDat[ , unemp_wave := max(unemp), by=list(wave,id)]
simDat[ wchng>-99, wchng_wave := shift(wchng,type="lead"), by=list(id)]
simDat[ , lrwage_wave := max(lwave_rwage), by=list(wave,id)]
simDat[ , lwage_wave := max(lwave_wage), by=list(wave,id)]
simSeam <- simDat[ seam==T, ]
simSeam[ , truncweight:=1]

simSeam[ , ch:=  (EE_wave==1 | UE_wave==1 | EU_wave==1)]
simSeam[ , sw:= sw_wave]

GMskewregs_mod = array(data=NA, dim=c(4,2))
MadDispregs_mod = array(data=NA, dim=c(4,2))
Kskewregs_mod = array(data=NA, dim=c(4,2))
DecDispregs_mod = array(data=NA, dim=c(4,2))


# Move / no move, Switch / no switch
for( chI in c(T,F)  ){
  for( swI in c(T,F)){
    
    #the sensitivity of GM to recession
    GMskewregs_mod[chI*2 + swI + 1,1]  = simSeam[ rec==F & ch==chI  & sw==swI & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) ]
    GMskewregs_mod[chI*2 + swI + 1,2]  = simSeam[ rec==T & ch==chI  & sw==swI & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) ]
    
    MadDispregs_mod[chI*2 + swI + 1,1] = simSeam[ rec==F & ch==chI  & sw==swI & truncweight>0, wtd.mad(wagechange_anan,wt =truncweight) ]
    MadDispregs_mod[chI*2 + swI + 1,2] = simSeam[ rec==T & ch==chI  & sw==swI & truncweight>0, wtd.mad(wagechange_anan,wt =truncweight) ]
    
    Kskewregs_mod[chI *2 + swI + 1,1] = simSeam[ rec==F & ch==chI  & sw==swI & truncweight>0, wtd.Kelley(wagechange_anan,wt =truncweight) ]
    Kskewregs_mod[chI *2 + swI + 1,2] = simSeam[ rec==T & ch==chI  & sw==swI & truncweight>0, wtd.Kelley(wagechange_anan,wt =truncweight) ]
    
    DecDispregs_mod[chI *2 + swI + 1,1] = simSeam[ rec==F & ch==chI  & sw==swI & truncweight>0, wtd.9010(wagechange_anan,wt =truncweight) ]
    DecDispregs_mod[chI *2 + swI + 1,2] = simSeam[ rec==T & ch==chI  & sw==swI & truncweight>0, wtd.9010(wagechange_anan,wt =truncweight) ]
    
  }
}

rownames(GMskewregs_mod)  = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(Kskewregs_mod)   = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(MadDispregs_mod) = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")
rownames(DecDispregs_mod) = c("ch0sw0", "ch0sw1", "ch1sw0" ,"ch1sw1")

(GMskewregs_mod[,1] - GMskewregs_mod[,2])/abs(GMskewregs_mod[,1])
(Kskewregs_mod[,1] - Kskewregs_mod[,2])/abs(Kskewregs_mod[,1])

(MadDispregs_mod[,1] - MadDispregs_mod[,2])/abs(MadDispregs_mod[,1])
(DecDispregs_mod[,1] - DecDispregs_mod[,2])/abs(DecDispregs_mod[,1])
