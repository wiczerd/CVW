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


wcSkew_ts <- DTseam[ truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) , by=list(month,year)]
unrate_ts <- DTseam[ truncweight>0 , wtd.mean(unrt, weights =truncweight)/100, by=list(month,year)]
funrate_ts <- DTseam[ truncweight>0 , wtd.mean(filunrt, weights =truncweight), by=list(month,year)]
wts_ts <- DTseam[ truncweight>0 , sum(truncweight), by=list(month,year)]
wcSkewUn_ts <- merge(wcSkew_ts,unrate_ts, by=c("month", "year"))
wcSkewUn_ts <- merge(wcSkewUn_ts,funrate_ts, by=c("month", "year"))
names(wcSkewUn_ts) <- c("month","year","GMskewness","unrate","HPunrate")
wcSkewUn_ts <- merge(wcSkewUn_ts,wts_ts, by=c("month", "year"))
names(wcSkewUn_ts) <- c("month","year","GMskewness","unrate","HPunrate","truncweight")

summary(lm(GMskewness~HPunrate, weights = truncweight , data=wcSkewUn_ts))


wcSkew_mvts <- DTseam[ ch==T & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) , by=list(month,year)]
unrate_mvts <- DTseam[ ch==T & truncweight>0 , wtd.mean(unrt, weights =truncweight)/100, by=list(month,year)]
funrate_mvts <- DTseam[ ch==T & truncweight>0 , wtd.mean(filunrt, weights =truncweight), by=list(month,year)]
wts_mvts <- DTseam[ ch==T & truncweight>0 , sum(truncweight), by=list(month,year)]
wcSkewUn_mvts <- merge(wcSkew_mvts,unrate_mvts, by=c("month", "year"))
wcSkewUn_mvts <- merge(wcSkewUn_mvts,funrate_mvts, by=c("month", "year"))
names(wcSkewUn_mvts) <- c("month","year","GMskewness","unrate","HPunrate")
wcSkewUn_mvts <- merge(wcSkewUn_mvts,wts_mvts, by=c("month", "year"))
names(wcSkewUn_mvts) <- c("month","year","GMskewness","unrate","HPunrate","truncweight")

summary(lm(GMskewness~HPunrate, weights = truncweight , data=wcSkewUn_mvts))

wcSkew_swts <- DTseam[ sw==T & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) , by=list(month,year)]
unrate_swts <- DTseam[ sw==T & truncweight>0 , wtd.mean(unrt, weights =truncweight)/100, by=list(month,year)]
funrate_swts <- DTseam[ sw==T & truncweight>0 , wtd.mean(filunrt, weights =truncweight), by=list(month,year)]
wts_swts <- DTseam[ sw==T & truncweight>0 , sum(truncweight), by=list(month,year)]
wcSkewUn_swts <- merge(wcSkew_swts,unrate_swts, by=c("month", "year"))
wcSkewUn_swts <- merge(wcSkewUn_swts,funrate_swts, by=c("month", "year"))
names(wcSkewUn_swts) <- c("month","year","GMskewness","unrate","HPunrate")
wcSkewUn_swts <- merge(wcSkewUn_swts,wts_swts, by=c("month", "year"))
names(wcSkewUn_swts) <- c("month","year","GMskewness","unrate","HPunrate","truncweight")

summary(lm(GMskewness~HPunrate, weights = truncweight , data=wcSkewUn_swts))


wcSkew_stts <- DTseam[ st==T & sw==F & truncweight>0, wtd.GroenveldMeeden(wagechange_anan,wt =truncweight) , by=list(month,year)]
unrate_stts <- DTseam[ st==T & sw==F & truncweight>0 , wtd.mean(unrt, weights =truncweight)/100, by=list(month,year)]
funrate_stts <- DTseam[ st==T  & sw==F& truncweight>0 , wtd.mean(filunrt, weights =truncweight), by=list(month,year)]
wts_stts <- DTseam[ st==T  & sw==F& truncweight>0 , sum(truncweight), by=list(month,year)]
wcSkewUn_stts <- merge(wcSkew_stts,unrate_stts, by=c("month", "year"))
wcSkewUn_stts <- merge(wcSkewUn_stts,funrate_stts, by=c("month", "year"))
names(wcSkewUn_stts) <- c("month","year","GMskewness","unrate","HPunrate")
wcSkewUn_stts <- merge(wcSkewUn_stts,wts_stts, by=c("month", "year"))
names(wcSkewUn_stts) <- c("month","year","GMskewness","unrate","HPunrate","truncweight")
wcSkewUn_stts[GMskewness==-Inf, GMskewness:=NaN]

summary(lm(GMskewness~HPunrate, weights = truncweight , data=wcSkewUn_stts))