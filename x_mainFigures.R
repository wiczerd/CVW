# x_mainFigures.R


library(data.table)
library(zoo)
library(stats)
library(Hmisc)
library(ggplot2)
library(xtable)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outdir = "~/workspace/CVW/R/Figures"
setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}
Max_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		max(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA_real_
		}
	}
}
Min_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		min(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA_real_
		}
	}
}
Any_narm <- function(x) {
	ux <- unique(x[!is.na(x)])
	if(length(ux)>0){
		any(ux)
	}else{
		if(is.integer(x)){
			NA_integer_
		}else{
			NA
		}
	}
}
 
stackedDist <- function( DT,gg,wc,pwc ){
	Ng <- DT[ is.finite(eval(as.name(gg)))==T , length(unique( eval(as.name(gg)) ))]
	stackedCcount <- array(0.,dim=c(1001*Ng,3))
	stackedCdist <- array(0.,dim=c(1001*Ng,3))
	for(wi in seq(0,1000) ){
		stackedCdist[wi*Ng+seq(1,Ng),1] <- seq(1,Ng)
		stackedCcount[wi*Ng+seq(1,Ng),1] <- seq(1,Ng)
		meanwagebin <-DT[ eval(as.name(pwc)) == wi, mean(eval(as.name(wc)))]
		stackedCdist[wi*Ng+seq(1,Ng),2] <- meanwagebin
		stackedCcount[wi*Ng+seq(1,Ng),2] <- meanwagebin
		tmp<-DT[ eval(as.name(pwc))== wi & is.finite( eval(as.name(gg)) ), sum(is.finite(eval(as.name(pwc)))),by = eval(as.name(gg))]
		names(tmp) <- c("gg","V1")
		setkey(tmp,gg)
		stackedCcount[as.integer(tmp$gg)+wi*Ng,3] <- tmp$V1
		stackedCdist[as.integer(tmp$gg)+wi*Ng,3] <- tmp$V1/sum(tmp$V1)
	}
	stackedCdist <- as.data.table(stackedCdist)
	names(stackedCdist) <- c("g","WageChange","Pct")
	return( stackedCdist )
}

stackedDens <- function( DT,gg,wc,pwc ){
	Ng <- DT[ is.finite(eval(as.name(gg)))==T , length(unique( eval(as.name(gg)) ))]
	stackedCcount <- array(0.,dim=c(1001*Ng,3))
	stackedCdist <- array(0.,dim=c(1001*Ng,3))
	for(wi in seq(0,1000) ){
		stackedCdist[wi*Ng+seq(1,Ng),1] <- seq(1,Ng)
		stackedCcount[wi*Ng+seq(1,Ng),1] <- seq(1,Ng)
		meanwagebin <-DT[ eval(as.name(pwc)) == wi, mean(eval(as.name(wc)))]
		stackedCdist[wi*Ng+seq(1,Ng),2] <- meanwagebin
		stackedCcount[wi*Ng+seq(1,Ng),2] <- meanwagebin
		tmp<-DT[ eval(as.name(pwc))== wi & is.finite( eval(as.name(gg)) ), sum(is.finite(eval(as.name(pwc)))),by = eval(as.name(gg))]
		names(tmp) <- c("gg","V1")
		setkey(tmp,gg)
		stackedCcount[as.integer(tmp$gg)+wi*Ng,3] <- tmp$V1
		stackedCdist[as.integer(tmp$gg)+wi*Ng,3] <- tmp$V1/sum(tmp$V1)
	}
	stackedCcount <- as.data.table(stackedCcount)
	names(stackedCcount) <- c("g","WageChange","Cnt")
	stackedCcount[ , incr:= WageChange - shift(WageChange),by=g]
	stackedCcount[ , incr:= na.locf(incr,fromLast=T)]
	stackedCcount[ , Pct:= Cnt/sum(Cnt)*incr]
	return( stackedCcount )
}

DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))


DTseam <- subset(DTseam, changer==T|stayer==T)

toKeep <- c("truncweight","cycweight","wpfinwgt","EU_wave","UE_wave","EE_wave","switchedOcc_wave","wagechange_wave","wagechangeEUE_wave","recIndic_wave","date")
DTseam <- subset(DTseam, is.finite(wagechange_wave) & is.finite(EU_wave) & is.finite(UE_wave)& is.finite(EE_wave))

# Occupation switchers and not--------------

DTseam[  DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F , g := 1]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==T & DTseam$EU_wave ==F , g := 2]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T , g := 3]
DTseam[!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T), g := 4]
DTseam[!is.na(DTseam$g), g1 := ifelse(g==1,1,0)]
DTseam[!is.na(DTseam$g), g2 := ifelse(g==2,1,0)]
DTseam[!is.na(DTseam$g), g3 := ifelse(g==3,1,0)]
DTseam[!is.na(DTseam$g), g4 := ifelse(g==4,1,0)]

mid99 <-DTseam[ , quantile(wagechange_wave,probs = c(.005,.995))]
DTseam[ wagechange_wave>mid99[1] & wagechange_wave<mid99[2], rnk_wagechange_wave  := frank(wagechange_wave)]
DTseam[ , pct1000_wagechange_wave  := round(1000*rnk_wagechange_wave/max(rnk_wagechange_wave,na.rm=T),digits=0)]
mid99 <- DTseam[ , quantile(rawwgchange_wave,probs = c(.005,.995),na.rm = T)]
DTseam[ rawwgchange_wave>mid99[1] & rawwgchange_wave<mid99[2], rnk_rawwgchange_wave := frank(rawwgchange_wave)]
DTseam[ , pct1000_rawwgchange_wave := as.integer(round(1000*rnk_rawwgchange_wave/max(rnk_rawwgchange_wave,na.rm=T),digits=0))]

stackedCdist <- stackedDist(DTseam,"g","wagechange_wave","pct1000_wagechange_wave")
stackedlogDens <- stackedDens(DTseam,"g","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
  scale_fill_manual( values = c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]) ,
                    labels=c("EE","UE","EU","Stay"))+ xlab("Residual Earnings Change")+ylab("Fraction by Type")+
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_wagechange.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange.png"),height=5,width=10)

stackedCdist <- stackedDist(DTseam,"g","rawwgchange_wave","pct1000_rawwgchange_wave")

ggplot(stackedCdist, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
  scale_fill_manual( values = c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]) ,
                      labels=c("EE","UE","EU","Stay"))+xlab("Nominal Earnings Change")+ylab("Fraction by Type")+
  theme(legend.title = element_blank(),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(linetype = "solid",color = "white"))

ggsave(paste0(outdir,"/stacked_rawwgchange.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_rawwgchange.png"),height=5,width=10)


# recession or not recession:
stackedCdist_rec <- stackedDist(subset(DTseam,recIndic_wave==T),"g","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist_rec, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
	scale_fill_manual( values = c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]) ,
					   labels=c("EE","UE","EU","Stay"))+ xlab("Residual Earnings Change, Recession")+ylab("Fraction by Type")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_wagechange_rec.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange_rec.png"),height=5,width=10)

stackedCdist_exp <- stackedDist(subset(DTseam, recIndic_wave==F),"g","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist_exp, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
	scale_fill_manual( values = c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]) ,
					   labels=c("EE","UE","EU","Stay"))+ xlab("Residual Earnings Change, Expansion")+ylab("Fraction by Type")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_wagechange_exp.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange_exp.png"),height=5,width=10)


# Occupation switchers and not--------------

DTseam[DTseam$switchedOcc_wave==T & DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F , s := 1]
DTseam[DTseam$switchedOcc_wave==T & DTseam$EE_wave==F & DTseam$UE_wave==T & DTseam$EU_wave ==F , s := 3]
DTseam[DTseam$switchedOcc_wave==T & DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T , s := 5]
DTseam[DTseam$switchedOcc_wave==F & DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F , s := 2]
DTseam[DTseam$switchedOcc_wave==F & DTseam$EE_wave==F & DTseam$UE_wave==T & DTseam$EU_wave ==F , s := 4]
DTseam[DTseam$switchedOcc_wave==F & DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T , s := 6]
DTseam[DTseam$switchedOcc_wave==T &!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T), s := 7]
DTseam[DTseam$switchedOcc_wave==F &!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T), s := 8]

stackedCdist <- stackedDist(DTseam,"s","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
	scale_fill_manual( values = c(hcl(h=seq(15, 375, length=9), l=50, c=100)[c(1:8)]) ,
					   labels=c("EE, Sw","EE, no Sw","UE, Sw","UE, no Sw","EU, Sw","EU, no Sw","Stay Sw", "Stay, no Sw"))+ 
	xlab("Residual Earnings Change")+ylab("Fraction by Type")+
	theme(legend.title = element_blank(),
		  #legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_occsw_wagechange.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_occsw_wagechange.png"),height=5,width=10)

stackedCdist_rec <- stackedDist(subset(DTseam,recIndic_wave==T),"s","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist_rec, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
	scale_fill_manual( values = c(hcl(h=seq(15, 375, length=9), l=50, c=100)[c(1:8)]) ,
					   labels=c("EE, Sw","EE, no Sw","UE, Sw","UE, no Sw","EU, Sw","EU, no Sw","Stay Sw", "Stay, no Sw"))+ 
	xlab("Residual Earnings Change, Recession")+ylab("Fraction by Type")+
	theme(legend.title = element_blank(),
		  #legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_occsw_wagechange_rec.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_occsw_wagechange_rec.png"),height=5,width=10)

stackedCdist_exp <- stackedDist(subset(DTseam,recIndic_wave==F),"s","wagechange_wave","pct1000_wagechange_wave")

ggplot(stackedCdist_exp, aes(x=WageChange,y=Pct, fill = as.factor(g))) + geom_area() + theme_bw()+
	scale_fill_manual( values = c(hcl(h=seq(15, 375, length=9), l=50, c=100)[c(1:8)]) ,
					   labels=c("EE, Sw","EE, no Sw","UE, Sw","UE, no Sw","EU, Sw","EU, no Sw","Stay Sw", "Stay, no Sw"))+ 
	xlab("Residual Earnings Change, Expansion")+ylab("Fraction by Type")+
	theme(legend.title = element_blank(),
		  #legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/stacked_occsw_wagechange_exp.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_occsw_wagechange_exp.png"),height=5,width=10)


# Figures with monthly data ----------------------------------

DTall <- readRDS( paste0(datadir,"/DTall_5.RData"))

#drop some un-need variables

#DTall[ , nwhite := black==T|hispanic==T]
#DTall <- DTall[ , c("eyear","emonth","esr","next.earnm","logearnm","black","hispanic","EE_ann","EU_ann","UE_ann",
#					"occ90","occwage","midUE_ann","midEE_ann","midEU_ann","estlemp","jobchng","jobchng_ann"):=NULL]

setkey(DTall, id, date)


#Transitions ------------------------

#set up EN, NE to compar to EU, UE
DTall[lfstat == 1, EN := (next.lfstat == 3)]
DTall[lfstat == 3, NE := (next.lfstat == 1)]
DTall[ , next.ustintid:= shift(ustintid,type = "lead"),by=id]
DTall[ EN==T, ustintid:= next.ustintid]

DTall[ mis==1, ltrunc := (lfstat>1)]
DTall[       , ltrunc := any(ltrunc,na.rm = F), by=list(id,ustintid)]
DTall[is.na(ltrunc), ltrunc := F]
DTall[ ustintid>0 & ltrunc==F, nempdur := seq_len(.N)-1, by=list(id,ustintid)]
DTall[ ustintid>0 & ltrunc==F, max.nempdur := max(nempdur,na.rm = T), by=list(id,ustintid)]
DTall[ , max.nempdur_wave := Max_narm(max.nempdur), by=list(id,wave)]


#Single variable for the type of transition ------
DTseam <- readRDS(paste0(datadir,"/DTall_5.RData"))

#Set up wage change for NE-------------------------
DTall[ , NE_max := any(NE,na.rm=T), by=list(id,wave)]
# DTall[ , EN_max := any(EN,na.rm=T), by=list(id,wave)]
#  
# DTseam <- DTall[ seam==T,]
# DTseam[ , next.wavewage := shift(wavewage,1,type="lead"), by=id]
# DTseam[ , last.wavewage := shift(wavewage,1,type="lag"), by=id]
# DTseam[ , next.ustintid_wave := shift(ustintid_wave,1,type="lead"), by=id]
# DTseam[ EN_max==T, ustintid_wave:= next.ustintid_wave]
# DTseam[ EN_max==T & is.na( ustintid_wave) & next.lfstat_wave==1, ustintid_wave := 10L]
# 
# 
# DTseam[last.lfstat_wave==1 & EN_max == T, wageAtEN := last.wavewage]
# DTseam[, wageAtEN := na.locf(wageAtEN, na.rm = F),by=list(ustintid_wave, id)]
# DTseam[next.lfstat_wave==1 & NE_max == T, wageAfterNE :=  next.wavewage]
# DTseam[NE_max == T, wagechangeENE_wave := wageAfterNE - wageAtEN]
# DTseam[, wagechangeENE_wave:= Mode(wagechangeENE_wave), by=list(ustintid_wave, id)]
# DTseam[ EE_wave==T, wagechangeENE_wave := wagechange_wave]
# DTseam[ !(EN_max |NE_max), wagechangeENE_wave := wagechange_wave]
# 
# DTseam <- subset(DTseam, select=c("id","wave","wagechangeENE_wave"))
# 
# DTall <- merge(DTall,DTseam, by=c("id","wave"), all=T)
# rm(DTseam)

#Transition rates and wages---------------------------

frateUE <- lm(formula = UE ~ nempdur + I(nempdur^2) + I(nempdur^3), data = subset(DTall, lfstat == 2))
# frateNE <- lm(formula = NE ~ nempdur + I(nempdur^2) + I(nempdur^3), data = subset(DTall, lfstat == 3))

frateDat <- array(NA, dim=c(24,2))
frateDat[,1] <- seq_len(24)
for( t in seq(1,24)){
	frateDat[t,2] <- DTall[ nempdur==t+1 & lfstat==2, mean(UE,na.rm = T)]
	# frateDat[t,3] <- DTall[ nempdur==t+1 & lfstat==3, mean(NE,na.rm = T)]
}
frateFit <- array(NA, dim=c(24,2))
frateFit[,1] <- seq_len(24)
frateFit[,2] <- frateFit[,1]*frateUE$coefficients[2] +frateFit[,1]^2*frateUE$coefficients[3] + frateFit[,1]^3*frateUE$coefficients[4] + frateUE$coefficients[1]
# frateFit[,3] <- frateFit[,1]*frateNE$coefficients[2] +frateFit[,1]^2*frateNE$coefficients[3] + frateFit[,1]^3*frateNE$coefficients[4] + frateNE$coefficients[1]

frateNN <- DTall[ lfstat==3 & ltrunc==T, mean(NE,na.rm = T)]

# plot these:
frateFitPlot <- data.table(frateFit)
names(frateFitPlot) <- c("Duration","EUE") #,"ENE"
frateDatPlot <- data.table(frateDat)
names(frateDatPlot) <- c("Duration","EUE")
#frateFitPlot[ , NNE := as.numeric(frateNN) ]
frateFitMelt <- melt(frateFitPlot, id.vars = "Duration")
frateDatMelt <- melt(frateDatPlot, id.vars = "Duration")

ggplot(data=frateFitMelt, aes(x=Duration,y=value)) + geom_line(size=1.5)+
	theme_bw()+xlab("Duration (mo)")+ylab("Finding Rate")+ylim(c(0,.23))+
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3,1,2)]) ,
						labels=c("Thru U","Thru N","Never E","",""))+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))+
	geom_point(data=frateDatMelt, aes(x=Duration,y=value))

ggsave(paste0(outdir,"/frateFit.eps"),height=5,width=10)
ggsave(paste0(outdir,"/frateFit.png"),height=5,width=10)

# ggplot(data=subset(frateFitMelt, variable%in%c("EUE","NNE")), aes(x=Duration,y=value,color=variable)) + geom_line(size=1.5)+
# 	theme_bw()+xlab("Duration (mo)")+ylab("Finding Rate")+ylim(c(0,.23))+
# 	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1,3)]) ,
# 						labels=c("Thru U","Never E"))+
# 	theme(legend.title = element_blank(),
# 		  legend.position = c(0.8,0.8),
# 		  legend.background = element_rect(linetype = "solid",color = "white"))+
# 	geom_point(data=subset(frateDatMelt, variable%in%c("EUE")), aes(x=Duration,y=value,color=variable))	
# 
# ggsave(paste0(outdir,"/frateFit_EUE_NNE.eps"),height=5,width=10)
# ggsave(paste0(outdir,"/frateFit_EUE_NNE.png"),height=5,width=10)

wchngUE<- lm( wagechangeEUE_wave ~max.nempdur_wave + I(max.nempdur_wave^2)+ I(max.nempdur_wave^3), data=subset(DTall, UE_wave==T & midUE==F & seam==T) )

wchngUE_Fit<- array(NA,dim=c(24,2))
wchngUE_Dat<- array(NA,dim=c(24,2))
wchngUE_Dat[ , 1]<- seq_len(24)
wchngUE_Fit[ , 1]<- seq_len(24)
for(t in seq(1,24)){
	wchngUE_Dat[t,2] <- DTall[ max.nempdur_wave==t & UE_wave==T & midUE==F & seam==T, mean(wagechangeEUE_wave,na.rm = T)]
}
wchngUE_Fit[,2] = wchngUE_Fit[,1]*wchngUE$coefficients[2] +wchngUE_Fit[,1]^2*wchngUE$coefficients[3] + wchngUE_Fit[,1]^3*wchngUE$coefficients[4] + wchngUE$coefficients[1]


# wchngNE<- lm( wagechangeENE_wave ~max.nempdur_wave + I(max.nempdur_wave^2)+ I(max.nempdur_wave^3), data=subset(DTall, NE_max==T & seam==T ) )

# wchngNE_Fit<- array(NA,dim=c(24,2))
# wchngNE_Dat<- array(NA,dim=c(24,2))
# wchngNE_Dat[ , 1]<- seq_len(24)
# wchngNE_Fit[ , 1]<- seq_len(24)
# for(t in seq(1,24)){
# 	wchngNE_Dat[t,2] <- DTall[ max.nempdur_wave==t & seam==T & NE_max==T , mean(wagechangeENE_wave,na.rm = T)]
# }
# wchngNE_Fit[,2] = wchngNE_Fit[,1]*wchngNE$coefficients[2] +wchngNE_Fit[,1]^2*wchngNE$coefficients[3] + wchngNE_Fit[,1]^3*wchngNE$coefficients[4] + wchngNE$coefficients[1]

# plot these:
wchngFitPlot <- data.table(wchngUE_Fit)
names(wchngFitPlot) <- c("Duration","EUE")
EEchng <- DTall[ EE_wave==T, mean(wagechange_wave,na.rm = T)]
wchngFitPlot[ , EE:= EEchng]
wchngDatPlot <- data.table(wchngUE_Dat)
names(wchngDatPlot) <- c("Duration","EUE")
wchngDatPlot[ , EE:= EEchng]
wchngFitMelt <- melt(wchngFitPlot, id.vars = "Duration")
wchngDatMelt <- melt(wchngDatPlot, id.vars = "Duration")

ggplot(data=wchngFitMelt, aes(x=Duration,y=value,color=variable)) + geom_line(size=1.5)+
	theme_bw()+xlab("Completed Duration (mo)")+ylab("Earnings Change (pct)")+ylim(c(-.4,.15))+
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1,2,1,2)]) ,
						labels=c("EUE","EE","",""))+
	theme(legend.title = element_blank(),
		  legend.position = c(0.2,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))+
	geom_point(data=wchngDatMelt, aes(x=Duration,y=value,color=variable))

ggsave(paste0(outdir,"/wchngFit.eps"),height=5,width=10)
ggsave(paste0(outdir,"/wchngFit.png"),height=5,width=10)



DTall[ seam==T & last.lfstat_wave==1, last.wavewage := shift(wavewage), by=id]

whois <- array(NA, dim = c(2,6))
rownames(whois)<-c("EN","EU")
colnames(whois)<-c("Mean","25","50","75","Female","<=30")
wmean <- DTall[ seam==T & lfstat_wave==1 & last.wavewage>0, wtd.mean(last.wavewage , weights= wpfinwgt, na.rm = T)]
whois[1,1]   <- DTall[ seam==T & NE_max==T           & last.wavewage>0, wtd.mean(last.wavewage     , weights= wpfinwgt, na.rm = T)] - wmean
whois[1,2:4] <- DTall[ seam==T & NE_max==T           & last.wavewage>0, wtd.quantile(last.wavewage , weights= wpfinwgt, probs = c(.25,.5,.75),na.rm = T)] - wmean
whois[2,1]   <- DTall[ seam==T & EU_wave==T&midEU==F & last.wavewage>0, wtd.mean(last.wavewage     , weights= wpfinwgt, na.rm = T)] -wmean
whois[2,2:4] <- DTall[ seam==T & EU_wave==T&midEU==F & last.wavewage>0, wtd.quantile(last.wavewage , weights= wpfinwgt, probs = c(.25,.5,.75),na.rm = T)] -wmean

whois[1,5]   <- DTall[ seam==T & NE_max==T           , wtd.mean(female, weights= wpfinwgt, na.rm = T)]
whois[2,5]   <- DTall[ seam==T & EU_wave==T&midEU==F , wtd.mean(female, weights= wpfinwgt, na.rm = T)]
whois[1,6]   <- DTall[ seam==T & NE_max==T           , wtd.mean(Young, weights= wpfinwgt, na.rm = T)]
whois[2,6]   <- DTall[ seam==T & EU_wave==T&midEU==F , wtd.mean(Young, weights= wpfinwgt, na.rm = T)]


outputtable <- xtable(whois, digits=2, 
					  align="l|cccc|ll", caption=paste0("Earnings, Relative to Mean, Prior to Separating"))
print(outputtable,include.rownames=T, hline.after= c(0,0,1,nrow(outputtable)), file=paste0(outdir,"/worker_chars.tex"))

