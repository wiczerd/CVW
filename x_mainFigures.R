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


#Set up wage change for NE-------------------------
DTall[ , NE_max := any(NE,na.rm=T), by=list(id,wave)]
DTall[ , EN_max := any(EN,na.rm=T), by=list(id,wave)]

DTseam <- DTall[ seam==T,]
DTseam[ , next.wavewage := shift(wavewage,1,type="lead"), by=id]
DTseam[ , last.wavewage := shift(wavewage,1,type="lag"), by=id]
DTseam[ , next.ustintid_wave := shift(ustintid_wave,1,type="lead"), by=id]
DTseam[ EN_max==T, ustintid_wave:= next.ustintid_wave]
DTseam[ EN_max==T & is.na( ustintid_wave) & next.lfstat_wave==1, ustintid_wave := 10L]


DTseam[last.lfstat_wave==1 & EN_max == T, wageAtEN := last.wavewage]
DTseam[, wageAtEN := na.locf(wageAtEN, na.rm = F),by=list(ustintid_wave, id)]
DTseam[next.lfstat_wave==1 & NE_max == T, wageAfterNE :=  next.wavewage]
DTseam[NE_max == T, wagechangeENE_wave := wageAfterNE - wageAtEN]
DTseam[, wagechangeENE_wave:= Mode(wagechangeENE_wave), by=list(ustintid_wave, id)]
DTseam[ EE_wave==T, wagechangeENE_wave := wagechange_wave]
DTseam[ !(EN_max |NE_max), wagechangeENE_wave := wagechange_wave]

DTseam <- subset(DTseam, select=c("id","wave","wagechangeENE_wave"))

DTall <- merge(DTall,DTseam, by=c("id","wave"), all=T)
rm(DTseam)

#Transition rates and wages---------------------------

frateUE <- lm(formula = UE ~ nempdur + I(nempdur^2) + I(nempdur^3), data = subset(DTall, lfstat == 2))
frateNE <- lm(formula = NE ~ nempdur + I(nempdur^2) + I(nempdur^3), data = subset(DTall, lfstat == 3))

frateDat <- array(NA, dim=c(24,3))
frateDat[,1] <- seq_len(24)
for( t in seq(1,24)){
	frateDat[t,2] <- DTall[ nempdur==t+1 & lfstat==2, mean(UE,na.rm = T)]
	frateDat[t,3] <- DTall[ nempdur==t+1 & lfstat==3, mean(NE,na.rm = T)]
}
frateFit <- array(NA, dim=c(24,3))
frateFit[,1] <- seq_len(24)
frateFit[,2] <- frateFit[,1]*frateUE$coefficients[2] +frateFit[,1]^2*frateUE$coefficients[3] + frateFit[,1]^3*frateUE$coefficients[4] + frateUE$coefficients[1]
frateFit[,3] <- frateFit[,1]*frateNE$coefficients[2] +frateFit[,1]^2*frateNE$coefficients[3] + frateFit[,1]^3*frateNE$coefficients[4] + frateNE$coefficients[1]

frateNN <- DTall[ lfstat==3 & ltrunc==T, mean(NE,na.rm = T)]

# plot these:
frateFitPlot <- data.table(frateFit)
names(frateFitPlot) <- c("Duration","EUE","ENE")
frateDatPlot <- data.table(frateDat)
names(frateDatPlot) <- c("Duration","EUE","ENE")
frateFitPlot[ , NNE := as.numeric(frateNN) ]
frateFitMelt <- melt(frateFitPlot, id.vars = "Duration")
frateDatMelt <- melt(frateDatPlot, id.vars = "Duration")

ggplot(data=frateFitMelt, aes(x=Duration,y=value,color=variable)) + geom_line(size=1.5)+
	theme_bw()+xlab("Duration (mo)")+ylab("Finding Rate")+ylim(c(0,.23))+
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3,1,2)]) ,
						labels=c("Thru U","Thru N","Never E","",""))+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))+
	geom_point(data=frateDatMelt, aes(x=Duration,y=value,color=variable))

ggsave(paste0(outdir,"/frateFit.eps"),height=5,width=10)
ggsave(paste0(outdir,"/frateFit.png"),height=5,width=10)

ggplot(data=subset(frateFitMelt, variable%in%c("EUE","NNE")), aes(x=Duration,y=value,color=variable)) + geom_line(size=1.5)+
	theme_bw()+xlab("Duration (mo)")+ylab("Finding Rate")+ylim(c(0,.23))+
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1,3)]) ,
						labels=c("Thru U","Never E"))+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))+
	geom_point(data=subset(frateDatMelt, variable%in%c("EUE")), aes(x=Duration,y=value,color=variable))	
ggsave(paste0(outdir,"/frateFit_EUE_NNE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/frateFit_EUE_NNE.png"),height=5,width=10)

wchngUE<- lm( wagechangeEUE_wave ~max.nempdur_wave + I(max.nempdur_wave^2)+ I(max.nempdur_wave^3), data=subset(DTall, UE_wave==T & midUE==F & seam==T) )

wchngUE_Fit<- array(NA,dim=c(24,2))
wchngUE_Dat<- array(NA,dim=c(24,2))
wchngUE_Dat[ , 1]<- seq_len(24)
wchngUE_Fit[ , 1]<- seq_len(24)
for(t in seq(1,24)){
	wchngUE_Dat[t,2] <- DTall[ max.nempdur_wave==t & UE_wave==T & midUE==F & seam==T, mean(wagechangeEUE_wave,na.rm = T)]
}
wchngUE_Fit[,2] = wchngUE_Fit[,1]*wchngUE$coefficients[2] +wchngUE_Fit[,1]^2*wchngUE$coefficients[3] + wchngUE_Fit[,1]^3*wchngUE$coefficients[4] + wchngUE$coefficients[1]


wchngNE<- lm( wagechangeENE_wave ~max.nempdur_wave + I(max.nempdur_wave^2)+ I(max.nempdur_wave^3), data=subset(DTall, NE_max==T & seam==T ) )

wchngNE_Fit<- array(NA,dim=c(24,2))
wchngNE_Dat<- array(NA,dim=c(24,2))
wchngNE_Dat[ , 1]<- seq_len(24)
wchngNE_Fit[ , 1]<- seq_len(24)
for(t in seq(1,24)){
	wchngNE_Dat[t,2] <- DTall[ max.nempdur_wave==t & seam==T & NE_max==T , mean(wagechangeENE_wave,na.rm = T)]
}
wchngNE_Fit[,2] = wchngNE_Fit[,1]*wchngNE$coefficients[2] +wchngNE_Fit[,1]^2*wchngNE$coefficients[3] + wchngNE_Fit[,1]^3*wchngNE$coefficients[4] + wchngNE$coefficients[1]

# plot these:
wchngFitPlot <- data.table(cbind(wchngUE_Fit,wchngNE_Fit[,2]))
names(wchngFitPlot) <- c("Duration","EUE","ENE")
wchngDatPlot <- data.table(cbind(wchngUE_Dat,wchngNE_Dat[,2]))
names(wchngDatPlot) <- c("Duration","EUE","ENE")
wchngFitMelt <- melt(wchngFitPlot, id.vars = "Duration")
wchngDatMelt <- melt(wchngDatPlot, id.vars = "Duration")

ggplot(data=wchngFitMelt, aes(x=Duration,y=value,color=variable)) + geom_line(size=1.5)+geom_hline(yintercept = 0.13,size=1.5)+
	theme_bw()+xlab("Completed Duration (mo)")+ylab("Earnings Change (pct)")+ylim(c(-.4,.15))+
	scale_color_manual( values = c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1,2,1,2)]) ,
						labels=c("Thru U","Thru N","",""))+
	theme(legend.title = element_blank(),
		  legend.position = c(0.2,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))+
	geom_point(data=wchngDatMelt, aes(x=Duration,y=value,color=variable))+geom_hline(yintercept = 0.13,size=1.5)

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

