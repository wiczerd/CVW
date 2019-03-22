# x_mainFigures.R


library(data.table)
library(zoo)
library(stats)
library(Hmisc)
library(scales)
library(ggplot2)

library(xtable)

minEarn = 1040 
minLEarn = log(2*minEarn)

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
 

stackedDens <- function( DT,gg,wc,wt=NULL ){
	
	Ng <- DT[ is.finite(eval(as.name(gg)))==T , length(unique( eval(as.name(gg)) ))]
	#fraction in each group:
	#fracg <- DT[ is.finite(eval(as.name(gg)))==T, ftable( eval(as.name(gg)) )/sum(is.finite(eval(as.name(gg))) ) ]
	if( !is.null(wt) ){
		fracg <- DT[ is.finite(eval(as.name(gg)))==T, sum(eval(as.name(wt)) ), by = gg ]
		setkeyv(fracg, cols=c(gg) )
		fracg <- fracg$V1/sum(fracg$V1)
	}else{
		fracg <- DT[ is.finite(eval(as.name(gg)))==T, ftable( eval(as.name(gg)) )/sum(is.finite(eval(as.name(gg))) ) ]
	}
	
	wcgidensity <- array(0,dim=c(512,Ng*2))
	for( gi in seq(1,Ng)){
		if(!is.null(wt)){
			wtsum <- DT[ eval(as.name(gg)) == gi & is.finite(eval(as.name(wc))) , sum( eval(as.name(wt)))]
			tmp <- DT[ eval(as.name(gg)) == gi & is.finite(eval(as.name(wc))) , density( eval(as.name(wc)), weights = eval(as.name(wt))/wtsum ,n=512)]
		}else{
			tmp <- DT[ eval(as.name(gg)) == gi & is.finite(eval(as.name(wc))) , density( eval(as.name(wc)) ,n=512)]
		}
		wcgidensity[,(gi-1)*2+1] <- tmp$x
		wcgidensity[,(gi-1)*2+2] <- tmp$y*fracg[gi]
		
	}
	wc_condldensity <- array(0,dim=c(512*Ng,Ng+1))
	wc_condldensity[,1] <- rbind( wcgidensity[,seq(1,Ng*2-1,2)] )
	wc_condldensity[,1] <- sort(wc_condldensity[,1])
	for( gi in seq(1,Ng)){
		densinterp <- approx(wcgidensity[,((gi-1)*2+1)],wcgidensity[,((gi-1)*2+2)],wc_condldensity[,1],yleft=0,yright=0)
		wc_condldensity[,gi+1] <- densinterp$y
	}
	wc_totdensity <- wc_condldensity
	if(Ng>1){
		for( gi in seq(2,Ng)){
			wc_totdensity[,gi+1] <-wc_totdensity[,gi] +wc_condldensity[,gi+1]
		}
	}
	wc_dtable <- wc_totdensity#cbind(wc_totdensity,wc_condldensity[,2:(Ng+1)])
	wc_dtable <- data.table(wc_dtable)
	names(wc_dtable)[1] <-"WageChange"
	for(gi in seq(1,Ng)){
		names(wc_dtable)[gi+1] <- paste0("TotPct",gi)
	#	names(wc_dtable)[Ng + gi+1] <- paste0("CondlPct",gi)
	}
	return(wc_dtable)
}

DTseam <- readRDS(paste0(datadir,"/DTseam.RData"))


DTseam[ ,wtEUE:= truncweight]
DTseam[ UE_wave==T,wtEUE:= 0.]
DTseam[ EU_wave==T,wtEUE:= truncweight]

DTseam[ , last2.stable_emp:= shift(last.stable_emp), by=id]
DTseam[ , last3.stable_emp:= shift(last2.stable_emp), by=id]
DTseam[ wave-2 != shift(wave,2), last2.stable_emp := NA]
DTseam[ wave-3 != shift(wave,3), last3.stable_emp := NA]

DTseam[ EU_anan& changer_anan & is.na(matched_EUUE_anan),changer_anan:= NA]
DTseam[ UE_anan& changer_anan & is.na(matched_EUUE_anan),changer_anan:= NA]
DTseam[ EU_anan& changer_anan & is.na(matched_EUUE_anan),EU_anan:= NA]
DTseam[ UE_anan& changer_anan & is.na(matched_EUUE_anan),UE_anan:= NA]

DTseam[ , nextann.wavewage := levwage + shift(levwage,type="lead") + shift(levwage,2,type="lead")  , by=id]
DTseam[ wave+1!=shift(wave,type = "lead") | wave+2!=shift(wave,2,type = "lead") , nextann.wavewage := NA]
DTseam[ , nextann.wavewage := log(nextann.wavewage + (1+nextann.wavewage^2)^.5) ]


DTseam[  EE_anan==T & UE_anan==F & EU_anan ==F , g := 2]
DTseam[  EE_anan==F & UE_anan==T | EU_anan ==T , g := 1]
DTseam[!(EE_anan==T | UE_anan==T | EU_anan ==T), g := 3]
DTseam[ !(changer_anan==T|stayer_anan==T| (EU_wave==F & nextann.wavewage<= 0) ), g:=NA]
#DTseam[ g==3 & (!last.stable_emp==T | !last2.stable_emp | !last3.stable_emp), g:=NA]
DTseam[ lastann.wavewage<minLEarn , g:=NA]
DTseam[!is.na(DTseam$g), g1 := ifelse(g==1,1,0)]
DTseam[!is.na(DTseam$g), g2 := ifelse(g==2,1,0)]
DTseam[!is.na(DTseam$g), g3 := ifelse(g==3,1,0)]

wcananDens <- stackedDens(DTseam,"g","wagechange_anan", wt="truncweight")
wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[value>exp(-6) , logValue := log(value)]
wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=logValue,ymin=-6,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Annual-Annual Log Earnings Change (6-wave view)")+ylab("log density")+xlim(c(-3,3)) + 
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]), name="",label=c("stay","EE","EU,UE") ) +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/stacked_wagechange_anan.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange_anan.png"),height=5,width=10)

DTseam[ is.finite(EE_anan)& is.finite(EU_anan) & is.finite(UE_anan), K:= 1]
DTseam[ !(changer_anan==T|stayer_anan==T|nextann.wavewage<=0), K:=NA]
DTseam[ lastann.wavewage<minLEarn , K:=NA]
wcananDens <- stackedDens(DTseam,"K","wagechange_anan", wt="truncweight")
wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[value>exp(-6) , logValue := value]
#wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=logValue,ymin=0,x=WageChange))+geom_ribbon(color="darkred",fill="darkred")+
	theme_bw()+xlab("Annual-Annual Log Earnings Change")+ylab("density")+xlim(c(-3,3)) + 
	#scale_fill_manual(values= "darkred") +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/all_wagechange_anan.eps"),height=5,width=10)
ggsave(paste0(outdir,"/all_wagechange_anan.png"),height=5,width=10)


DTseam[ , g:=NULL]
DTseam[  EE_wave==T & UE_wave==F & EU_wave ==F , g := 2]
DTseam[  EE_wave==F & UE_wave==T | EU_wave ==T , g := 1]
DTseam[  EE_wave==T &(UE_wave==T | EU_wave ==T), g := id%%2 + 1]
DTseam[!(EE_anan==T | UE_anan==T | EU_anan ==T), g := 3]
DTseam[ !(changer_anan==T|stayer_anan==T|(EU_wave==F & nextann.wavewage<= 0)), g:=NA]
DTseam[ g==3 & (EU_anan | UE_anan | EE_anan), g:=NA]
DTseam[ lastann.wavewage<minLEarn , g:=NA]
DTseam[!is.na(DTseam$g), g1 := ifelse(g==1,1,0)]
DTseam[!is.na(DTseam$g), g2 := ifelse(g==2,1,0)]
DTseam[!is.na(DTseam$g), g3 := ifelse(g==3,1,0)]

wcananDens <- stackedDens(DTseam,"g","wagechange_anan", wt="truncweight")
wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[value>exp(-6) , logValue := log(value)]
wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=logValue,ymin=-6,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Annual-Annual Log Earnings Change ")+ylab("log density")+xlim(c(-3,3)) + 
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]), name="",label=c("stay","EE","EU,UE") ) +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/stacked_wagechange_anan_wv.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange_anan_wv.png"),height=5,width=10)

wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[value>exp(-6) , Value := value]
wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=Value,ymin=0,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Annual-Annual Log Earnings Change ")+ylab("Density")+xlim(c(-2.75,2.75)) + coord_cartesian(ylim=c(0,.3))+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:3)]), name="",label=c("stay","EE","EU,UE") ) +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/stacked_LevWagechange_anan_wv.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_LevWagechange_anan_wv.png"),height=5,width=10)



DTseam[ , g:=NULL]
DTseam[ (EE_wave==T & UE_wave==F & EU_wave ==F)&changer , g := 2]
DTseam[ (EE_wave==F & UE_wave==T | EU_wave ==T)&changer , g := 1]
DTseam[!(EU_wave==T | UE_wave==T | EE_wave ==T)&stayer  , g := 3]
DTseam[ !(changer==T|stayer==T|(EU_wave==F & nextann.wavewage<= 0)), g:=NA]
DTseam[ g==3 & (EU_anan | UE_anan | EE_anan), g:=NA]
DTseam[ lastann.wavewage<minLEarn , g:=NA]
DTseam[!is.na(DTseam$g), g1 := ifelse(g==1,1,0)]
DTseam[!is.na(DTseam$g), g2 := ifelse(g==2,1,0)]
DTseam[!is.na(DTseam$g), g3 := ifelse(g==3,1,0)]

wcananDens <- stackedDens(DTseam,"g","wagechangeEUE_wave", wt="wtEUE")
wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[value>exp(-6) , logValue := log(value)]
wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=logValue,ymin=-6,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wave-Wave Log Earnings Change")+ylab("Log Density")+xlim(c(-3,3)) + 
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]), name="",label=c("stay","EE","EU,UE") ) +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/stacked_wagechangeEUE_wave.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechangeEUE_wave.png"),height=5,width=10)

wcananMelt <- melt(wcananDens, id.vars = "WageChange")
wcananMelt[ , logValue := value]
wcananMelt[ , g:=4L-as.integer(variable)]
ggplot(subset(wcananMelt,is.finite(value)) ,aes(ymax=logValue,ymin=0,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wave-Wave Log Earnings Change")+ylab("Density")+xlim(c(-3,3)) + coord_cartesian(ylim=c(0,.3))+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]), name="",label=c("stay","EE","EU,UE") ) +
	geom_vline(xintercept = 0.) 
ggsave(paste0(outdir,"/stacked_LevWagechangeEUE_wave.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_LevWagechangeEUE_wave.png"),height=5,width=10)


#do it among stayers
DTseam[ , c("rank_w_tm1","pct_w_tm1") :=NULL]
DTseam[ truncweight>0 & (stayer_anan==T) & is.finite(lastann.wavewage)& nextann.wavewage>0, rank_w_tm1 := frank(lastann.wavewage) ]
DTseam[ truncweight>0 & (stayer_anan==T) & is.finite(lastann.wavewage)& nextann.wavewage>0, rank_w_tm1 := rank_w_tm1/max(rank_w_tm1,na.rm=T)]
DTseam[ truncweight>0 & is.finite(rank_w_tm1), pct_w_tm1 :=as.integer( round(100*rank_w_tm1))]
pct_w_tm1 <- data.table(DTseam[ , wtd.mean( wagechange_anan,weights=truncweight,na.rm=T), by=pct_w_tm1])
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.10,weights=truncweight,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.50,weights=truncweight,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.90,weights=truncweight,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")

names(pct_w_tm1) <- c("pct","mn_w_chng","P10_w_chng","P50_w_chng","P90_w_chng")
pct_w_tm1 <- melt(pct_w_tm1,id.vars="pct")
ggplot(pct_w_tm1, aes( x=pct,y=value, color=variable ))+theme_bw()+geom_smooth(span=.1,se=F)+
	scale_color_manual(labels=c("Mean Wage Change","10pct Wage Change","50pct Wage Change","90pct Wage Change"),
					   values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]))+ #ggtitle("Only 6-wave stable")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white"))+xlim(c(2,98))+ylim(c(-.7,2.1))
ggsave(paste0(outdir,"/pctwtm1_wc_stable.eps"),height=5,width=10)
ggsave(paste0(outdir,"/pctwtm1_wc_stable.png"),height=5,width=10)

#do it among anyone with high-enough earnings
DTseam[ , c("rank_w_tm1","pct_w_tm1") :=NULL]
DTseam[truncweight>0 & lastann.wavewage>minLEarn & is.finite(lastann.wavewage) &(changer_anan|stayer_anan) & (EU_wave==T|nextann.wavewage>0), rank_w_tm1 := frank(lastann.wavewage) ]
DTseam[truncweight>0 & lastann.wavewage>minLEarn & is.finite(lastann.wavewage) &(changer_anan|stayer_anan) & (EU_wave==T|nextann.wavewage>0), rank_w_tm1 := rank_w_tm1/max(rank_w_tm1,na.rm=T)]
DTseam[truncweight>0 & is.finite(rank_w_tm1), pct_w_tm1 :=as.integer( round(100*rank_w_tm1))]
pct_w_tm1 <- data.table(DTseam[ , wtd.mean( wagechange_anan,weights=perwt,na.rm=T), by=pct_w_tm1])
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.10,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.50,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.90,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")

names(pct_w_tm1) <- c("pct","mn_w_chng","P10_w_chng","P50_w_chng","P90_w_chng")
pct_w_tm1 <- melt(pct_w_tm1,id.vars="pct")
ggplot(pct_w_tm1, aes( x=pct,y=value, color=variable ))+theme_bw()+geom_smooth(span=.1,se=F)+ylab("Log Earnings Change") + xlab("Percentile")+
	scale_color_manual(labels=c("Mean Wage Change","10pct Wage Change","50pct Wage Change","90pct Wage Change"),
					   values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]))+#ggtitle("With earnings high enough")+
theme(legend.title = element_blank(),
	  legend.position = c(0.8,0.8),
	  legend.background = element_rect(linetype = "solid",color = "white")) + xlim(c(2,98))+ylim(c(-.7,2.1))
ggsave(paste0(outdir,"/pctwtm1_wc_earntr.eps"),height=5,width=10)
ggsave(paste0(outdir,"/pctwtm1_wc_earntr.png"),height=5,width=10)

#stayers on the overall percentiles 

DTseam[ , c("rank_w_tm1","pct_w_tm1") :=NULL]
DTseam[truncweight>0 & lastann.wavewage>minLEarn & is.finite(lastann.wavewage) &(changer_anan|(stayer_anan)  ) & nextann.wavewage>0, rank_w_tm1 := frank(lastann.wavewage) ]
DTseam[truncweight>0 & lastann.wavewage>minLEarn & is.finite(lastann.wavewage) &(changer_anan|(stayer_anan)  ) & nextann.wavewage>0, rank_w_tm1 := rank_w_tm1/max(rank_w_tm1,na.rm=T)]
DTseam[truncweight>0 & is.finite(rank_w_tm1), pct_w_tm1 :=as.integer( round(100*rank_w_tm1))]
pct_w_tm1 <- data.table(DTseam[ , wtd.mean( wagechange_anan,weights=perwt,na.rm=T), by=pct_w_tm1])
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.10,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.50,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")
pct_w_tm1 <- merge(pct_w_tm1,data.table(DTseam[ , wtd.quantile(wagechange_anan,prob=0.90,weights=perwt,na.rm=T), by=pct_w_tm1]), by = "pct_w_tm1")

names(pct_w_tm1) <- c("pct","mn_w_chng","P10_w_chng","P50_w_chng","P90_w_chng")
pct_w_tm1 <- melt(pct_w_tm1,id.vars="pct")
ggplot(pct_w_tm1, aes( x=pct,y=value, color=variable ))+theme_bw()+geom_smooth(span=.1,se=F)+
	scale_color_manual(labels=c("Mean Wage Change","10pct Wage Change","50pct Wage Change","90pct Wage Change"),
					   values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]))+#ggtitle("With earnings high enough")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.8),
		  legend.background = element_rect(linetype = "solid",color = "white")) + xlim(c(2,98))
ggsave(paste0(outdir,"/pctwtm1_wc_st_fullpct.eps"),height=5,width=10)
ggsave(paste0(outdir,"/pctwtm1_wc_st_fullpct.png"),height=5,width=10)


# Hamburger plots for occupational-wage ladder (lack-thereof)

# EE *********************
DTseam[changer ==T & is.finite(wagechangeEUE_wave) & (EE_wave) & switched_wave==T, rank_wcEUE_wave := frank(wagechangeEUE_wave)/.N ]
DTseam[, pct_wcEUE_wave :=as.integer( round(100*rank_wcEUE_wave),5)]

wc_occwc <- data.table(DTseam[ , wtd.mean(wagechangeEUE_wave,weights=truncweight,na.rm=T), by=pct_wcEUE_wave])
names(wc_occwc) <- c("pct_wcEUE_wave","wchng")
wc_occwc <- merge(wc_occwc,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.10,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwc) <- c("pct_wcEUE_wave","wchng","P10occhcng")
wc_occwc <- merge(wc_occwc,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.5,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwc) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng")
wc_occwc <- merge(wc_occwc,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.90,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwc) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng","P90occchng")
wc_occwc <- merge(wc_occwc,data.table(DTseam[ switched_wave==T, wtd.mean    (occwagechangeEUE_wave,          weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwc) <- c("pct","wchng","P10occhcng","P50occchng","P90occchng","Meanoccchng")

wc_occwc_meltEE <- melt(wc_occwc,id.vars="pct")
ggplot(wc_occwc_meltEE, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.15,se=F)+theme_bw()+
	scale_color_manual(labels=c("Earnings Change","10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c("grey",hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgEE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgEE.png"),height=5,width=10)
wc_occwc[ , wchng:=NULL]
wc_occwc_meltEE <- melt(wc_occwc,id.vars="pct")
ggplot(wc_occwc_meltEE, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.15,se=F)+theme_bw()+
	scale_color_manual(labels=c("10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgnoearnEE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgnoearnEE.png"),height=5,width=10)

DTseam[ , c("pct_wcEUE_wave","rank_wcEUE_wave"):=NULL]
rm(wc_occwc)
rm(wc_occwc_meltEE)
#EUE *********************

DTseam[is.finite(wagechangeEUE_wave) & EU_wave & changer & switched_wave==T , rank_wcEUE_wave := frank(wagechangeEUE_wave)/.N ]
DTseam[, pct_wcEUE_wave :=as.integer( round(100*rank_wcEUE_wave),5)]

wc_occwcEUE <- data.table(DTseam[EU_wave==T , wtd.mean(wagechangeEUE_wave,weights=truncweight,na.rm=T), by=pct_wcEUE_wave])
names(wc_occwcEUE) <- c("pct_wcEUE_wave","wchng")
wc_occwcEUE <- merge(wc_occwcEUE,data.table(DTseam[ EU_wave==T &switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.10,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcEUE) <- c("pct_wcEUE_wave","wchng","P10occhcng")
wc_occwcEUE <- merge(wc_occwcEUE,data.table(DTseam[ EU_wave==T &switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.5,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcEUE) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng")
wc_occwcEUE <- merge(wc_occwcEUE,data.table(DTseam[ EU_wave==T &switched_wave==T, wtd.quantile(occwagechangeEUE_wave,prob=0.90,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcEUE) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng","P90occchng")
wc_occwcEUE <- merge(wc_occwcEUE,data.table(DTseam[ EU_wave==T &switched_wave==T, wtd.mean    (occwagechangeEUE_wave,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcEUE) <- c("pct","wchng","P10occhcng","P50occchng","P90occchng","Meanoccchng")

wc_occwc_meltEUE <- melt(wc_occwcEUE,id.vars="pct")
ggplot(wc_occwc_meltEUE, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.15,se=F)+theme_bw()+
	scale_color_manual(labels=c("Earnings Change","10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c("grey",hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgEU.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgEU.png"),height=5,width=10)

wc_occwcEUE[ , wchng:=NULL]
wc_occwc_meltEUE <- melt(wc_occwcEUE,id.vars="pct")
ggplot(wc_occwc_meltEUE, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.15,se=F)+theme_bw()+
	scale_color_manual(labels=c("10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.2,0.85),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgnoearnEU.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgnoearnEU.png"),height=5,width=10)

DTseam[ , c("pct_wcEUE_wave","rank_wcEUE_wave"):=NULL]
rm(wc_occwcEUE)
rm(wc_occwc_meltEUE)

#stayers *********************
DTseam[stayer_anan ==T & stayer==T & switched_wave==T& is.finite(wagechangeEUE_wave) , rank_wcEUE_wave := frank(wagechangeEUE_wave)/.N ]
DTseam[, pct_wcEUE_wave :=as.integer( round(100*rank_wcEUE_wave),5)]

wc_occwcst <- data.table(DTseam[ , wtd.mean(wagechangeEUE_wave,weights=truncweight,na.rm=T), by=pct_wcEUE_wave])
names(wc_occwcst) <- c("pct_wcEUE_wave","wchng")
wc_occwcst <- merge(wc_occwcst,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechange_wave,prob=0.10,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcst) <- c("pct_wcEUE_wave","wchng","P10occhcng")
wc_occwcst <- merge(wc_occwcst,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechange_wave,prob=0.5,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcst) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng")
wc_occwcst <- merge(wc_occwcst,data.table(DTseam[ switched_wave==T, wtd.quantile(occwagechange_wave,prob=0.90,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcst) <- c("pct_wcEUE_wave","wchng","P10occhcng","P50occchng","P90occchng")
wc_occwcst <- merge(wc_occwcst,data.table(DTseam[ switched_wave==T, wtd.mean(occwagechange_wave,weights=truncweight,na.rm=T), by=pct_wcEUE_wave]), by = "pct_wcEUE_wave")
names(wc_occwcst) <- c("pct","wchng","P10occhcng","P50occchng","P90occchng","mean")

wc_occwc_meltst <- melt(wc_occwcst,id.vars="pct")
ggplot(wc_occwc_meltst, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.1,se=F)+theme_bw()+
	scale_color_manual(labels=c("Earnings Change","10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c("grey",hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.2),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgst.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgst.png"),height=5,width=10)


wc_occwcst[ , wchng:=NULL]
wc_occwc_meltst <- melt(wc_occwcst,id.vars="pct")
ggplot(wc_occwc_meltst, aes( x=pct,y=value, color=variable ))+geom_smooth(span=.15,se=F)+theme_bw()+
	scale_color_manual(labels=c("10pct Occ Earnings Change","50pct Occ Earnings Change","90pct Occ Earnings Change","Mean Occ Earnings Change","Mean Occ Earnings Change"),
					   values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)],"black"))+ylab("")+xlab("Earnings Growth Percentile")+
	theme(legend.title = element_blank(),
		  legend.position = c(0.8,0.13),
		  legend.background = element_rect(linetype = "solid",color = "white"))
ggsave(paste0(outdir,"/occwgnoearnst.eps"),height=5,width=10)
ggsave(paste0(outdir,"/occwgnoearnst.png"),height=5,width=10)

DTseam[ , c("pct_wcEUE_wave","rank_wcEUE_wave"):=NULL]
rm(wc_occwcst)
rm(wc_occwc_meltst)



toKeep <- c("truncweight","cycweight","wpfinwgt","EU_wave","UE_wave","EE_wave","switchedOcc_wave","wagechange_wave","wagechangeEUE_wave", "wagechange_wvan","recIndic2_wave","date")
DTseam <- subset(DTseam, is.finite(EU_wave) & is.finite(UE_wave)& is.finite(EE_wave))

# Occupation switchers and not--------------

DTseam[  DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F , g := 1]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==T & DTseam$EU_wave ==F , g := 2]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T , g := 3]
DTseam[!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T), g := 4]
DTseam[!is.na(DTseam$g), g1 := ifelse(g==1,1,0)]
DTseam[!is.na(DTseam$g), g2 := ifelse(g==2,1,0)]
DTseam[!is.na(DTseam$g), g3 := ifelse(g==3,1,0)]
DTseam[!is.na(DTseam$g), g4 := ifelse(g==4,1,0)]

#stackedCdist <- stackedDist(DTseam,"g","wagechange_wave","pct1000_wagechange_wave")
#DTseam[ g==2 & wagechange_wave< -5, g:=NA]
stayerSwitcherDens <- stackedDens(DTseam,"g","wagechange_wvan")
stayerSwitcherMelt <- melt(stayerSwitcherDens, id.vars = "WageChange")
stayerSwitcherMelt[value>exp(-15) , logValue := log(value)]
stayerSwitcherMelt[ , g:=5L-as.integer(variable)]
ggplot(subset(stayerSwitcherMelt,is.finite(logValue)) ,aes(ymax=logValue,ymin=-15,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wage Change")+ylab("Stacked log density")+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=5), l=50, c=100)[c(1:4)]), name="")
ggsave(paste0(outdir,"/stacked_wagechange.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechange.png"),height=5,width=10)


#EE and stayer:
DTseam[ g==1 , gEE := 1]
DTseam[ g==4 , gEE := 2]
stayerEEDens <- stackedDens(subset(DTseam, !(EU_wave|UE_wave)),"gEE","wagechange_wave")
stayerEEmelt <- melt(stayerEEDens,id.vars = "WageChange")
stayerEEmelt[ WageChange>-2 & WageChange<2, logValue:=log(value)]
stayerEEmelt[ , g:=3L-as.integer(variable)]
ggplot(subset(stayerEEmelt,is.finite(logValue)) ,aes(ymax=logValue,ymin=-4,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wage Change")+ylab("Stacked log density")+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:2)]), name="", label=c("Stay","EE"))
ggsave(paste0(outdir,"/stacked_wagechangeEE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechangeEE.png"),height=5,width=10)

#EE, EUE & stayer
DTseam[ , g:=NULL]
DTseam[  DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F , g := 1]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==T & DTseam$EU_wave ==T , g := 2]
DTseam[!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T), g := 3]
stayerEUEDens <- stackedDens(subset(DTseam, !UE_wave==T & is.finite(wagechangeEUE_wave)),"g","wagechangeEUE_wave")
stayerEUEmelt <- melt(stayerEUEDens,id.vars = "WageChange")
stayerEUEmelt[ WageChange>-3 & WageChange<3, logValue:=log(value)]
stayerEUEmelt[ , g:=4L-as.integer(variable)]
ggplot(subset(stayerEUEmelt,is.finite(logValue)) ,aes(ymax=logValue,ymin=-8,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wage Change")+ylab("Stacked log density")+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)]), name="", label=c("Stay","EE","EUE"))
ggsave(paste0(outdir,"/stacked_wagechangeEUE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stacked_wagechangeEUE.png"),height=5,width=10)

ggplot(subset(stayerEUEmelt,is.finite(logValue)) ,aes(ymax=value,ymin=0,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wage Change")+ylab("Stacked density")+
	scale_fill_manual(values=c(hcl(h=seq(15, 375, length=4), l=50, c=100)[c(1:3)]), name="", label=c("Stay","EE","EUE"))
ggsave(paste0(outdir,"/stackedlev_wagechangeEUE.eps"),height=5,width=10)
ggsave(paste0(outdir,"/stackedlev_wagechangeEUE.png"),height=5,width=10)


# Recession and expansion:
DTseam[  DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F & recIndic2_wave==F, g := 1]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T & recIndic2_wave==F, g := 2]
DTseam[!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T)& recIndic2_wave==F, g := 3]
DTseam[  DTseam$EE_wave==T & DTseam$UE_wave==F & DTseam$EU_wave ==F & recIndic2_wave==T, g := 4]
DTseam[  DTseam$EE_wave==F & DTseam$UE_wave==F & DTseam$EU_wave ==T & recIndic2_wave==T, g := 5]
DTseam[!(DTseam$EE_wave==T | DTseam$UE_wave==T | DTseam$EU_wave ==T)& recIndic2_wave==T, g := 6]
stayerEUEDens <- stackedDens(subset(DTseam, !UE_wave==T & is.finite(wagechangeEUE_wave)),"g","wagechangeEUE_wave")
stayerEUEmelt <- melt(stayerEUEDens,id.vars = "WageChange")
stayerEUEmelt[ WageChange>-4 & WageChange<4, logValue:=log(value)]
stayerEUEmelt[ , g:=7L-as.integer(variable)]
ggplot(subset(stayerEUEmelt,is.finite(logValue)) ,aes(ymax=logValue,ymin=-8,x=WageChange,fill=as.factor(g)))+geom_ribbon()+
	theme_bw()+xlab("Wage Change")+ylab("Stacked log density")+
	scale_fill_manual(values=hcl(h=seq(15, 3/4*(375-15)-15, length=3), l=c(50,50,50,80,80,80), c=100), name="", label=c("Stay, expansion","EE, expansion","EUE, expansion",
																									 "Stay, recession","EE, recession","EUE, recession"))


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

