
# May 26, 2015
# Calculate wage change statistics:
# 1) mean, standard deviation, median wage change for non-occupation-switch job changes
# 2) mean, standard deviation, median wage change for EE, UE, pooled occupation switches
# 3) quantile regressions
# 4) fraction of workers with positive and negative wage changes
# 5) correlation between unemployment rate and fraction of positive changes
library(Hmisc)
library(dplyr)
library(ggplot2)
library(xlsx)
library(stargazer)
library(quantreg)
library(McSpatial)
library(oaxaca)
library(reshape2)

setwd("~/workspace/CVW/R")

# Use 1 digit occupations from CEPR? (soc2d)
useSoc2d <- T
useRegResid <- T

# Read unemployment data
haver <- read.xlsx("./Data/unrate.xlsx", sheetName = "data", 
				   startRow = 2, colIndex = 2:4)
# Change date to first of the month for merging
haver <- haver %>%
	mutate(month = format(date, "%m"),
		   year = format(date, "%Y"),
		   date = as.Date(paste(month, "/1/", year, sep=""), "%m/%d/%Y")) %>%
	select(-year, -month)

toKeep <- c("wpfinwgt", "switchedOcc", "EE", "UE",  
			"wageChange", "wageChange_wU", "lfStat", "date",
			"wageChange_q", "wageChange_q_wU")

detach("package:xlsx")
detach("package:xlsxjars")
detach("package:rJava")


# Wage change distributions ---------------------------------------------

if(useSoc2d & useRegResid) {
	setwd("./Data/soc2d/RegResid")
} else if(useSoc2d & !useRegResid){
	setwd("./Data/soc2d/Raw")
} else if(!useSoc2d & useRegResid){
	setwd("./Data/occ/RegResid")
} else if(!useSoc2d & !useRegResid){
	setwd("./Data/occ/Raw")
}
wageChanges<-readRDS("wageChanges.RData")   


# exclude outliers
resChangeOutlier_rec <- quantile(wageChanges$wageChange_wU[wageChanges$recIndic == 1],probs=c(0.025,0.975),na.rm=T)
resChangeOutlier_exp <- quantile(wageChanges$wageChange_wU[wageChanges$recIndic == 0],probs=c(0.025,0.975),na.rm=T)
#wageChanges$Out <- (wageChanges$wageChange_wU<resChangeOutlier[1] |
#						wageChanges$wageChange_wU>resChangeOutlier[2] |
#						is.na(wageChanges$wageChange_wU) )
wageChanges$Out <- (  (wageChanges$wageChange_wU<resChangeOutlier_rec[1] |
					   	wageChanges$wageChange_wU>resChangeOutlier_rec[2]) 
					  & wageChanges$Rec == 1 ) |
	(  (wageChanges$wageChange_wU<resChangeOutlier_exp[1] |
			wageChanges$wageChange_wU>resChangeOutlier_exp[2]) 
			& wageChanges$Rec == 0 )
resChangeOutlier_q_rec <- quantile(wageChanges$wageChange_q_wU[wageChanges$Rec == 1],probs=c(0.025,0.975),na.rm=T)
resChangeOutlier_q_exp <- quantile(wageChanges$wageChange_q_wU[wageChanges$Rec == 0],probs=c(0.025,0.975),na.rm=T)
wageChanges$Out_q <- (  (wageChanges$wageChange_q_wU<resChangeOutlier_q_rec[1] |
						 	wageChanges$wageChange_q_wU>resChangeOutlier_q_rec[2]) 
						& wageChanges$Rec == 1 ) |
	(  (wageChanges$wageChange_q_wU<resChangeOutlier_q_exp[1] |
			wageChanges$wageChange_q_wU>resChangeOutlier_q_exp[2]) 
			& wageChanges$Rec == 0 )
wageChanges_In <- subset(wageChanges,!wageChanges$Out | is.na(!wageChanges$Out))


wageChangesLong <- melt(wageChanges_In, id.vars = c("id", "date", "Rec"),
						measure.vars = c("wageChange_wU", "wageChange_q_wU"))
wageChangesLong <- subset(wageChangesLong, !is.na(Rec))

# format variables for better plotting
wageChangesLong$Rec <- as.factor(wageChangesLong$Rec)
levels(wageChangesLong$Rec) <- c("Expansion", "Recession")
levels(wageChangesLong$variable) <- c("Monthly", "Quarterly")

# exp/rec panes
#postscript(file = "./Figures/wageChangeDensityRecExp.eps", width = 782, height = 569)
#ggplot(wageChangesLong, aes(value, fill = variable)) +
#        geom_density(alpha = 0.5) +
#        facet_grid(. ~ Rec) +
#        xlim(c(-3.75, 3.75)) +
#        ggtitle("Wage change distribution in expansion & recession") +
#        labs(fill = "Aggregation")
#dev.off()

# monthly/quarterly panes
png(file="./Figures/wageChangeDensityMthQtr.png", width = 782, height = 569)
ggplot(wageChangesLong, aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	facet_grid(. ~ variable) +
	xlim(c(-3.75, 3.75)) +
	ggtitle("Wage change distribution by aggregation") +
	labs(fill = "Business cycle")
dev.off()

# quarterly pane
png(file="./Figures/wageChangeDensityQtr.png",onefile=FALSE, width = 782, height = 569)
ggplot(subset(wageChangesLong,wageChangesLong$variable=="Quarterly"), aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	xlim(c(-3.75, 3.75)) +
	ggtitle("Wage change distribution, Quarterly") +
	labs(fill = "Business cycle")
dev.off()

# monthly pane
png(file="./Figures/wageChangeDensityMth.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Monthly"), aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	xlim(c(-2.1, 2.1)) +
	ggtitle("Wage change distribution, Monthly") +
	labs(fill = "Business cycle")
dev.off()

#box plot of the same thing
png(file="./Figures/wageChangeBoxQtr.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Quarterly"), aes(y = value, x = Rec)) +
	geom_boxplot() +
	guides(fill=F) +
	stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
	ggtitle("Wage change distribution, Quarterly") +
	labs(fill = "Business cycle")
dev.off()

png(file="./Figures/wageChangeBoxMth.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Monthly"), aes(y = value, x = Rec)) +
	geom_boxplot() +
	guides(fill=F) +
	stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
	ggtitle("Wage change distribution, Monthly") +
	labs(fill = "Business cycle")
dev.off()



wageChangesLong <- melt(wageChanges_In, id.vars = c("id", "date", "Rec"),
						measure.vars = c("wageChange", "wageChange_q"))
wageChangesLong <- subset(wageChangesLong, !is.na(Rec))

# format variables for better plotting
wageChangesLong$Rec <- as.factor(wageChangesLong$Rec)
levels(wageChangesLong$Rec) <- c("Expansion", "Recession")
levels(wageChangesLong$variable) <- c("Monthly", "Quarterly")


# monthly pane
png(file="./Figures/wageChangeDensityMth_woU.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Monthly"), aes(value, fill = Rec)) +
	geom_density(alpha = 0.5) +
	xlim(c(-2.1, 2.1)) +
	ggtitle("Wage change distribution  excluding EU, Monthly") +
	labs(fill = "Business cycle")
dev.off()


png(file="./Figures/wageChangeBoxMth_woU.png", width = 782, height = 569)
ggplot( subset(wageChangesLong,wageChangesLong$variable=="Monthly" & wageChangesLong$value >-1 & wageChangesLong$value<1), aes(y = value, x = Rec)) +
	geom_boxplot() +
	guides(fill=F) +
	stat_summary(fun.y=mean, geom="point", shape=5, size=4) +
	ggtitle("Wage change distribution excluding EU, Monthly") +
	labs(fill = "Business cycle")
dev.off()

## Decompositions:
# # subset for recessions and expansions
# wageChangesRec <- subset(wageChanges,wageChanges$Rec)
# wageChangesExp <- subset(wageChanges,!wageChanges$Rec)
# 
# wageChangesRecLong <- melt(wageChangesRec, id.vars = c("id", "date"),
#                          measure.vars = c("wageChange_wU", "wageChange_q_wU"))
# wageChangesExpLong <- melt(wageChangesExp, id.vars = c("id", "date"),
#                            measure.vars = c("wageChange_wU", "wageChange_q_wU"))

resChangeOutlier <- quantile(wageChanges$wageChange_wU,probs=c(0.025,0.975),na.rm=T)
wageChanges$Out <- (wageChanges$wageChange_wU<resChangeOutlier[1] |
						wageChanges$wageChange_wU>resChangeOutlier[2] |
						is.na(wageChanges$wageChange_wU) )
wageChangesIn <- subset(wageChanges,!Out)
wageChangesIn.kde <- density(wageChangesIn$wageChange_wU,na.rm=T)
plot(wageChangesIn.kde)

wageChangesIn.rec <- qplot(wageChange_wU,na.rm=T,data=wageChangesIn, colour = Rec, geom="density")
plot(wageChangesIn.rec)
wageChangesIn.kde <- density(wageChangesIn$wageChange_wU,na.rm=T)
plot(wageChangesIn.kde)

bp_wU<-boxplot(wageChange_wU~Rec,data=wageChanges,names=c("Expansion","Recession"))
title("Wage change distribution when changing jobs or into unemployment")

bp<-boxplot(wageChange~Rec,data=wageChanges,names=c("Expansion","Recession"))
title("Wage change distribution, excluding unemployment")

# Machado - Mata Decomposition ----------------------------------------

# do the regressions I'll use
#wageChanges <- subset(wageChanges, !wageChanges$Out)
wageChangesRec <- subset(wageChanges,recIndic & (UE|EE) & is.finite(wageChange) & is.finite(switchedOcc) &  is.finite(UE))
wageChangesExp <- subset(wageChanges,!recIndic & (UE|EE) & is.finite(wageChange) & is.finite(switchedOcc) & is.finite(UE))

qtl_delw <- c(0.1,0.25,0.5,0.75,0.9)
Nqtl = length(qtl_delw)
qr.exp <- list()
for(ti in seq(1:Nqtl)){
	qr.exp[[ti]] <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw[ti], weights= wpfinwgt,  data=wageChangesExp)
}
qr.rec <- list()
for(ti in seq(1:Nqtl)){
	qr.rec[[ti]] <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw[ti], weights= wpfinwgt,  data=wageChangesRec )
}

#output these regressions:
tab<-stargazer(qr.exp[1:Nqtl],  out="qr_exp.tex",
			   dep.var.labels.include = FALSE, no.space=T,  omit.stat=c("rsq","f","adj.rsq"),
			   title            = "Earnings change quantile regression over the cycle",
			   digits=2, column.sep.width= "1pt",df=F, flip=T,
			   covariate.labels = c("Const", "Switched","Unep Indic"))

tab<-stargazer(qr.rec[1:Nqtl],  out="qr_rec.tex",
			   dep.var.labels.include = FALSE, no.space=T,  omit.stat=c("rsq","f","adj.rsq"),
			   title            = "Earnings change quantile regression over the cycle",
			   digits=2, column.sep.width= "1pt",df=F, flip=T,
			   covariate.labels = c("Const", "Switched","Unep Indic"))

qr.exp.tog <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesExp)
qr.rec.tog <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesRec)
qr_rec_exp<-cbind(t(qr.exp.tog$coefficients),t(qr.rec.tog$coefficients))

rownames(qr_rec_exp)<- c("0.10","0.25","0.50","0.75","0.9")
colnames(qr_rec_exp)<- c("Const","Occ Sw", "Unemp Indic","Const","Occ Sw", "Unemp Indic")


qr_rec_exp.xt<-xtable(qr_rec_exp,label="tab:qr_rec_exp",digits=3,caption="Earnings change quantile regression coefficients")
print(qr_rec_exp.xt,file="qr_rec_exp.tex",hline.after=c(-1,-1,0,nrow(qr_rec_exp.xt)) )


#test joint equality of coefficients
wCE_0 <- with(wageChangesExp, data.frame(wageChange, switchedOcc, UE, wpfinwgt ))
wCR_0 <- with(wageChangesRec, data.frame(wageChange, switchedOcc, UE, wpfinwgt ))
nobsE <- nrow(wCE_0)
nobsR <- nrow(wCR_0)
Nsim = 500
coSwR_ni <- matrix(0.,nrow=Nsim,ncol=length(qtl_delw) )
coSwE_ni <- matrix(0.,nrow=Nsim,ncol=length(qtl_delw) )
# initialize the count of test
N_EleR_all = 0 # contrapositive of \exists \tau s.t. R<E
N_EleR_any = 0 # contrapositive of \forall \tau R< E
N_EleR_pt  = matrix(0.,nrow=1,ncol=length(qtl_delw))
for(ni in 1:Nsim){
	# sample from the two distributions (with replacement) run the quantile regression, then run Qu's test
	wCE_ni <- wCE_0[sample(nobsE, nobsE, replace=T ), ]
	wCR_ni <- wCR_0[sample(nobsR, nobsR, replace=T ), ]
	wReg.e_ni <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wCE_ni, ci=F)
	wReg.r_ni <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wCR_ni, ci=F)
	coSwR_ni[ni,] <- wReg.r_ni$coefficients[2,]
	coSwE_ni[ni,] <- wReg.e_ni$coefficients[2,]
	
	# test of equality
	test_ni = as.integer(coSwE_ni[ni,]-coSwR_ni[ni,] < 0)
	if( sum(test_ni)==length(qtl_delw) ){
		N_EleR_all = N_EleR_all+1	
	}
	test_ni = as.integer( coSwE_ni[ni,]-coSwR_ni[ni,] < 0 )
	if( sum(test_ni)>0 ){
		N_EleR_any = N_EleR_any+1
	}
	N_EleR_pt = N_EleR_pt + test_ni
	
	# Machado Mata decomposition
	fit.e_ni <- wReg.e_ni$fitted.values
	fit.r_ni <- wReg.r_ni$fitted.values
	
}
save.image()

# these are the quantile regressions we will be decomposing:
qtl_delw <- seq(0.1, .9, .1)
mm_rq.expansion <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesExp)
mm_rq.recession <- rq(wageChange ~ switchedOcc +UE,tau = qtl_delw, weights= wpfinwgt,  data=wageChangesRec )
mm_lm.expansion <- lm(wageChange ~ switchedOcc +UE,weights= wpfinwgt,  data=wageChangesExp)
mm_lm.recession <- lm(wageChange ~ switchedOcc +UE,weights= wpfinwgt,  data=wageChangesRec)


# setup regressions, but use weighted regression:
mm_bmat.exp <-rq(wageChange ~ switchedOcc + UE , taumat=qtl_delw, data=wageChangesExp, graphb = F)
mm_bmat.rec <-rq(wageChange ~ switchedOcc + UE , taumat=qtl_delw, data=wageChangesRec, graphb = F)
#(mm_bmat.exp - t(mm_rq.expansion$coefficients))/t(mm_rq.expansion$coefficients)
#(mm_bmat.rec - t(mm_rq.recession$coefficients))/t(mm_rq.recession$coefficients)
for(ki in 1:ncol(mm_bmat.exp)){
	for(ti in 1:nrow(mm_bmat.exp)){
		mm_bmat.exp[ti,ki] = mm_rq.expansion$coefficients[ki,ti]
		mm_bmat.rec[ti,ki] = mm_rq.recession$coefficients[ki,ti]
	}
}
qr_coefs <-data.frame(qtl_delw,mm_bmat.exp,mm_bmat.rec)
ggA <- ggplot( qr_coefs, aes(y = switchedOccTRUE, x = qtl_delw)) +
	geom_line(size = 2) + geom_point() + geom_hline(aes(yintercept=mm_lm.expansion$coefficients[2]), size = 2)
ggA <- ggA + geom_line( aes(y = switchedOccTRUE.1, x = qtl_delw), colour = "red", size = 2) +
	geom_hline(aes(yintercept=mm_lm.recession$coefficients[2]), colour="red", size = 2) +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of Occupation Switching"))
ggsave("./qtl_swocc.eps",width = 6, height = 4)
ggA
ggB <- ggplot( qr_coefs, aes(y = UETRUE, x = qtl_delw)) +
	geom_line( size = 2) + geom_point() + geom_hline(aes(yintercept=mm_lm.expansion$coefficients[3]), size = 2)
ggB <- ggB + geom_line( aes(y = UETRUE.1, x = qtl_delw), colour = "red", size = 2) +
	geom_hline(aes(yintercept=mm_lm.recession$coefficients[3]), colour="red", size = 2) +
	labs(list(x="Quantile", y="Wage Change Effect", title="Effect of EUE"))
ggsave("qtl_eue.eps",width = 6, height = 4)
ggB
set.seed(941987)
mm_wageChanges <- qregsim2(wageChange ~ switchedOcc + UE , ~ switchedOcc + UE , wageChangesExp, wageChangesRec,
						   mm_bmat.exp, mm_bmat.rec, timenames=c("Expansion","Recession"),nsim=14700) # sample size is ~ same as data size 
dat_mm <- data.frame(mm_wageChanges)
dat_mm <- dat_mm %>%
	mutate(pct1211 = d1211/d2211 ) %>%
	mutate(pct2212 = d2212/d2211 ) %>%
	mutate(pct1211 = as.numeric(ifelse(abs(pct1211) > 2,NA,pct1211)) ) %>%
	mutate(pct2212 = as.numeric(ifelse(abs(pct2212) > 2,NA,pct2212)) )
med.pct1211 <- quantile(dat_mm$pct2212,na.rm=T, probs=.5)
ggC <- ggplot(dat_mm, aes(x=ytarget,y=pct2212) )	 + geom_smooth(size=2) +
	labs(list(x="Wage Change",y="Flows' contribution",title="The change in distributions due to changed worker flows"))
#		geom_hline(aes(yintercept = med.pct1211), size=2)
ggC
ggsave("./Figures/mm_pctflows.png",width = 5, height = 5)
ggsave("./Figures/mm_pctflows.eps",width = 5, height = 5)
for(ki in 1:7){
	mm_wageChanges_swonly <- qregsim2(wageChange ~ switchedOcc + UE , ~ switchedOcc , wageChangesRec, wageChangesExp,
								  mm_bmat.rec, mm_bmat.exp, timenames=c("Recession","Expansion"),nsim=14700) 
}
dat_mm_swonly <- data.frame(mm_wageChanges_swonly)
dat_mm_swonly <- dat_mm_swonly %>%
	mutate(pct1211 = d1211/d2211 ) %>%
	mutate(pct2212 = d2212/d2211 ) %>%
	mutate(pct1211 = as.numeric(ifelse(abs(pct1211) > 2,NA,pct1211)) ) %>%
	mutate(pct2212 = as.numeric(ifelse(abs(pct2212) > 2,NA,pct2212)) )
med.pct1211 <- quantile(dat_mm_swonly$pct2212,na.rm=T, probs=.5)
ggD <- ggplot(dat_mm_swonly, aes(x=ytarget,y=pct2212) )	 + geom_smooth(size=2) +
	labs(list(x="Wage Change",y="Switch flows' contribution",title="Change in distributions due to changed switching probabilities"))
#		geom_hline(aes(yintercept = med.pct1211), size=2)
ggD
ggsave("./Figures/mm_pctswflows1.png",width = 5, height = 5)
ggsave("./Figures/mm_pctswflows1.eps",width = 5, height = 5)

#oaxaca decomp
wageChanges_oa <-subset(wageChanges, (UE | EE) ) 
#wageChanges_oa$LRec <- as.logical(ifelse(wageChanges_oa$LRec != F | wageChanges_oa$LRec != F, NA, wageChanges_oa$LRec))
wageChanges_oa <-subset(wageChanges_oa,  is.finite(Rec) & (Rec==1 | Rec==0) & is.finite(wageChange) & is.finite(switchedOcc) & is.finite(UE))
oa_wageChanges <- oaxaca(wageChange ~ switchedOcc + 
						 	UE | Rec,  data=wageChanges_oa, R=30)

oaxaca.results.1 <- oaxaca(ln.real.wage ~ age + female + LTHS + some.college +
						   	college + advanced.degree | foreign.born,
						   data = chicago, R = 30)