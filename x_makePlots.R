# January 27, 2016
# Make occupation switching plots
# 1) wage change distribution for whole sample
# 2) EE wage change distribution
# 3) EUE wage change distribution
library(data.table)
library(ggplot2)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching/")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

wagechanges <- readRDS("./Data/balancedwagechanges.RData")

toKeep <- c("switchedOcc", 
	    "Young",
	    "HSCol",
		"recIndic",
	    "wagechange",
	    "wagechange_EUE", 
	    "wagechange_all", 
	    "balanceweight", 
	    "EE","EU","UE")

# select toKeep columns only
wagechanges <- wagechanges[, toKeep, with = FALSE]

# create group weight totals so each total within group will sum to 1
wagechanges <- wagechanges[!is.na(switchedOcc)]
wagechanges[, groupBalanceWeight_all := sum(balanceweight), by = switchedOcc]
wagechanges[!is.na(wagechange_EUE) & (EE | EU), groupBalanceWeight_EUE := sum(balanceweight), by = switchedOcc]
wagechanges[, groupBalanceWeight_rec := sum(balanceweight), by = "switchedOcc,recIndic"]
wagechanges[!is.na(wagechange_EUE) & (EE | EU), groupBalanceWeight_EUErec := sum(balanceweight), by = "switchedOcc,recIndic"]
wagechanges[EE == TRUE, groupBalanceWeight_EE := sum(balanceweight), by = switchedOcc]
wagechanges[EE == TRUE, groupBalanceWeight_EErec := sum(balanceweight), by = "switchedOcc,recIndic"]


# plot all earnings changes
ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_all, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,.8)) +
	xlim(c(-10,10)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_job.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_job.eps",height= 5,width=10)
#cdf
ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_all, linetype = switchedOcc)) +
	stat_ecdf()+
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	theme_bw() +
	theme(legend.position = c(.2,0.85)) +
	ylab("") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_job_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_job_cdf.eps",height= 5,width=10)	

ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_rec, color = recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,0.8)) +
	xlim(c(-10,10)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_joboccrec.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_joboccrec.eps",height= 5,width=10)

ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_rec, color = recIndic)) +
	geom_line(stat ="density",size = 0.5) +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,1.4)) +
	xlim(c(-10,10)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_jobrec.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_jobrec.eps",height= 5,width=10)

ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_rec, color = recIndic)) +
	stat_ecdf() +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	xlim(c(-10,10)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_jobrec_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_jobrec_cdf.eps",height= 5,width=10)


ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_rec, color = recIndic, linetype = switchedOcc)) +
	stat_ecdf() +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	xlab("Earnings change")+ ylab("")
ggsave(filename = "Figures/wagechanges_allrec_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_allrec_cdf.eps",height= 5,width=10)


# plot EUE wage changes
ggplot(subset(wagechanges, EE==T|EU==T), aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUE, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,1.5)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EUE.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EUE.eps",height= 5,width=10)

ggplot(subset(wagechanges, EE==T|EU==T), aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUE, linetype = switchedOcc)) +
	stat_ecdf() +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab(" ") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EUE_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EUE_cdf.eps",height= 5,width=10)


# plot EUE wage changes
ggplot(subset(wagechanges, EE==T|EU==T), aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUErec, color = recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,1.5)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EUErec.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EUErec.eps",height= 5,width=10)

ggplot(subset(wagechanges, EE==T|EU==T), aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUErec, color = recIndic, linetype = switchedOcc)) +
	stat_ecdf() +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EUErec_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EUErec_cdf.eps",height= 5,width=10)

# plot EE wage changes
ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EE, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,1.5)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EE.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EE.eps",height= 5,width=10)
# 
ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EE, linetype = switchedOcc)) +
	stat_ecdf() +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	theme_bw() +
	xlim(c(-4,4)) +
	theme(legend.position = c(0.2,0.85)) +
	ylab("") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EE_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EE_cdf.eps",height= 5,width=10)

ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EErec, color=recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,1.5)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EErec.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EErec.eps",height= 5,width=10)
ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EErec, color=recIndic, linetype = switchedOcc)) +
	stat_ecdf() +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("") +
	xlab("Earnings change")
ggsave(filename = "Figures/wagechanges_EErec_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_EErec_cdf.eps",height= 5,width=10)

## Now use full distribution ##---------------------------------
rm(wagechanges)

DTall <- readRDS("./Data/DTall_6.RData")

toKeep <- c(toKeep,"wpfinwgt")

# select toKeep columns only
DTall <- DTall[, toKeep, with = FALSE]
DTall <- subset(DTall, is.finite(wpfinwgt) & is.finite(wagechange_all))
# plot all wage changes
ggplot(DTall, aes(wagechange_all, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	ylim(c(0,.8)) +
	xlim(c(-10,10)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_unc.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_unc.eps",height= 5,width=10)
#cdf
ggplot(DTall, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_all, linetype = switchedOcc)) +
	stat_ecdf()+
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	theme_bw() +
	theme(legend.position = c(.2,0.85)) +
	ylab("") +
	xlab("Earnings change") 
ggsave(filename = "Figures/wagechanges_unc_cdf.png",height= 5,width=10)
ggsave(filename = "Figures/wagechanges_unc_cdf.eps",height= 5,width=10)	


rm(DTall)

# Now use wave-level data----------------------------------
DTseam <- readRDS("./Data/DTseam.RData")

ggplot(  subset(DTseam, (EE_wave|EU_wave|UE_wave)& !is.na(switchedOcc_wave)), aes(wagechange_wave, linetype = switchedOcc_wave)) +
	geom_density(stat ="density",size = 1.25) +
	ylim(c(0,1.2)) +
	xlim(c(-3,3)) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	theme_bw()+
	theme(legend.position = c(.2,0.85)) +
	xlab("Earnings growth")+ylab("")
ggsave(filename = "Figures/wave_distchng.png",height= 5,width=10)
ggsave(filename = "Figures/wave_distchng.eps",height= 5,width=10)	


ggplot(  subset(DTseam, (EE_wave|EU_wave|UE_wave) & switchedOcc_wave==T), aes(wagechange_wave, color = recIndic_wave)) +
	geom_density(stat ="density",size = 1.25,linetype="dashed") +
	ylim(c(0,1.2)) +
	xlim(c(-3,3)) +
	scale_color_manual(values = c(hcl(h=seq(30, 390, length=4), l=50, c=100)[1:2]) ,
						  labels = c("Expansion", "Recession"),
						  name = "") +
	theme_bw()+
	theme(legend.position = c(.2,0.85)) +
	xlab("Earnings growth")+ylab("")
ggsave(filename = "Figures/wave_distchng_swrec.png",height= 5,width=10)
ggsave(filename = "Figures/wave_distchng_swrec.eps",height= 5,width=10)	


ggplot(  subset(DTseam, (EE_wave|EU_wave|UE_wave) & switchedOcc_wave==F), aes(wagechange_wave, color = recIndic_wave)) +
	geom_density(stat ="density",size = 1.25) +
	ylim(c(0,1.2)) +
	xlim(c(-3,3)) +
	scale_color_manual(values = c(hcl(h=seq(30, 390, length=4), l=50, c=100)[1:2]) ,
					   labels = c("Expansion", "Recession"),
					   name = "") +
	theme_bw()+
	theme(legend.position = c(.2,0.85)) +
	xlab("Earnings growth")+ylab("")
ggsave(filename = "Figures/wave_distchng_noswrec.png",height= 5,width=10)
ggsave(filename = "Figures/wave_distchng_noswrec.eps",height= 5,width=10)	
