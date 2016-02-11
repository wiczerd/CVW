# January 27, 2016
# Make occupation switching plots
# 1) wage change distribution for whole sample
# 2) EE wage change distribution
# 3) EUE wage change distribution
library(data.table)
library(ggplot2)

setwd("G:/Research_Analyst/Eubanks/Occupation Switching/")

wagechanges <- readRDS("./Data/balancedwagechanges.RData")

toKeep <- c("switchedOcc", 
	    "Young",
	    "HSCol",
		"recIndic",
	    "wagechange",
	    "wagechange_EUE", 
	    "wagechange_all", 
	    "balanceweight", 
	    "EE")

# select toKeep columns only
wagechanges <- wagechanges[, toKeep, with = FALSE]

# create group weight totals so each total within group will sum to 1
wagechanges <- wagechanges[!is.na(switchedOcc)]
wagechanges[, groupBalanceWeight_all := sum(balanceweight), by = switchedOcc]
wagechanges[!is.na(wagechange_EUE), groupBalanceWeight_EUE := sum(balanceweight), by = switchedOcc]
wagechanges[, groupBalanceWeight_rec := sum(balanceweight), by = "switchedOcc,recIndic"]
wagechanges[!is.na(wagechange_EUE), groupBalanceWeight_EUErec := sum(balanceweight), by = "switchedOcc,recIndic"]
wagechanges[EE == TRUE, groupBalanceWeight_EE := sum(balanceweight), by = switchedOcc]
wagechanges[EE == TRUE, groupBalanceWeight_EErec := sum(balanceweight), by = "switchedOcc,recIndic"]


# plot all wage changes
ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_all, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,0.7)) +
	xlim(c(-12,12)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, all") 
ggsave(filename = "Figures/wagechanges_all.png")
ggsave(filename = "Figures/wagechanges_all.eps")

ggplot(wagechanges, aes(wagechange_all, weight = balanceweight/groupBalanceWeight_rec, color = recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,0.7)) +
	xlim(c(-12,12)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, all") 
ggsave(filename = "Figures/wagechanges_allrec.png")
ggsave(filename = "Figures/wagechanges_allrec.eps")

# plot EUE wage changes
ggplot(wagechanges, aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUE, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,1.27)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, EUE")
ggsave(filename = "Figures/wagechanges_EUE.png")
ggsave(filename = "Figures/wagechanges_EUE.eps")

# plot EUE wage changes
ggplot(wagechanges, aes(wagechange_EUE, weight = balanceweight/groupBalanceWeight_EUErec, color = recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,1.27)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, EUE")
ggsave(filename = "Figures/wagechanges_EUErec.png")
ggsave(filename = "Figures/wagechanges_EUErec.eps")

# plot EE wage changes
ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EE, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
			      labels = c("Non-switchers", "Switchers"),
			      name = "") +
	ylim(c(0,1.27)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, EE")
ggsave(filename = "Figures/wagechanges_EE.png")
ggsave(filename = "Figures/wagechanges_EE.eps")

ggplot(subset(wagechanges,EE), aes(wagechange, weight = balanceweight/groupBalanceWeight_EErec, color=recIndic, linetype = switchedOcc)) +
	geom_line(stat ="density",size = 0.5) +
	scale_linetype_manual(values = c("solid", "dashed"),
						  labels = c("Non-switchers", "Switchers"),
						  name = "") +
	scale_color_manual(values = c("black", "red"),
					   labels = c("Expansion", "Recession"),
					   name = "") +
	ylim(c(0,1.27)) +
	xlim(c(-4,4)) +
	theme_bw() +
	theme(legend.position = c(0.2,0.85)) +
	ylab("Density") +
	xlab("Wage change, EE")
ggsave(filename = "Figures/wagechanges_EErec.png")
ggsave(filename = "Figures/wagechanges_EErec.eps")