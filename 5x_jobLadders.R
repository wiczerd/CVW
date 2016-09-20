library(data.table)
library(zoo)
library(stats)
library(ggplot2)

wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
rootdir <- wd0
figuredir <- paste0(rootdir, "/Figures")
setwd(rootdir)

DTall <- readRDS(paste0(datadir,"/DTall_5.RData"))

# Define climbing
DTall[, climbingEE := wagechange_wave > 0.05]
DTall[, climbingEUE := wagechangeEUE_wave > 0.05]

# Probability of climbing: EE
EE <- DTall[EE_wave == TRUE, .(P_climb = weighted.mean(climbingEE, wpfinwgt, na.rm = TRUE)), 
	    by = list(switchedOcc_wave, ageGrp, HSCol)]
EE <- EE[!is.na(switchedOcc_wave),]
EE <- melt(EE, id.vars = c("switchedOcc_wave", "ageGrp", "HSCol"))
EE$ageGrp <- factor(EE$ageGrp)
levels(EE$variable) <- ""
levels(EE$ageGrp) <- c("Young", "Prime", "Old")
EE$HSCol <- factor(EE$HSCol)
levels(EE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EE, aes(x = variable, y = value, fill = factor(switchedOcc_wave))) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_grid(ageGrp ~ HSCol) +
	coord_cartesian(ylim=c(0.3,.7)) +
	xlab("") +
	ylab("Probability of Earnings Gain") +
	theme_bw() +
	geom_abline(intercept = 0.5, slope=0.)+
	scale_fill_discrete(name = "Occupation Switch",
			    labels = c("Did not switch occupation", "Switched occupation")) +
	ggtitle("E-E: Probability of Wage Gain")
ggsave(filename = paste0(figuredir, "/EEladder.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/EEladder.png"), width = 10,height = 5)

# Probability of climbing: EUE
EUE <- DTall[!is.na(wagechangeEUE_wave), .(P_climb = weighted.mean(climbingEUE, wpfinwgt, na.rm = TRUE)), 
	     by = list(switchedOcc_wave, ageGrp, HSCol)]
EUE <- EUE[!is.na(switchedOcc_wave),]
EUE <- melt(EUE, id.vars = c("switchedOcc_wave", "ageGrp", "HSCol"))
EUE$ageGrp <- factor(EUE$ageGrp)
levels(EUE$ageGrp) <- c("Young", "Prime", "Old")
levels(EUE$variable) <- ""
EUE$HSCol <- factor(EUE$HSCol)
levels(EUE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EUE, aes(x = variable, y = value, fill = factor(switchedOcc_wave))) +
	geom_bar(stat = "identity", position = "dodge") +
	facet_grid(ageGrp ~ HSCol) +
	coord_cartesian(ylim=c(0.3,.7)) +
	xlab("") +
	ylab("Probability of Earnings Gain") +
	theme_bw() +
	geom_abline(intercept = 0.5, slope=0.)+
	scale_fill_discrete(name = "Occupation Switch",
			    labels = c("Did not switch occupation", "Switched occupation")) +
	ggtitle("E-U-E: Probability of Wage Gain")
ggsave(filename = paste0(figuredir, "/EUEladder.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/EUEladder.png"), width = 10,height = 5)

# Probability of occupation switch AND EE
EE <- DTall[, .(P_EEswitch = weighted.mean(switchedOcc_wave & EE_wave, wpfinwgt, na.rm = TRUE)), 
	    by = list(ageGrp, HSCol)]
EE$ageGrp <- factor(EE$ageGrp)
levels(EE$ageGrp) <- c("Young", "Prime", "Old")
EE$HSCol <- factor(EE$HSCol)
levels(EE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EE, aes(x=ageGrp, y=P_EEswitch, fill=HSCol)) +
	geom_bar(stat = "identity", position = "dodge") +
	xlab("Age Group") +
	ylab("Probability") +
	ggtitle("Pr[Occ. Switch & EE]") +
	theme_bw() +
	scale_fill_discrete(name = "Education")
ggsave(filename = paste0(figuredir, "/switchAndEE.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/switchAndEE.png"), width = 10,height = 5)

# Probability of occupation switch AND EUE
EUE <- DTall[, .(P_EUEswitch = weighted.mean(switchedOcc_wave & !is.na(wagechangeEUE_wave), wpfinwgt, na.rm = TRUE)), 
	     by = list(ageGrp, HSCol)]
EUE$ageGrp <- factor(EUE$ageGrp)
levels(EUE$ageGrp) <- c("Young", "Prime", "Old")
EUE$HSCol <- factor(EUE$HSCol)
levels(EUE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EUE, aes(x=ageGrp, y=P_EUEswitch, fill=HSCol)) +
	geom_bar(stat = "identity", position = "dodge") +
	xlab("Age Group") +
	ylab("Probability") +
	ggtitle("Pr[Occ. Switch & EUE]") +
	theme_bw() +
	scale_fill_discrete(name = "Education")
ggsave(filename = paste0(figuredir, "/switchAndEUE.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/switchAndEUE.png"), width = 10,height = 5)

# Probability of occupation switch GIVEN EE
EE <- DTall[EE_wave == TRUE, .(P_EEswitch = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), 
	    by = list(ageGrp, HSCol)]
EE$ageGrp <- factor(EE$ageGrp)
levels(EE$ageGrp) <- c("Young", "Prime", "Old")
EE$HSCol <- factor(EE$HSCol)
levels(EE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EE, aes(x=ageGrp, y=P_EEswitch, fill=HSCol)) +
	geom_bar(stat = "identity", position = "dodge") +
	xlab("Age Group") +
	ylab("Probability") +
	ggtitle("Pr[Occ. Switch | EE]") +
	theme_bw() +
	scale_fill_discrete(name = "Education")
ggsave(filename = paste0(figuredir, "/switchGivenEE.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/switchGivenEE.png"), width = 10,height = 5)

# Probability of occupation switch GIVEN EUE
EUE <- DTall[!is.na(wagechangeEUE_wave), .(P_EUEswitch = weighted.mean(switchedOcc_wave, wpfinwgt, na.rm = TRUE)), 
	     by = list(ageGrp, HSCol)]
EUE$ageGrp <- factor(EUE$ageGrp)
levels(EUE$ageGrp) <- c("Young", "Prime", "Old")
EUE$HSCol <- factor(EUE$HSCol)
levels(EUE$HSCol) <- c("Less than High School", "High School", "College")
ggplot(EUE, aes(x=ageGrp, y=P_EUEswitch, fill=HSCol)) +
	geom_bar(stat = "identity", position = "dodge") +
	xlab("Age Group") +
	ylab("Probability") +
	ggtitle("Pr[Occ. Switch | EUE]") +
	theme_bw() +
	scale_fill_discrete(name = "Education")
ggsave(filename = paste0(figuredir, "/switchGivenEUE.eps"), width = 10,height = 5)
ggsave(filename = paste0(figuredir, "/switchGivenEUE.png"), width = 10,height = 5)




