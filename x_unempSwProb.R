#%%%%%%%
# Run LPM for switching probability on different transitions (EU,UE and EE)

library(data.table)
library(zoo)
library(stats)
library(texreg)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Figures"
setwd(wd0)


regressors <- c("age", 
				"educ", 
				"female", 
				"race", 
				"yearsschool",
				"experience",
				"black", 
				"hispanic")

sipp <- readRDS(paste0(datadir,"/DTall_3.RData"))


sipp[, switchedOcc_wave := switchedOcc_wave *100]
LPM_NBER_EUUE <- lm(switchedOcc_wave~recIndic_wave+
								max.unempdur_wave+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validEUUE==T) )

LPM_NBER_EUUE_alldisp <- lm(switchedOcc_wave~recIndic_wave*displaced+
								max.unempdur_wave+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validEUUE==T) )

LPM_NBER_EUUE_splitdisp <- lm(switchedOcc_wave~recIndic_wave*displaced_layoff+recIndic_wave*displaced_empclosed+recIndic_wave*displaced_slackbiz+
							  	max.unempdur_wave+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validEUUE==T) )

LPM_NBER_UE_alldisp <-lm(switchedOcc_wave~recIndic_UE*displaced+
						   	max.unempdur_wave+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validEUUE==T) )

LPM_NBER_EU_alldisp <-lm(switchedOcc_wave~recIndic_EU*displaced+
   	max.unempdur_wave+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validEUUE==T) )

LPM_NBER_EE <-lm(switchedOcc_wave~recIndic_wave+
						 	factor(ageGrp)+female+factor(HSCol), data=subset(sipp,EE_wave==T & midEE==F) )

texreg(list(LPM_NBER_EE,LPM_NBER_EUUE,LPM_NBER_EUUE_alldisp), label="LPM_all",digits=2,file=paste0(outputdir,"/LPM_all.tex"),
	 label = "Occupational switches---conditional on job transition---are pro-cyclical",
	 custom.model.names = c("EE", "EU,UE","EU,UE"),
	 custom.coef.names= c("Base","Recession","30-55",">55","Female","HS grad","Col grad","Unemp Dur","Displaced","Rec X Displaced"),
	 reorder.coef=c(2,10,9,8,seq(3,7),1))

texreg(list(LPM_NBER_EE,LPM_NBER_EUUE,LPM_NBER_EUUE_alldisp,LPM_NBER_UE_alldisp,LPM_NBER_EU_alldisp), label="LPM_all_EU",digits=3,file=paste0(outputdir,"/LPM_all_EU.tex"), 
	   label= "Whether counted at separation or finding, switches are pro-cyclical",
	   custom.model.names = c("EE", "EU,UE","EU,UE","UE","EU"),
	   custom.coef.names= c("Base","Recession","30-55",">55","Female","HS grad","Col grad","Unemp Dur","Displaced","Rec X Displaced",
	   					 "Recession @ UE","Rec UE X Disp","Recession @ EU","Rec EU X Disp"),
	   reorder.coef=c(2,11,13,10,12,14,9,8,seq(3,7),1))

texreg(list(LPM_NBER_EE,LPM_NBER_EUUE,LPM_NBER_EUUE_alldisp,LPM_NBER_EUUE_splitdisp), label="LPM_splitdisp" ,digits=3,file=paste0(outputdir,"/LPM_splitdisp"),
	   label="Switching is less cyclical for displacement due to firm-level shocks",
	   custom.model.names= c("EE","EU,UE","EU,UE","EU,UE"),
	   custom.coef.names= c(),
	   reorder.coef=c())