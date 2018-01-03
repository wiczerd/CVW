#%%%%%%%
# Run LPM for switching probability on different transitions (EU,UE and EE)

library(data.table)
library(zoo)
library(stats)
library(texreg)
library(mFilter)
library(ggplot2)

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

wd0 = "~/workspace/CVW/R"
xwalkdir = paste0(wd0, "/Crosswalks")
inputdir <- paste0(wd0, "/InputDataDE")
datadir = paste0(wd0, "/Results")
outputdir = paste0(wd0, "/Figures")
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

sipp <- subset(sipp, seam==T)

rec_dates <- read.table(textConnection(
"Peak, Trough
2001-03-01, 2001-11-01
2007-12-01, 2009-06-01"), sep=',',header=TRUE, 
colClasses=c('Date', 'Date'))


sipp[, switchedOcc_wave := switchedOcc_wave *100]
sipp[ mis==1, ltrunc := (lfstat>1)]
sipp[       , ltrunc := any(ltrunc,na.rm = F), by=list(id,ustintid)]
sipp[is.na(ltrunc), ltrunc := F]
sipp[EE_wave==T, nempdur := 0]
sipp[ ustintid>0 & ltrunc==F, nempdur := seq_len(.N)-1, by=list(id,ustintid)]
sipp[ ustintid>0 & ltrunc==F, max.nempdur := max(nempdur,na.rm = T), by=list(id,ustintid)]
sipp[ , max.nempdur_wave := Max_narm(max.nempdur), by=list(id,wave)]
sipp[ EE_wave==T, max.nempdur_wave := 0]
sipp[ EE_wave==T, max.unempdur_wave := 0]

sipp[ validEUUE==T | (EE_wave==T&midEE==F) , validTransX := T]
sipp[ (validEUUE==T & EU_wave==T & midEE==F) | (EE_wave==T&midEE==F) , validTransX_sep := T]

sipp[ , fil_unrt := `filtered.unrate$cycle`]

fil_unrt_ts <- sipp[ validTransX_sep==T, mean(fil_unrt) , by=date]
ggplot(fil_unrt_ts, aes(date,V1))+geom_line()+ylab("HP Filtered Unemployment")+  geom_rect(data = rec_dates, 
	mapping = aes(xmin = Peak,xmax=Trough ,ymin = -Inf, ymax = +Inf),inherit.aes = FALSE, fill = "red", alpha = "0.2")
ggsave(filename = paste0(outputdir,"/fil_unrt.eps"), height = 5,width=10)
ggsave(filename = paste0(outputdir,"/fil_unrt.png"), height = 5,width=10)

LPM_urt_ndur_coefsd <- array(NA,dim=c(13,2))
LPM_urt_ndur_disp_coefsd <- array(NA,dim=c(13,2))
LPM_urt_udur_coefsd <- array(NA,dim=c(13,2))
LPM_urt_udur_disp_coefsd <- array(NA,dim=c(13,2))


for( d in seq(0,12) ){
	LPM_urt_here <- summary(lm(switchedOcc_wave~fil_unrt+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & max.nempdur_wave >= d) ))
	LPM_urt_ndur_coefsd[d+1,1] <-LPM_urt_here$coefficients["fil_unrt",1]
	LPM_urt_ndur_coefsd[d+1,2] <-LPM_urt_here$coefficients["fil_unrt",2]
	LPM_urt_disp_here <- summary(lm(switchedOcc_wave~fil_unrt+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validTransX_sep==T & displaced==T & max.nempdur_wave >= d) ))
	LPM_urt_ndur_disp_coefsd[d+1,1] <-LPM_urt_disp_here$coefficients["fil_unrt",1]
	LPM_urt_ndur_disp_coefsd[d+1,2] <-LPM_urt_disp_here$coefficients["fil_unrt",2]
	
	LPM_urt_here <- summary(lm(switchedOcc_wave~fil_unrt+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validTransX_sep==T & max.unempdur_wave >= d) ))
	LPM_urt_udur_coefsd[d+1,1] <-LPM_urt_here$coefficients["fil_unrt",1]
	LPM_urt_udur_coefsd[d+1,2] <-LPM_urt_here$coefficients["fil_unrt",2]
	LPM_urt_disp_here <- summary(lm(switchedOcc_wave~fil_unrt+factor(ageGrp)+female+factor(HSCol), data=subset(sipp,validTransX_sep==T & displaced==T & max.unempdur_wave >= d) ))
	LPM_urt_udur_disp_coefsd[d+1,1] <-LPM_urt_disp_here$coefficients["fil_unrt",1]
	LPM_urt_udur_disp_coefsd[d+1,2] <-LPM_urt_disp_here$coefficients["fil_unrt",2]
}


LPM_urt_EU <- lm(switchedOcc_wave~max.nempdur_wave*fil_unrt+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & EU_wave==T) )
summary(LPM_urt_EU)
LPM_urt_disp_EU <- lm(switchedOcc_wave~fil_unrt*max.nempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & displaced==T & EU_wave==T) )
summary(LPM_urt_disp_EU)
#interact with displaced
LPM_urt_disp_EU <- lm(switchedOcc_wave~fil_unrt*max.nempdur_wave*displaced+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T& EU_wave==T) )
summary(LPM_urt_disp_EU)
#include EE ######################################################################
LPM_urt_EU_EE <- lm(switchedOcc_wave~max.nempdur_wave*fil_unrt+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T) )
summary(LPM_urt_EU_EE)
LPM_urt_disp_EU_EE <- lm(switchedOcc_wave~fil_unrt*max.nempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & displaced==T) )
summary(LPM_urt_disp_EU_EE)
#interact with displaced
LPM_urt_disp_EU_EE <- lm(switchedOcc_wave~fil_unrt*max.nempdur_wave*displaced+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T) )
summary(LPM_urt_disp_EU_EE)


# By NBER indicator: ######################################################################
LPM_NBER_EU <- lm(switchedOcc_wave~max.nempdur_wave*recIndic_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & EU_wave==T) )
summary(LPM_NBER_EU)
LPM_NBER_disp_EU <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & displaced==T & EU_wave==T) )
summary(LPM_NBER_disp_EU)
#interact with displaced
LPM_NBER_disp_EU <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave*displaced+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T& EU_wave==T) )
summary(LPM_NBER_disp_EU)
#EU and UE ##################################################################################
LPM_NBER_EUUE <- lm(switchedOcc_wave~max.nempdur_wave*recIndic_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T) )
summary(LPM_NBER_EUUE)
LPM_NBER_disp_EUUE <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T & displaced==T) )
summary(LPM_NBER_disp_EUUE)
#interact with displaced
LPM_NBER_disp_EUUE <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave*displaced+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T) )
summary(LPM_NBER_disp_EUUE)

#include EE ######################################################################
LPM_NBER_EU_EE_add <-lm(switchedOcc_wave~max.nempdur_wave + recIndic_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T) )
summary(LPM_NBER_EU_EE_add)
LPM_NBER_EU_EE <- lm(switchedOcc_wave~max.nempdur_wave*recIndic_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T) )
summary(LPM_NBER_EU_EE)
LPM_NBER_disp_EU_EE <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T & displaced==T) )
summary(LPM_NBER_disp_EU_EE)
#interact with displaced
LPM_NBER_disp_EU_EE <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave*displaced+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validTransX_sep==T) )
summary(LPM_NBER_disp_EU_EE)




LPM_NBER_EUUE <- lm(switchedOcc_wave~recIndic_wave + max.nempdur_wave+
						factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T) )
summary(LPM_NBER_EUUE)

LPM_NBER_EUUE <- lm(switchedOcc_wave~recIndic_wave*max.nempdur_wave+
								factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T) )
summary(LPM_NBER_EUUE)

LPM_NBER_disp_EUUE <- lm(switchedOcc_wave~recIndic_wave*displaced+
								max.unempdur_wave+factor(ageGrp)+female+factor(HSCol)+factor(occ), data=subset(sipp,validEUUE==T) )
summary(LPM_NBER_disp_EUUE)

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