## calibration states for Guido/Victoria project


# Read in SIPP data and prepare it for analysis
library(foreign)
library(readstata13)
library(data.table)
library(zoo)
library(lubridate)
library(ggplot2)
library(mFilter)
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

#rootdir <- "G:/Research_Analyst/Eubanks/Occupation Switching/SIPP to go"
rootdir <- "~/workspace/CVW/R/"
datadir <- paste(rootdir, "InputDataDE", sep = "")
outputdir <- paste(rootdir, "Results", sep = "")
figuredir <- paste(rootdir, "Figures", sep = "")
intermed_plots = T
final_plots = F
max_wavefreq = 1 # controls whether take max over months in wave or wave-frequency observation
midTrans = 0 # controls whether to mark w/in wave transitions in the prior wave too
setwd(rootdir)

sipp <- readRDS( paste0(outputdir,"/DTall_3.RData"))

#create job tenure:
sipp[ , estintid:=NULL] #I don't trust the job id at a monthly level
sipp[ , last.lfstat := shift(lfstat), by=id]
sipp[, newemp := lfstat == 1 & (last.lfstat >= 2 | mis == 1 | (shift(EE)==T & id==shift(id)) )]
sipp[lfstat == 1 & is.na(newemp), newemp := FALSE]
sipp[newemp == TRUE, estintid := cumsum(newemp), by = id]
sipp[lfstat>1, estintid:= -1] #ensure I don't carry over a non-employment spell
sipp[lfstat == 1, estintid := na.locf(estintid, na.rm = FALSE), by = id]
sipp[, newemp := NULL]

sipp[ lfstat==1, obsten := seq_len(.N) , by=list(id,estintid)]
sipp[ lfstat==1, maxobsten := max(obsten) , by=list(id,estintid)]
sipp[ lfstat==1 , estintend:= obsten == maxobsten]
sipp[ lfstat==1 & estintend==T, endswith:= ifelse(EE==T, 1,2)]
sipp[ lfstat==1 & estintend==T & endswith==2, endswith:= ifelse(EU==T, 2,3)]

sipp[ , ltrunc:= any(mis==1 & obsten==1), by=list(id,estintid)]
sipp[ , rtrunc:= any(mis>=maxmis-23 & obsten==1), by=list(id,estintid)]
sipp[ lfstat>=2 & (shift(EU)==T & shift(id)==id), rtrunc:= any(mis>=maxmis-23 ), by=list(id,ustintid)]

#add in left-truncated time:
sipp[ obsten==1, trunctime := ifelse(year-syear>=0 ,(year-syear)*12  ,0.  ) ]
sipp[ obsten==1, trunctime := ifelse(year-syear>=0 ,(year-syear)*12  ,0.  ) ]
sipp[ lfstat==1 & ltrunc==1 & is.finite(trunctime), etenure :=obsten + trunctime]


#OCCUPATIONAL CATEGORIES
#sipp[ occ_1d==1, occ:= 1 ] # non-routine cognitive
#sipp[ occ_1d==2, occ:= 2 ] # routine cognitive
#sipp[ occ_1d==3, occ:= 3 ] # non-routine manual
#sipp[ occ_1d>=4 & occ_1d<=6, occ:= 4 ] # routine manual

Nocc = 4
Nstats = 4
tenure_dist = array(NA,dim=c(Nocc,Nstats*3))
for( oi in seq(1,Nocc) ){
	sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, summary(maxobsten)]
	tenure_dist[oi,1] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten<=3)]
	tenure_dist[oi,2] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>3 & maxobsten<=12)]
	tenure_dist[oi,3] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>12 & maxobsten<=24)]

	for(ui in c(1,2)){
		tenure_dist[oi,1+ui*Nstats] = sipp[ occ ==oi & endswith==ui & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten<=3)]
		tenure_dist[oi,2+ui*Nstats] = sipp[ occ ==oi & endswith==ui & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>3 & maxobsten<=12)]
		tenure_dist[oi,3+ui*Nstats] = sipp[ occ ==oi & endswith==ui & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>12 & maxobsten<=24)]
		tenure_dist[oi,4+ui*Nstats] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(endswith==ui, na.rm=T)]
	}
	
}

Nocc = 4
Nstats = 4
dur_dist = array(NA,dim=c(Nocc,Nstats))
for( oi in seq(1,Nocc) ){
	sipp[ occ ==oi & (panel==1996 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, summary(max.unempdur)]
	dur_dist[oi,1] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==2 & UE==T & rtrunc==F, mean(max.unempdur<=3)]
	dur_dist[oi,2] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==2 & UE==T & rtrunc==F, mean(max.unempdur>3 & max.unempdur<=12)]
	dur_dist[oi,3] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==2 & UE==T & rtrunc==F, mean(max.unempdur>12 )]
	dur_dist[oi,4] = sipp[ occ ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat<=2         & rtrunc==F, mean(esr>2,na.rm=T)]
}

## Compute tenure distribution by 2-digit occupation to match Adams-Prassl facts
Nocc = length(sipp[ lfstat==1, ftable(soc2d)])
Nstats = 4 #this has be greater than 2 in order to solve for the number of each type in each occupation:
ii =1
for( oi in seq(11,55, by=2) ){
	tenure_dist[ii,1] = sipp[ soc2d ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten<=3)]
	tenure_dist[ii,2] = sipp[ soc2d ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>3 & maxobsten<=12)]
	tenure_dist[ii,3] = sipp[ soc2d ==oi & (panel==1996 | panel==2001 | panel==2004)& lfstat==1 & obsten==maxobsten & rtrunc==F, mean(maxobsten>12 & maxobsten<=24)]
	ii=ii+1
}


