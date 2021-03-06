# January 21, 2016
# Calculate wage changes
# 1) create wagechange, wagechange_stayer, wagechangeEUE, and wagechange_all variables
# 2) create occwagechange variable
# 3) save intermediate result, DTall_5.RData

library(data.table)
library(zoo)
library(stats)

#the minimum 3 month earnings allowable
minEarn <- 1040 #10 hours, $8/hour, 13-weeks

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

DTall <- readRDS(paste0(datadir,"/DTall_4.RData"))
#drop some un-need variables
DTall[ , nwhite := black==T|hispanic==T]
DTall <- DTall[ , c("eyear","logearnm","black","hispanic"):=NULL]
#I don't want to drop these ->  DTall <- DTall[ lfstat_wave<3, ]

setkey(DTall, id, date)


#***********************************************************************
#***********************************************************************
# compute seam-to-seam and year-to-year wage change variable ----------------------------
#***********************************************************************
#***********************************************************************

# residual wages
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]

DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
DTall[ , levwage := mean(levwage,na.rm=T)*nmo_lf1, by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , wavewage := log(levwage + (1+levwage^2)^.5) ]
DTall[ , nawavewage:= all(is.na(usewage)) ,by= list(id,wave)]
DTall[ nawavewage==T, wavewage:=NA_real_]
#drop the lowest resid wages, implies working less than $80 /month:
DTall[ lfstat_wave==1 & wavewage<log(minEarn+(1+minEarn^2)^.5), wavewage:=NA]

#put back levwage for esr ~=1
DTall[!is.na(usewage) & esr==1, esr1wg := 1/2*(exp(usewage)-exp(-usewage))]
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
DTall[ , esr1wg := mean(esr1wg,na.rm=T)*nmo_lf1, by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , empwg := log(esr1wg + (1+esr1wg^2)^.5) ]

DTall[ , c("nawavewage","esr1wg"):=NULL]


# raw wages
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,wave)]
DTall[ , waveearnm := mean(earnm,na.rm=T)*nmo_lf1, by= list(id,wave)] #if one month is missing, give it the average of the other 3
DTall[ , waverawwg := log(waveearnm + (1+waveearnm^2)^.5) ]
DTall[ , nawavewage:= all(is.na(usewage) ),by= list(id,wave)]
DTall[ nawavewage==T, waverawwg:=NA_real_]
#drop the lowest wages, implies working less than $80 /month:
DTall[ lfstat_wave==1 & waverawwg<log(minEarn+(1+minEarn^2)^.5), waverawwg:=NA]
DTall[ lfstat_wave==1, earn_imp_wave := sum(earn_imp==1,na.rm=T), by=list(id,wave)]

DTall[ , c("nawavewage"):=NULL]

#************************************************************************************
#wave frequency changes --------------------------------------------
#************************************************************************************
DTseam <- DTall[ seam==T,]

DTseam[ , next2.lfstat_wave := shift(lfstat_wave,2,type="lead"), by=id] 
DTseam[ , next3.lfstat_wave := shift(lfstat_wave,3,type="lead"), by=id] 

	
#need to add change across waves (use wavewage)

DTseam[                                    , next.empwg     := shift(empwg,1,type="lead"), by=id]
DTseam[ wave+1!=shift(wave,type = "lead")  , next.empwg     := NA]
DTseam[                                    , next.wavewage     := shift(wavewage,1,type="lead"), by=id]
DTseam[ wave+1!=shift(wave,type = "lead")  , next.wavewage     := NA]
DTseam[                                    , next2.wavewage    := shift(wavewage,2,type="lead"), by=id]
DTseam[ wave+2!=shift(wave,2,type = "lead"), next2.wavewage    := NA]
DTseam[                                    , next3.wavewage    := shift(wavewage,3,type="lead"), by=id]
DTseam[ wave+3!=shift(wave,3,type = "lead"), next3.wavewage    := NA]

DTseam[ , nextann.wavewage := levwage + shift(levwage,type="lead") + shift(levwage,2,type="lead")  , by=id]
DTseam[ wave+1!=shift(wave,type = "lead") | wave+2!=shift(wave,2,type = "lead") , nextann.wavewage := NA_real_ ]
DTseam[ , nextann.wavewage := log(nextann.wavewage + (1+nextann.wavewage^2)^.5) ]

DTseam[                                  , last.empwg    := shift(empwg,1,type="lag") , by=id]
DTseam[ wave-1!=shift(wave,type = "lag" ), last.empwg    := NA]
DTseam[                                  , last.wavewage    := shift(wavewage,1,type="lag") , by=id]
DTseam[ wave-1!=shift(wave,type = "lag" ), last.wavewage    := NA]
DTseam[                                  , last.lfstat_wave := shift(lfstat_wave,1,type="lag") , by=id]
DTseam[ wave-1!=shift(wave,type = "lag" ), last.lfstat_wave := NA]
DTseam[ next.lfstat_wave==1  &  next.wavewage<log(minEarn+(1+minEarn^2)^.5) , next.wavewage  := NA_real_]
DTseam[ next2.lfstat_wave==1 & next2.wavewage<log(minEarn+(1+minEarn^2)^.5) , next2.wavewage := NA_real_]
DTseam[ next3.lfstat_wave==1 & next3.wavewage<log(minEarn+(1+minEarn^2)^.5) , next3.wavewage := NA_real_]
DTseam[ last.lfstat_wave==1  & last.wavewage <log(minEarn+(1+minEarn^2)^.5) , last.wavewage := NA_real_]
DTseam[ , next.esr_max := shift(esr_max, type="lead")]

DTseam[ , wagechange_wave  := next.wavewage   - wavewage] #just using straight, period-wise wage changes. 
DTseam[ , wagechange2_wave := next2.wavewage  - wavewage] 
DTseam[ , wagechange3_wave := next3.wavewage  - wavewage] 

DTseam[ , lastann.wavewage:= shift(levwage) + shift(levwage,2) + shift(levwage,3), by=id]
DTseam[ wave-1!=shift(wave)|wave-2!=shift(wave,2)|wave-3!=shift(wave,3), lastann.wavewage := NA_real_ ]
DTseam[ , lastann.wavewage:= log(lastann.wavewage + (1+lastann.wavewage^2)^.5 )]

DTseam[ , wagechange_anan := nextann.wavewage - lastann.wavewage]

DTseam[ , next.wagechange_wave := shift(wagechange_wave, type="lead"),by=id]
DTseam[ , last.wagechange_wave := shift(wagechange_wave, type="lag" ),by=id]

# create wave-level EUE wage change-----------------------------------------
# find wage in period before an EU and one period after UE
DTseam[ , EU_wave_first := EU_wave==T , by=id]
DTseam[ , UE_wave_last  := UE_wave==T , by=id]

DTseam[ , last.lfstat_wave:= shift(lfstat_wave), by=id]
DTseam[last.lfstat_wave==1 & EU_wave_first == T & last.stable_emp==T, wageAtEU     := last.empwg]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wageAfterUE  := next.empwg]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wage2AfterUE := next2.wavewage]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wage3AfterUE := next3.wavewage]

DTseam[, wageAtEU := na.locf(wageAtEU, na.rm = F),by=list(ustintid_wave, id)]

DTseam[UE_wave_last == T, wagechangeEUE_wave  := wageAfterUE - wageAtEU]
DTseam[UE_wave_last == T, wagechange2EUE_wave := wage2AfterUE - wageAtEU]
DTseam[UE_wave_last == T, wagechange3EUE_wave := wage3AfterUE - wageAtEU]

DTseam[, wagechangeEUE_wave := Mode(wagechangeEUE_wave) , by=list(ustintid_wave, id)]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, wagechangeEUE_wave  := next.empwg  - last.empwg]
#DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, wagechangeEUE_wave  := next.wavewage  - last.wavewage]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, wagechange2EUE_wave := next2.wavewage - last.wavewage]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, wagechange3EUE_wave := next3.wavewage - last.wavewage]
#DTseam[ !(EU_wave|UE_wave|EE_wave), wagechangeEUE_wave := wagechange_wave]
DTseam[ !(EU_wave|UE_wave|EE_wave) , wagechangeEUE_wave := next.empwg - empwg]


DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, esrAfterUE := next.esr_max]
DTseam[ustintid_wave>0 , esrAfterUE := Mode(esrAfterUE),by = list(ustintid_wave,id)]
DTseam[EU_wave==T ,    next.esr_max := esrAfterUE]

DTseam[ , c("wageAtEU","wageAfterUE","wage2AfterUE","wage3AfterUE","EU_wave_first", "UE_wave_last"):=NULL]


#************************************************************************************
#adding the raw-earnings changes:-----------------------------------------------------
#************************************************************************************

DTseam[                                    , next.waverawwg  := shift(waverawwg,1,type="lead"), by=id]
DTseam[ wave+1!=shift(wave,type = "lead")  , next.waverawwg  := NA]
DTseam[                                    , next2.waverawwg := shift(waverawwg,2,type="lead"), by=id]
DTseam[ wave+2!=shift(wave,2,type = "lead"), next2.waverawwg := NA]
DTseam[                                    , next3.waverawwg := shift(waverawwg,3,type="lead"), by=id]
DTseam[ wave+3!=shift(wave,3,type = "lead"), next3.waverawwg := NA]
DTseam[                                    , last.waverawwg  := shift(waverawwg,1,type="lag"), by=id]
DTseam[ wave-1!=shift(wave,type = "lag" )  , last.waverawwg  := NA]
DTseam[ next.lfstat_wave==1 & next.waverawwg< log(minEarn+(1+minEarn^2)^.5)  , next.waverawwg := NA_real_]
DTseam[ next.lfstat_wave==1 & next3.waverawwg<log(minEarn+(1+minEarn^2)^.5) , next2.waverawwg := NA_real_]
DTseam[ next.lfstat_wave==1 & next2.waverawwg<log(minEarn+(1+minEarn^2)^.5) , next3.waverawwg := NA_real_]
DTseam[ last.lfstat_wave==1 & last.waverawwg< log(minEarn+(1+minEarn^2)^.5)  , last.waverawwg := NA_real_]

DTseam[ , rawwgchange_wave := next.waverawwg - waverawwg]

DTseam[ , next.rawwgchange_wave := shift(rawwgchange_wave, type="lead"),by=id]
DTseam[ , last.rawwgchange_wave := shift(rawwgchange_wave, type="lag" ),by=id]

DTseam[ , lastann.rawwg:= shift(earnm) + shift(earnm,2) + shift(earnm,3), by=id]
DTseam[ wave-1!=shift(wave)|wave-2!=shift(wave,2)|wave-3!=shift(wave,3), lastann.rawwg := NA_real_ ]
DTseam[ , lastann.rawwg:= log(lastann.rawwg + (1+lastann.rawwg^2)^.5 )]

DTseam[ , nextann.rawwg := earnm + shift(earnm,type="lead") + shift(earnm,2,type="lead")  , by=id]
DTseam[ wave+1!=shift(wave,type = "lead") | wave+2!=shift(wave,2,type = "lead") , nextann.rawwg := NA_real_ ]
DTseam[ , nextann.rawwg := log(nextann.rawwg + (1+nextann.rawwg^2)^.5) ]


DTseam[ , rawwgchange_anan := nextann.rawwg - lastann.rawwg]

# create wave-level EUE wage change-----------------------------------------
# find wage in period before an EU and one period after UE
DTseam[ , EU_wave_first := EU_wave==T , by=id]
DTseam[ , UE_wave_last  := UE_wave==T , by=id]

DTseam[last.lfstat_wave==1 & EU_wave_first == T & last.stable_emp==T, wageAtEU     := last.waverawwg]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wageAfterUE  := next.waverawwg]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wage2AfterUE := next2.waverawwg]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, wage3AfterUE := next3.waverawwg]


DTseam[, wageAtEU := na.locf(wageAtEU, na.rm = F),by=list(ustintid_wave, id)]

DTseam[UE_wave_last == T, rawwgchangeEUE_wave  := wageAfterUE  - wageAtEU]
DTseam[UE_wave_last == T, rawwgchange2EUE_wave := wage2AfterUE - wageAtEU]
DTseam[UE_wave_last == T, rawwgchange3EUE_wave := wage3AfterUE - wageAtEU]
DTseam[, rawwgchangeEUE_wave := Mode(rawwgchangeEUE_wave) , by=list(ustintid_wave, id)]
DTseam[, rawwgchange2EUE_wave:= Mode(rawwgchange2EUE_wave), by=list(ustintid_wave, id)]
DTseam[, rawwgchange3EUE_wave:= Mode(rawwgchange3EUE_wave), by=list(ustintid_wave, id)]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, rawwgchangeEUE_wave  := next.waverawwg  - last.waverawwg]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, rawwgchange2EUE_wave := next2.waverawwg - last.waverawwg]
DTseam[ EE_wave==T & last.stable_emp==T & next.stable_emp==T, rawwgchange3EUE_wave := next3.waverawwg - last.waverawwg]
DTseam[ !(EU_wave|UE_wave|EE_wave), rawwgchangeEUE_wave := rawwgchange_wave]

DTseam[ , c("wageAtEU","wageAfterUE","wage2AfterUE","wage3AfterUE","EU_wave_first", "UE_wave_last"):=NULL]


#cleaning the stayers:-----------------------------------------------------------------
# wagechange wave =NA for large gains or losses that revert

DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & (last.lfstat_wave==1) & (next.lfstat_wave==1) & lfstat_wave==T, 
	   wagechange_notransbad := (rawwgchange_wave>2|rawwgchange_wave<(-2))  & abs(next.waverawwg - last.waverawwg)<.1 ] 


minwaverawwg <- DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & lfstat_wave==1, min(waverawwg,na.rm=T)]
DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & lfstat_wave==1 & next.lfstat_wave==1 ,
	   wagechange_ENbad := (next.waverawwg<minwaverawwg)] 


DTseam[lfstat_wave>=2 & next.lfstat_wave>=2 & !(EU_wave==T|UE_wave==T|EE_wave==T) , wagechange_notransbad:= T] 

DTseam[ is.na(wagechange_notransbad)  , wagechange_notransbad :=F] 
DTseam[ is.na(wagechange_ENbad)       , wagechange_ENbad :=F] 


DTseam[ !(EU_anan|UE_anan|EE_anan) & nextann.wavewage<log(minEarn+(1+minEarn^2)^.5), wagechange_ENbad_anan := T]
DTseam[                              lastann.wavewage<log(minEarn+(1+minEarn^2)^.5), wagechange_ENbad_anan := T]

DTseam[ is.na(wagechange_ENbad_anan), wagechange_ENbad_anan :=F] 

DTseam[ !(EU_anan|UE_anan|EE_anan), lastann.notransbad := shift(wagechange_notransbad)|shift(wagechange_notransbad,2)|shift(wagechange_notransbad,3) ,by=id]
DTseam[ !(EU_anan|UE_anan|EE_anan) & is.na(lastann.notransbad) & !is.na(shift(wagechange_notransbad)), lastann.notransbad := shift(wagechange_notransbad)|shift(wagechange_notransbad,2) ,by=id]

DTseam[ , lfstat_wave3 := shift(lfstat_wave)==1 & lfstat_wave==1 & shift(lfstat_wave,type="lead")==1 & shift(lfstat_wave,2,type="lead")==1 & shift(lfstat_wave,3,type="lead")==1,by=id]
DTseam[ lfstat_wave3==T, nextann.notransbad := (wagechange_notransbad)|shift(wagechange_notransbad,type="lead")|shift(wagechange_notransbad,type="lead",2) ,by=id]
DTseam[ ,lfstat_wave3 := NULL]
DTseam[ , lfstat_wave3 := (shift(lfstat_wave)==1 & lfstat_wave==1 & shift(lfstat_wave,type="lead")==1 & shift(lfstat_wave,2,type="lead")==1) ,by=id]
DTseam[ lfstat_wave3==T & is.na(nextann.notransbad) & !is.na(wagechange_notransbad), nextann.notransbad := (wagechange_notransbad)|shift(wagechange_notransbad,type="lead") ,by=id]
DTseam[ , wagechange_notrbad_anan := nextann.notransbad | lastann.notransbad]
DTseam[ is.na(wagechange_notrbad_anan), wagechange_notrbad_anan:= F]
DTseam[ ,lfstat_wave3 := NULL]

DTseam[!(EU_wave==T|UE_wave==T|EE_wave==T) & jobchng_wave==T , wagechange_wave_jcbad := T]
DTseam[ , next.job_wave := shift(job_wave,type="lead"),by=id]
DTseam[ , last.job_wave := shift(job_wave,type="lag"),by=id]
DTseam[ lfstat_wave==1 & next.lfstat==1, no_jobchng_wave:= (job_wave==next.job_wave)&(job_wave==last.job_wave)]
# DTseam[ EE_wave==T & no_jobchng_wave==T, wagechange_wave_jcbad :=T]
DTseam[is.na(wagechange_wave_jcbad )==T , wagechange_wave_jcbad := F]
DTseam[ , c("next.job_wave","last.job_wave"):=NULL]

#wagechanges in the 2004 months with many lost observations:
# ? DTseam[ wave>=7 & panel==2004, wagechange_wave_2004bad :=T]

#lowest/highest wages out:
lowwageqtls= DTseam[ lfstat_wave==1, quantile(waverawwg, na.rm = T, probs=c(.02))]
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_low :=waverawwg<lowwageqtls[1] | next.waverawwg<lowwageqtls[1] ]
highwageqtls= DTseam[ lfstat_wave==1, quantile(waverawwg, na.rm = T, probs=c(.98))]
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T)  , wagechange_wave_high :=waverawwg>highwageqtls[1] | next.waverawwg>highwageqtls[1] ]

# take out the early attrition
DTseam[ , panelmaxmis:= max(maxmis,na.rm = T), by=panel]
DTseam[ , pctmaxmis:= maxmis/panelmaxmis]

# take out imputed earnings:
DTseam[ !(EU_wave==T|UE_wave==T|EE_wave==T) & (earn_imp_wave>=1 | shift(earn_imp_wave,type="lead")>=1), wagechange_wave_imp := T, by=id]
DTseam[ is.na(wagechange_wave_imp)==T, wagechange_wave_imp := F]


DTseam<-subset(DTseam, select = c("wagechange_wave","wagechange_wave_jcbad","wagechange_notransbad","wagechange_ENbad","wagechange_ENbad_anan","wagechange_notrbad_anan"
								  ,"wagechange_wave_low","wagechange_wave_high","wagechange_wave_imp","pctmaxmis","next.esr_max"
								  ,"wagechangeEUE_wave","next.wavewage","rawwgchangeEUE_wave","rawwgchange_wave","wagechange_anan","lastann.wavewage","nextann.wavewage",
								  "rawwgchange_anan","lastann.rawwg","id","wave"))
DTall<- merge(DTall,DTseam, by=c("id","wave"),all.x=T)


#****************************************************
#adding the occupational wage change-------------
#***************************************************

DTall[ , occwage_wave:= sum(1/2*(exp(occwage)-exp(-occwage)),na.rm=T) , by=list(id,wave)] #sum in levels
DTall[ , occwage_wave:= log(occwage_wave + (1+occwage_wave^2)^.5)] #convert back to inv hyperbolic sine
DTall[ , nawavewage:= all(is.na(occwage) ),by= list(id,wave)]
DTall[ nawavewage==T, occwage_wave:=NA_real_]

DTseam <- DTall[ seam==T]
DTseam[ , next.occwage_wave := shift(occwage_wave,type="lead"), by=id]
DTseam[ , last.occwage_wave := shift(occwage_wave,type="lag") , by=id]
DTseam[ switchedOcc_wave==T, occwagechange_wave:= next.occwage_wave-occwage_wave]
DTseam[ switchedOcc_wave==F, occwagechange_wave:= 0.]

DTseam[ , EU_wave_first := EU_wave==T , by=id]
DTseam[ , UE_wave_last  := UE_wave==T , by=id]

DTseam[ , last.lfstat_wave:= shift(lfstat_wave), by=id]
DTseam[last.lfstat_wave==1 & EU_wave_first == T & last.stable_emp==T, occwageAtEU     := last.occwage_wave]
DTseam[next.lfstat_wave==1 & UE_wave_last  == T & next.stable_emp==T, occwageAfterUE  := next.occwage_wave]

DTseam[, occwageAtEU := na.locf(occwageAtEU, na.rm = F),by=list(ustintid_wave, id)]

DTseam[UE_wave_last == T, occwagechangeEUE_wave  := occwageAfterUE - occwageAtEU]

DTseam[, occwagechangeEUE_wave:= Mode(occwagechangeEUE_wave), by=list(ustintid_wave, id)]

DTseam[ EE_wave==T, occwagechangeEUE_wave:= occwagechange_wave]
DTseam[ switchedOcc_wave==F, occwagechangeEUE_wave:= 0.]

DTseam<-subset(DTseam, select = c("occwagechange_wave","occwagechangeEUE_wave","id","wave"))
DTall <-merge(DTall,DTseam,by=c("id","wave"),all.x=T)

saveRDS(DTall, paste0(datadir,"/DTall_5.RData"))


#************************************************************************
#Annual wage changes ---------------------------
#************************************************************************

# residual wages
DTall[!is.na(usewage) , levwage := 1/2*(exp(usewage)-exp(-usewage))]

DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,yri)]
DTall[ , annwage := mean(levwage,na.rm = T)*nmo_lf1, by=list(id,yri)]
DTall[ , annwage := log(annwage + (1+annwage^2)^0.5)]
DTall[ , naannwage:= all(is.na(usewage) ),by= list(id,yri)]
DTall[ naannwage==T, annwage:=NA_real_]
DTall[ , c("levwage","naannwage"):=NULL]

# raw wages
DTall[ , nmo_lf1 := sum(lfstat==1), by=list(id,yri)]
DTall[ , annrawwg := mean(earnm,na.rm = T)*nmo_lf1, by=list(id,yri)]
DTall[ , annrawwg := log(annrawwg + (1+annrawwg^2)^0.5)]
DTall[ , naannwage:= all(is.na(usewage) ),by= list(id,yri)]
DTall[ naannwage==T, annrawwg:=NA_real_]
DTall[ , c("naannwage","nmo_lf1"):=NULL]


#merge together transitions
sipp_ann <- readRDS(paste0(datadir,"/sipp_ann.RData"))
DTall<- merge(DTall,sipp_ann,by=c("id","yri"),all.x=T)

DTall[ , yrseam := (yri !=shift(yri,type="lead") | mis==maxmis), by=id]
DTann <- DTall[ yrseam==T,]
#need to add change across waves (use wavewage)
DTann[ , next.annwage := shift(annwage,1,type="lead"), by=id]
DTann[ , next2.annwage := shift(annwage,2,type="lead"), by=id]
DTann[ , next3.annwage := shift(annwage,3,type="lead"), by=id]
DTann[ , last.annwage := shift(annwage,1,type="lag"), by=id]

DTann[ , wagechange_ann := next.annwage - annwage]
DTann[ , wagechange2_ann := next2.annwage - annwage]
DTann[ , wagechange3_ann := next3.annwage - annwage]

DTann[ , next.wagechange_ann := shift(wagechange_ann, type="lead"),by=id]
DTann[ , last.wagechange_ann := shift(wagechange_ann, type="lag" ),by=id]

DTann[ annwage>0     , rank.annwage := rank(annwage,na.last = "keep",ties.method="average"), by=year]
DTann[ annwage>0     , rank.annwage := (rank.annwage - min(rank.annwage,na.rm = T))/(max(rank.annwage,na.rm = T) - min(rank.annwage,na.rm = T)), by=year]
DTann[ annwage>0     , pctile.annwage := round(rank.annwage*100)]
DTann[ last.annwage>0, last.rank.annwage := rank(last.annwage,na.last = "keep",ties.method="average"), by=year]
DTann[ last.annwage>0, last.rank.annwage := (last.rank.annwage - min(last.rank.annwage,na.rm = T))/(max(last.rank.annwage,na.rm = T) - min(last.rank.annwage,na.rm = T)), by=year]
DTann[ last.annwage>0, last.pctile.annwage := round(last.rank.annwage*100)]

#raw annual wage changes:
DTann[ , next.annrawwg := shift(annrawwg,1,type="lead"),by=id]
DTann[ , next2.annrawwg := shift(annrawwg,2,type="lead"),by=id]
DTann[ , next3.annrawwg := shift(annrawwg,3,type="lead"),by=id]
DTann[ , last.annrawwg := shift(annrawwg,1,type="lag"),by=id]

DTann[ , rawwgchange_ann := next.annrawwg - annrawwg]
DTann[ , rawwgchange2_ann := next2.annrawwg - annrawwg]
DTann[ , rawwgchange3_ann := next3.annrawwg - annrawwg]

DTann[ , next.rawwgchange_ann := shift(rawwgchange_ann)]

DTann[ annrawwg>0     , rank.annrawwg := rank(annrawwg,na.last = "keep",ties.method="average"), by=year]
DTann[ annrawwg>0     , rank.annrawwg := (rank.annrawwg - min(rank.annrawwg,na.rm = T))/(max(rank.annrawwg,na.rm = T) - min(rank.annrawwg,na.rm = T)), by=year]
DTann[ annrawwg>0     , pctile.annrawwg := round(rank.annrawwg*100)]
DTann[ last.annrawwg>0, last.rank.annrawwg := rank(last.annrawwg,na.last = "keep",ties.method="average"), by=year]
DTann[ last.annrawwg>0, last.rank.annrawwg := (last.rank.annrawwg - min(last.rank.annrawwg,na.rm = T))/(max(last.rank.annrawwg,na.rm = T) - min(last.rank.annrawwg,na.rm = T)), by=year]
DTann[ last.annrawwg>0, last.pctile.annrawwg := round(last.rank.annrawwg*100)]

DTann<-subset(DTann, select = c("rawwgchange_ann","wagechange_ann","pctile.annrawwg","pctile.annwage","id","yri"))
DTall<- merge(DTall,DTann,by=c("id","yri"),all.x=T)


saveRDS(DTall, paste0(datadir,"/DTall_5_ann.RData"))

#wc_wave <- DTall[seam==T & lfstat_wave==1 & next.lfstat_wave==1 & wagechange_wave_bad==F & wc_wave<0.2, .(wc_wave = weighted.mean(wagechange_wave, wpfinwgt, na.rm = TRUE)), by = list(panel, date)]
#ggplot(wc_wave, aes(date, wc_wave, color = panel, group = panel)) +
#	geom_point() +
#	geom_line() +xlab("") + ylab("mean wage change, stayers wave-frequency")

