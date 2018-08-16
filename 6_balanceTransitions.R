# January 22, 2016
# Balance transitions, knock out recalls and re-weight.
# 1) 
library(data.table)
library(zoo)
library(Hmisc)

wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
outputdir = "~/workspace/CVW/R/Results"

recall_adj = T
dur_adj = F

check_durcomp = T

setwd(wd0)
Mode <- function(x) {
	ux <- unique(x[!is.na(x)])
	ux[which.max(tabulate(match(x, ux)))]
}

# for getting the duration right
CPSunempdur <- readRDS("./InputData/CPSunempDurDist.RData")


recallRecodeJobID <- function(DF){
	#this function just uses job ID to identify recalls, and may miss especially long-term when job id can reset.  
	# DF[job> 0 , jobpos := job]
	# DF[lfstat==1 | UE==T , nextjob := shift(jobpos,1,type="lead"), by=id]
	# DF[ is.finite(ustintid) & (lfstat>=2 | EU==T) , nextjob := Mode(nextjob), by=list(id,ustintid)]
	# # will convert all recall stints as lfstat == NA_integer_
	# DF[EU==T , recalled := (jobpos == nextjob) ]
	# DF[is.finite(ustintid) &(EU==T | lfstat>=2 ), recalled:=any(recalled) , by=list(id,ustintid)]
	# DF[EU==T & recalled==T, EU:=NA]
	# DF[UE==T & recalled==T, UE:=NA]
	# DF[lfstat>=2 & recalled==T, lfstat:=NA]
	# DF[ , c("jobpos","nextjob"):=NULL]
	
	# do it at wave level:
	DFseam <- DF[seam==T, ]
	DFseam[job_wave> 0 , jobpos_wave := job_wave]
	DFseam[ , next.jobpos_wave := shift(jobpos_wave, type="lead"),by=id]
	DFseam[ ustintid_wave>0 & UE_wave!=T, next.jobpos_wave:=NA]
	DFseam[ ustintid_wave>0 , next.jobpos_wave:=Mode(next.jobpos_wave), by=list(id,ustintid_wave)]
	DFseam[ EU_wave==T & matched_EUUE_max ==T , recalled_wave := jobpos_wave== next.jobpos_wave]
	DFseam[ ustintid_wave>0, recalled_wave := any(recalled_wave,na.rm = T), by=list(id,ustintid_wave)]
	DFseam[ ustintid_wave>0 & is.na(recalled_wave), recalled_wave := F]

	DFseam <- subset(DFseam, select= c("recalled_wave","id","wave"))
	
	DF<- merge(DF,DFseam,by=c("id","wave"),all.x=T)
	DF[is.finite(ustintid) &(EU==T | lfstat>=2 ), recalled:=any(recalled_wave) , by=list(id,ustintid)]
	DF[EU==T & recalled==T, EU:=NA]
	DF[UE==T & recalled==T, UE:=NA]
	DF[lfstat>=2 & recalled==T, lfstat:=NA]
	
	DF[EU_wave==T & recalled_wave==T, EU_wave:=NA]
	DF[UE_wave==T & recalled_wave==T, UE_wave:=NA]
	DF[lfstat_wave>=2 & recalled_wave==T, lfstat_wave:=NA]
	#return(DF)
}

recallRecodeShorTerm <- function(DF){
	#this function closely follows Fujita & Moscarini for identification of short-term recalls
	DF[job> 0 , jobpos := job]
	DF[ lfstat==1 | UE==T, nextjob := shift(jobpos,1,type="lead"), by=id]
	DF[ is.finite(ustintid) & (lfstat>=2 | EU==T) , nextjob := Mode(nextjob), by=list(id,ustintid)]
	# will convert all recall stints as lfstat == NA_integer_
	DF[ , ENEnoseam := lfstat >=2 & ( shift(lfstat,1,type="lead")==1 & shift(lfstat,1,type="lag")==1
									  |  (sum(lfstat ==1, na.rm=T)==2 & (shift(lfstat,1,type="lead")==2 |shift(lfstat)==2)) )
	   , by=list(id,wave)]
	DF[ is.na(ENEnoseam)==T, ENEnoseam:=F ]
	DF[ , recalled := as.numeric( ENEnoseam==T & ( (shift(jobpos)==shift(jobpos,1,type="lead") & max.unempdur==1) |
												   	( max.unempdur==2 & 
												   	  	(shift(jobpos,2,type="lag")==shift(jobpos,1,type="lead")|
												   	  	shift(jobpos,1,type="lag")==shift(jobpos,2,type="lead")) ) )),  
	   by = list( id,wave ) ]
	#impute recalls if there is a seam
	DF[ , ENEwseam := ( max.unempdur<=2 & lfstat ==2 ) &
	   	(wave != shift(wave) | wave != shift(wave,1,type="lead")) ,by=id]
	DF[is.na(ENEwseam), ENEwseam:=F]
	predrecall <- glm( recalled ~ wagechangeEUE + I(wagechangeEUE^2) + recIndic + factor(esr) + union + 
					   	factor(HSCol) + factor(Young)  + I(wagechangeEUE>(-.5)) + I(wagechangeEUE<.5 & wagechangeEUE>(-.5))
					   + switchedOcc + switchedAddress, 
					   na.action = na.exclude, data= subset(DF,ENEnoseam),family=binomial(link="probit"))
	DF[ (ENEwseam | ENEnoseam), 
		fitted_recall := predict(predrecall, type= "response", newdata=subset(DF,(ENEwseam==T | ENEnoseam==T) ))]
	DF[ ENEwseam==T & !is.na(fitted_recall), recalled:= fitted_recall ]
	DF[ , recalled:= ifelse( ENEwseam & (shift(jobpos,1,type="lead") == shift(jobpos,1,type="lag"))
							 , 1,recalled ),  by= id]
	DF[ , recalled:= ifelse( ENEwseam & max.unempdur==2 & 
							 	((shift(jobpos,1,type="lead") == shift(jobpos,2,type="lag")) | (shift(jobpos,2,type="lead") == shift(jobpos,1,type="lag")))
							 , 1,recalled ),  by= id]
	DF[is.na(recalled), recalled := 0.]
	#get the remaining recalled:
	DF[EU==T , recalled := ifelse(jobpos == nextjob,1,recalled) ]
	DF[(lfstat>=2 | EU==T) & !is.na(ustintid) & ustintid>0, recalled := max(recalled,na.rm=T), by=list(id,ustintid)]
	DF[ , EUrecalled:= ( recalled > 0.75 & EU==T), by=id]
	DF[ , UErecalled:= ( recalled > 0.75 & UE==T), by=id]
	
	DF[ recalled > 0.75, EU:= F, by=id]
	DF[ recalled > 0.75, UE:= F, by=id]
	DF[ recalled > 0.75 & lfstat>=2, lfstat := 0L, by=id]
	DF[ , c("jobpos","ENEwseam","ENEnoseam") := NULL]
}

DTall <- readRDS(paste0(datadir,"/DTall_5.RData"))

if(with(DTall,exists("matched_EUUE_max"))==F){
	DTall[ , matched_EUUE_max:= any(matched_EUUE,na.rm=F), by=list(id,wave)]
}
# sum weights for UE, EU, and EE
UEreadweight <- DTall[UE==T & !is.na(wagechange_month) & is.finite(ustintid), sum(wpfinwgt, na.rm = TRUE)]
EUreadweight <- DTall[EU==T & !is.na(wagechange_month) & is.finite(ustintid), sum(wpfinwgt, na.rm = TRUE)]
EEreadweight <- DTall[EE==T & !is.na(wagechange_month), sum(wpfinwgt, na.rm = TRUE)]

UEreadweight_wave <- DTall[UE_wave==T & !is.na(wagechange_wave) & is.finite(ustintid_wave), sum(wpfinwgt, na.rm = TRUE)]
EUreadweight_wave <- DTall[EU_wave==T & !is.na(wagechange_wave) & is.finite(ustintid_wave), sum(wpfinwgt, na.rm = TRUE)]
EEreadweight_wave <- DTall[EE_wave==T & !is.na(wagechange_wave), sum(wpfinwgt, na.rm = TRUE)]

if( recall_adj == T){
	
	DFseam <- DTall[seam==T, ]
	DFseam[job_wave> 0 , jobpos_wave := job_wave]
	DFseam[ , next.jobpos_wave := shift(jobpos_wave, type="lead"),by=id]
	DFseam[ ustintid_wave>0 & UE_wave!=T, next.jobpos_wave:=NA]
	DFseam[ ustintid_wave>0 , next.jobpos_wave:=Mode(next.jobpos_wave), by=list(id,ustintid_wave)]
	DFseam[ EU_wave==T & matched_EUUE_max ==T , recalled_wave := jobpos_wave== next.jobpos_wave]
	DFseam[ ustintid_wave>0, recalled_wave := any(recalled_wave,na.rm = T), by=list(id,ustintid_wave)]
	DFseam[ ustintid_wave>0 & is.na(recalled_wave), recalled_wave := F]
		
	DFseam <- subset(DFseam, select= c("recalled_wave","id","wave"))
	
	DTall<- merge(DTall,DFseam,by=c("id","wave"),all.x=T)
	#DTall <- DTall[DFseam,on=c("id","wave")]
	DTall[is.finite(ustintid) &(EU==T | lfstat>=2 ), recalled:=any(recalled_wave) , by=list(id,ustintid)]
	DTall[EU==T & recalled==T, EU:=NA]
	DTall[UE==T & recalled==T, UE:=NA]
	DTall[lfstat>=2 & recalled==T, lfstat:=NA]
	
	DTall[EU_wave==T & recalled_wave==T, EU_wave:=NA]
	DTall[UE_wave==T & recalled_wave==T, UE_wave:=NA]
	DTall[lfstat_wave>=2 & recalled_wave==T, lfstat_wave:=NA]
	
	UEnorecallweight <- DTall[UE==T & !is.na(wagechange_month) & is.finite(ustintid), sum(wpfinwgt, na.rm = TRUE)]
	EUnorecallweight <- DTall[EU==T & !is.na(wagechange_month) & is.finite(ustintid), sum(wpfinwgt, na.rm = TRUE)]
	EEnorecallweight <- DTall[EE==T & !is.na(wagechange_month), sum(wpfinwgt, na.rm = TRUE)]
	
	UEnorecallweight_wave <- DTall[UE_wave==T & !is.na(wagechange_wave) & is.finite(ustintid_wave), sum(wpfinwgt, na.rm = TRUE)]
	EUnorecallweight_wave <- DTall[EU_wave==T & !is.na(wagechange_wave) & is.finite(ustintid_wave), sum(wpfinwgt, na.rm = TRUE)]
	EEnorecallweight_wave <- DTall[EE_wave==T & !is.na(wagechange_wave), sum(wpfinwgt, na.rm = TRUE)]
	
	#recall rate:
	1.-UEnorecallweight/UEreadweight
	1.-UEnorecallweight_wave/UEreadweight_wave
	
}

#now subset everyone:
DTall[wagechange_notransbad==F & wagechange_wave_low==F & wagechange_wave_high==F & wagechange_wave_jcbad==F &  wagechange_wave_imp==F &
	   	!(EU_wave==T|UE_wave==T|EE_wave==T) & 
	   	lfstat_wave==1 & next.lfstat_wave==1, stayer:= T]

DTall[wagechange_notransbad==F  & # wagechange_wave_jcbad==F & 
	   	(EU_wave==T|UE_wave==T|EE_wave==T)  , changer:= T]

DTall[changer==T, stayer:= F]
DTall[stayer ==T, changer:=F]

DTall[(changer ==T & !(midEE|midEU|midUE) )| stayer==F, changermo:=T]


DTall[ !(EU_wave==T|UE_wave==T|EE_wave==T) & lfstat_wave==1 & next.lfstat_wave==1, cleaning_wts:= sum(wpfinwgt>0,na.rm = T)/sum(stayer==T,na.rm=T)  , by=date]
DTall[  (EU_wave==T|UE_wave==T|EE_wave==T), cleaning_wts:= sum(wpfinwgt>0,na.rm = T)/sum(changer==T,na.rm=T) , by=date]
DTall[ is.na(cleaning_wts)==F, cleaning_wts:= cleaning_wts*wpfinwgt]


#save only the seams as a 'wave-change' set ---------------------
DTseam <- subset(DTall, seam==T)
 
DTseam[ , c("EE","EU","UE"):= NULL]

#construct duration composition by month
if(check_durcomp==T){
	maxdurUE_wave_t <- DTseam[ UE_wave==T & midUE==F ,  c(mean(max.unempdur_wave<=1),
													 mean(max.unempdur_wave >1 & max.unempdur_wave<=3 ),
													 mean(max.unempdur_wave >3 & max.unempdur_wave<=6 ),
													 mean(max.unempdur_wave >6 & max.unempdur_wave<=12 ),
													 mean(max.unempdur_wave >12) ),by=list(wave, panel)]
	maxdurUE_t <- DTall[UE==T, c(mean(max.unempdur <=1),
								 mean(max.unempdur  >1 & max.unempdur <=3 ),
								 mean(max.unempdur  >3 & max.unempdur <=6 ),
								 mean(max.unempdur  >6 & max.unempdur <=12 ),
								 mean(max.unempdur  >12) ),by=list(wave, panel)]
	
	names(maxdurUE_t) <- c("wave","panel","frac")
	maxdurUE_t[ , dur := c(0,1,2,3,4)]
	maxdurUE_t <- reshape(maxdurUE_t,idvar=c("wave","panel"),timevar = "dur",direction="wide")
	setkeyv(maxdurUE_t,c("panel","wave"))
	
	names(maxdurUE_wave_t) <- c("wave","panel","frac")
	maxdurUE_wave_t[ , dur := c(0,1,2,3,4)]
	maxdurUE_wave_t <- reshape(maxdurUE_wave_t,idvar=c("wave","panel"),timevar = "dur",direction="wide")
	setkeyv(maxdurUE_wave_t,c("panel","wave"))
}

##balance seams EU and UE
DTseam[ is.finite(ustintid_wave), u_matched := any(matched_EUUE_max,na.rm=F), by=list(id,ustintid_wave)]

#DTseam[ is.na(matched_EUUE_wave), matched_EUUE_wave:=F]
#DTseam[ ustintid_wave>0, EU_match:= any(UE_wave), by=list(id,ustintid_wave) ]
#DTseam[ ustintid_wave>0, UE_match:= any(EU_wave), by=list(id,ustintid_wave) ]
#DTseam[ ustintid_wave>0, matched_EUUE_wave:= UE_match & EU_match ]
#DTseam[ matched_EUUE_wave!=T & EU_wave==T, EU_nomatch:= T]
#DTseam[ is.na(EU_nomatch), EU_nomatch:= F]
#DTseam[ matched_EUUE_wave!=T & UE_wave==T, UE_nomatch:= T]
#DTseam[ is.na(UE_nomatch), UE_nomatch:= F]

# these are equivalent to find non-matches:
# DTseam[ (EU_wave ==T & midEU==F) | (UE_wave==T & midUE==F), EU_match := shift(UE_wave,1,type = "lead")==T & is.finite(shift(wagechange_wave,1,type = "lead")), by=id]
# DTseam[ (EU_wave ==T & midEU==F) | (UE_wave==T & midUE==F), UE_match := shift(EU_wave,1,type = "lag" )==T & is.finite(shift(wagechange_wave,1,type = "lag" )), by=id]
# 
# DTseam[ , EU_nomatch := ((EU_match ==F | is.na(EU_match)) & EU_wave==T)]
# DTseam[ , UE_nomatch := ((UE_match ==F | is.na(UE_match)) & UE_wave==T)]
# DTseam[ EU_match==T, EU_nomatch:= F]
# DTseam[ UE_match==T, UE_nomatch:= F]

#DTseam[ is.finite(ustintid_wave), u_nomatch := any(EU_nomatch==T|UE_nomatch==T), by=list(id,ustintid_wave)]

DTseam[ is.finite(ustintid_wave), u_nomatch := any((EU_wave==T & u_matched==F)|UE_wave==T & u_matched==F), by=list(id,ustintid_wave)]
DTseam[ is.finite(ustintid_wave), EU_nomatch := any((EU_wave==T & matched_EUUE_wave==F)), by=list(id,ustintid_wave)]
DTseam[ is.finite(ustintid_wave), UE_nomatch := any((UE_wave==T & matched_EUUE_wave==F)), by=list(id,ustintid_wave)]
DTseam[, misRemaining := max(mis), by=list(id, panel)]
DTseam[, misRemaining := misRemaining-mis , by=id]
DTseam[ , wis:=seq(1,.N), by=id]
DTseam[, wisRemaining := max(wis), by=list(id, panel)]
DTseam[, wisRemaining := wisRemaining-wis , by=id]

DTseam[ , perwt:= mean(wpfinwgt), by=id]

EUtruenomatchrt <- DTseam[((EU_wave==T&midEU==F)) & wisRemaining > 4  & is.finite(ustintid_wave), 
						  sum(EU_nomatch*perwt,na.rm=T)/sum(perwt)] #missing counts as not matched
UEtruenomatchrt <- DTseam[ lfstat_wave==2 & (UE_wave==T&midUE==F) & wis > 4  & is.finite(ustintid_wave), 
						  sum(UE_nomatch*perwt,na.rm=T)/sum(perwt)]
Utruenomatchrt  <- DTseam[wis > 4 & wisRemaining > 4 & is.finite(ustintid_wave) , 
						  sum(u_nomatch*perwt,na.rm=T)/sum(perwt)]

DTseam[ matched_EUUE_max==F & EU_wave ==T, EU_wave := NA]
DTseam[ matched_EUUE_max==F & UE_wave ==T, UE_wave := NA]

#do some reweighting for left- and right-truncation
DTseam[ , truncweight := perwt]
#reweight entire u stint
wtsdisp <- array(0.,dim=(6))
for (wi in seq(2,4)){
	
	EUtruncnomatchrt <- DTseam[((EU_wave==T&midEU==F))& (wisRemaining < wi ) & is.finite(ustintid_wave), 
							  sum(EU_nomatch*perwt,na.rm=T)/sum(perwt)]
	wtsdisp[wi-1] <- (1.-EUtruenomatchrt)/(1.-EUtruncnomatchrt)
	DTseam[ (wisRemaining < wi ) & is.finite(ustintid_wave), truncweight := perwt*wtsdisp[wi-1]]
	
	UEtruncnomatchrt <- DTseam[(lfstat_wave==2&(UE_wave==T&midUE==F))& (wis < wi ) & is.finite(ustintid_wave), 
							   sum(UE_nomatch*perwt,na.rm=T)/sum(perwt)]
	wtsdisp[wi+2] <- (1.-UEtruenomatchrt)/(1.-UEtruncnomatchrt)
	
	DTseam[ (wis < wi ) & is.finite(ustintid_wave), truncweight := perwt*wtsdisp[wi+2]]
	
	
}

DTseam[ , cycweight := perwt]
#do some re-weighting for the cycle
for(ri in c(T,F)){
	for(si in c(T,F)){
		wt0 = DTseam[ recIndic_wave==ri & EE_wave==T & switched_wave==si, sum(cycweight)]
		wt1 = DTseam[ recIndic_wave==ri & EE_wave==T & switched_wave==si & !is.na(wagechange_wave), sum(cycweight)]
		DTseam[ recIndic_wave==ri & EE_wave==T & switched_wave==si & !is.na(wagechange_wave), cycweight := cycweight*wt0/wt1]
		wt0 = DTseam[ recIndic_wave==ri & EU_wave==T & switched_wave==si, sum(cycweight)]
		wt1 = DTseam[ recIndic_wave==ri & EU_wave==T & switched_wave==si & !is.na(wagechange_wave), sum(cycweight)]
		DTseam[ recIndic_wave==ri & EU_wave==T & switched_wave==si & !is.na(wagechange_wave), cycweight := cycweight*wt0/wt1]
		wt0 = DTseam[ recIndic_wave==ri & UE_wave==T & switched_wave==si, sum(cycweight)]
		wt1 = DTseam[ recIndic_wave==ri & UE_wave==T & switched_wave==si & !is.na(wagechange_wave), sum(cycweight)]
		DTseam[ recIndic_wave==ri & UE_wave==T & switched_wave==si & !is.na(wagechange_wave), cycweight := cycweight*wt0/wt1]
	}
}
# need to update wages for the midEU, midUE to be 1/2 weight
# DTseam[ , last.midEU:= shift(midEU, type="lag"), by=id]
# DTseam[ , last.midEE:= shift(midEE, type="lag"), by=id]
# DTseam[ , last.midUE:= shift(midUE, type="lag"), by=id]
# chngwt1 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(truncweight,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , truncweight := 0.5*truncweight]
# chngwt2 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(truncweight,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , truncweight := chngwt1/chngwt2*truncweight]
# chngwt1 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(cycweight,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , cycweight   := 0.5*cycweight]
# chngwt2 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(cycweight,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , cycweight := chngwt1/chngwt2*cycweight]
# chngwt1 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(perwt,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , perwt   := 0.5*perwt]
# chngwt2 <- DTseam[ EU_wave==T|UE_wave==T|EE_wave==T, sum(perwt,na.rm = T)]
# DTseam[ (midEU|last.midEU|midEE|last.midEE|midUE|last.midUE) , perwt := chngwt1/chngwt2*perwt]


# correct for cleaning and truncation
DTseam[ is.na(cleaning_wts)==F, cleaningtruncweight:= truncweight*cleaning_wts]
DTseam[ is.na(cleaning_wts)==T, cleaningtruncweight:= truncweight]



if( dur_adj == T){
	DTseam[ , max.unempdur:= maxunempdur_wave]
	DTseam <- merge(DTseam, CPSunempdur, by = "date", all.x = TRUE)
	#this is the failure rate
	DTseam[,SIPPmax_LT15  :=wtd.mean( (max.unempdur< 15*12/52)                         , perwt,na.rm=T)]
	DTseam[,SIPPmax_15_26 :=wtd.mean( (max.unempdur>=15*12/52 & max.unempdur<=26*12/52) , perwt,na.rm=T)]
	DTseam[,SIPPmax_GT26  :=wtd.mean( (max.unempdur> 26*12/52)                         , perwt,na.rm=T)]
	#this is the cumulative survival rate
	DTseam[,CPSsurv_LT15  :=(mean(FRM5_14  ,na.rm=T)+ mean(LT5,na.rm=T))/100]
	DTseam[,CPSsurv_15_26 :=(mean(FRM15_26 ,na.rm=T))/100]
	DTseam[,CPSsurv_GT26  :=(mean(GT26     ,na.rm=T))/100]
	#1-durwt_LT15*SIPPmax_LT15 = CPSsurv_GT26+ CPSsurv_15_26
	DTseam[(max.unempdur<15*12/52), durwt :=  (1.-CPSsurv_GT26+ CPSsurv_15_26)/SIPPmax_LT15]
	#durwt_FRM1526*SIPPmax_15_26 = CPSsurv_15_26
	DTseam[(max.unempdur>=15*12/52 & max.unempdur<=26*12/52), durwt :=  CPSsurv_15_26/SIPPmax_15_26]
	#durwt_GT26*SIPPmax_GT26 = CPSsurv_GT26
	DTseam[(max.unempdur>26*12/52), durwt :=  CPSsurv_GT26/SIPPmax_GT26]
	
	DTseam[ (max.unempdur<15*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	DTseam[ (max.unempdur>=15*12/52 & max.unempdur<=26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	DTseam[ (max.unempdur>26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	
	DTseam[ , c("SIPPmax_LT15","SIPPmax_15_26","SIPPmax_GT26","CPSsurv_LT15","CPSsurv_15_26","CPSsurv_GT26")
	             := NULL]
	#do not increase the incidence of unemployment
	balwtEU <- DTseam[EU_wave==T, sum(waveweight, na.rm=T)]
	balwtUE <- DTseam[UE_wave==T, sum(waveweight, na.rm=T)]
	
	DTseam[EU_wave ==T | UE_wave ==T, waveweight:= waveweight*durwt] 
	DTseam[EU_wave==T, waveweight:= waveweight/sum(waveweight,na.rm=T)*balwtEU ]
	DTseam[UE_wave==T, waveweight:= waveweight/sum(waveweight,na.rm=T)*balwtUE ]
	
	DTseam[,durwt :=NULL]

	# scale weights back to original total
	wtscale <- DTseam[, sum(waveweight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE)]
	DTseam[, waveweight := waveweight/wtscale]
}

DTseam[ !(is.finite(EE_wave) & is.finite(EU_wave)&is.finite(UE_wave)), changer:=NA]

DTseam <- subset(DTseam, is.finite(wpfinwgt) & is.finite(wagechange_wave))
saveRDS(DTseam, paste0(outputdir,"/DTseam.RData"))



# create data set with defined wage changes only ------------------------------
wagechanges <- DTall[!is.infinite(wagechange_month) & !is.na(wagechange_month),]
setkey(wagechanges, id, date)


# keep only count EEs and balanced EUs and UEs
wagechanges[, balancedEU := EU & shift(UE, 1, type = "lead"), by = id]
wagechanges[, balancedUE := UE & shift(EU, 1, type = "lag"), by = id]
wagechanges <- wagechanges[EE | balancedEU | balancedUE,]


wagechanges[, max.unempdur:= ifelse(UE==T & is.na(max.unempdur), shift(max.unempdur), max.unempdur) , by = id]
wagechanges[, max.unempdur:= ifelse(UE==T & max.unempdur<shift(max.unempdur), shift(max.unempdur), max.unempdur) , by = id]
wagechanges[, max.unempdur:= ifelse(EU==T & max.unempdur<shift(max.unempdur,1,type="lead"), shift(max.unempdur,1,type="lead"), max.unempdur) , by = id]
wagechanges[, ustintid := ifelse( EU==T & is.na(ustintid), shift(ustintid,1,type="lead"), ustintid ), by=id]
wagechanges[, ustintid := ifelse( UE==T & is.na(ustintid), shift(ustintid,1,type="lag"), ustintid ), by=id]

# set HSCol and Young to max over panel to ensure balance
wagechanges[, HSCol := max(HSCol), by = id]
wagechanges[, Young := as.logical(max(Young)), by = id]

# create new person weights
wagechanges[, perwt := mean(wpfinwgt, na.rm = TRUE), by = id]
DTall[, perwt := mean(wpfinwgt, na.rm = TRUE), by = id]


# re-weight for total distribution
UEweight.balanced <- wagechanges[UE & !is.na(wagechange_month), sum(perwt, na.rm = TRUE)]
EUweight.balanced <- wagechanges[EU & !is.na(wagechange_month), sum(perwt, na.rm = TRUE)]
EEweight.balanced <- wagechanges[EE & !is.na(wagechange_month), sum(perwt, na.rm = TRUE)]

# re-inflate weights for workers who enter and leave as unemployed & divide by 2 for the whole transition
# this should take the UE, because many EU will leave sample by exit LF
wagechanges[, balanceweight := perwt]

if(dur_adj==T){
	
	# change weight among unemployed to match <15 weeks, 15-26 weeks, 27+ from CPS
	wagechanges <- merge(wagechanges, CPSunempdur, by = "date", all.x = TRUE)
	wagechanges$year <-as.integer(format(wagechanges$date, "%Y") )
	wagechanges[ UE ==T | EU==T, year := shift(year,1,type="lag")] # be sure UE is in the same year as EU
	
	#for( yi in seq(min(wagechanges$year, na.rm=T),max(wagechanges$year, na.rm=T)  )){
		#this is the failure rate
		wagechanges[,SIPPmax_LT15  :=wtd.mean( (max.unempdur< 15*12/52)                         , perwt,na.rm = T)]
		wagechanges[,SIPPmax_15_26 :=wtd.mean( (max.unempdur>=15*12/52 & max.unempdur<=26*12/52) , perwt,na.rm = T)]
		wagechanges[,SIPPmax_GT26  :=wtd.mean( (max.unempdur> 26*12/52)                         , perwt,na.rm = T)]
		#this is the cumulative survival rate
		wagechanges[,CPSsurv_LT15  :=(mean(FRM5_14  ,na.rm = T) + mean(LT5,na.rm = T))/100]
		wagechanges[,CPSsurv_15_26 :=(mean(FRM15_26 ,na.rm = T))/100]
		wagechanges[,CPSsurv_GT26  :=(mean(GT26     ,na.rm = T))/100]
		#1-durwt_LT15*SIPPmax_LT15 = CPSsurv_GT26+ CPSsurv_15_26
		wagechanges[(max.unempdur<15*12/52), durwt := (1.-CPSsurv_GT26+ CPSsurv_15_26)/SIPPmax_LT15]
		#durwt_FRM1526*SIPPmax_15_26 = CPSsurv_15_26
		wagechanges[(max.unempdur>=15*12/52 & max.unempdur<=26*12/52), durwt :=  CPSsurv_15_26/SIPPmax_15_26]
		#durwt_GT26*SIPPmax_GT26 = CPSsurv_GT26
		wagechanges[(max.unempdur>26*12/52), durwt :=  CPSsurv_GT26/SIPPmax_GT26]
	#}
	wagechanges[ (max.unempdur<15*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	wagechanges[ (max.unempdur>=15*12/52 & max.unempdur<=26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	wagechanges[ (max.unempdur>26*12/52), quantile(durwt, na.rm=T,probs = c(.1,.25,.5,.75,.9))]
	
	wagechanges[ , c("SIPPmax_LT15","SIPPmax_15_26","SIPPmax_GT26","CPSsurv_LT15","CPSsurv_15_26","CPSsurv_GT26")
				:= NULL]
	#do not increase the incidence of unemployment
	balwtEU <- sum(wagechanges$balanceweight[wagechanges$EU], na.rm=T)
	balwtUE <- sum(wagechanges$balanceweight[wagechanges$UE], na.rm=T)
	
	wagechanges$balanceweight <- wagechanges$balanceweight*wagechanges$durwt
	wagechanges$balanceweight[wagechanges$EU] <- wagechanges$balanceweight[wagechanges$EU]/
		sum(wagechanges$balanceweight[wagechanges$EU], na.rm=T)*balwtEU
	wagechanges$balanceweight[wagechanges$UE] <- wagechanges$balanceweight[wagechanges$UE]/
		sum(wagechanges$balanceweight[wagechanges$UE], na.rm=T)*balwtUE
	
	wagechanges[,c("year","durwt"):=NULL]
}


# scale weights back to original total
wtscale <- wagechanges[, sum(balanceweight, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE)]
wagechanges[, balanceweight := balanceweight/wtscale]

# create EUE balance weights that give double weight to the EU and zero to the UE
wagechanges[, balanceweightEUE := balanceweight]
wagechanges[EU==T, balanceweightEUE := balanceweight*2]
wagechanges[UE==T, balanceweightEUE := 0.]


wagechangesBalanced<-subset(wagechanges, select=c("id","date","balancedEU","balancedUE","max.unempdur","balanceweight"))

setkey(DTall,id,date)
DTall[is.finite(ustintid), balancedEU := any(UE,na.rm=T)==T & EU==T, by = list(id,ustintid)]
DTall[is.finite(ustintid), balancedUE := any(EU,na.rm=T)==T & UE==T, by = list(id,ustintid)]
DTall<- merge(DTall,wagechangesBalanced,by=c("id","date"),all.x=T)

#make the wage changes NA if not-balanced
DTall[UE==T & balancedUE.y==F, c("wagechangeEUE","wagechange_month") := NA_real_]
DTall[EU==T & balancedEU.y==F, c("wagechangeEUE","wagechange_month") := NA_real_]

DTall[UE==T, UE := balancedUE.y==T]
DTall[EU==T, EU := balancedEU.y==T]
#there are a few from x that don't have ustintid, fill these with y.

DTall[ is.finite(ustintid), max.unempdur := max(c(max.unempdur.y, max.unempdur.x), na.rm=T), by=list(id,ustintid)]
DTall[ !is.finite(balanceweight), balanceweight:= 0.]
DTall[, balanceweight := max(balanceweight, na.rm=T), by=list(id,ustintid)]

DTall[, c("max.unempdur.x","max.unempdur.y","balancedEU.x","balancedEU.y") := NULL ]

# create weights & EUE specific stuff

DTall[                 , allwt := wpfinwgt]
DTall[EU==T|UE==T|EE==T, allwt := balanceweight]
DTall[                 , allwtEUE := allwt]
DTall[EU==T            , allwtEUE := allwtEUE*2.]
DTall[UE==T            , allwtEUE := 0.]
DTall[                 , wagechangeEUE := ifelse(EU==T, wagechangeEUE,wagechange_month)]
DTall[UE==T            , wagechangeEUE := NA_real_]

wagechanges[EE==T      , balanceweightEUE:=balanceweight]
wagechanges[EU==T      , balanceweightEUE:=balanceweight*2]
wagechanges[UE==T      , balanceweightEUE:=0.]

saveRDS(DTall,paste0(outputdir,"/DTall_6.RData"))


#- Diagnostics

# check weights
UEweight <- wagechanges[UE & !is.na(wagechange_month), sum(balanceweight, na.rm = TRUE)]
EUweight <- wagechanges[EU & !is.na(wagechange_month), sum(balanceweight, na.rm = TRUE)]
EEweight <- wagechanges[EE & !is.na(wagechange_month), sum(balanceweight, na.rm = TRUE)]


# new weight ratio diff from Fallick Fleischmann UE/EE ratio
# Fallick/Fleischman CPS numbers, just for comparison
EEpop <- 1.3630541872/100
# halfway between EU and UE prob
UEpop <- 0.5*(0.9201970443 + 0.8162561576)/100
# includes transtions to and from N
SNpop <- (31.5270935961 + 0.8527093596 + 1.5384236453 + 0.8620689655 + 1.6463054187)/100
# probability conditional on being in the labor force both periods
EEprob <-EEpop/(1.-SNpop)
UEprob <-UEpop/(1.-SNpop)
round(UEweight/EEweight, 10) - round(UEprob/EEprob, 10)

# new weights sum to original weights
wagechanges[, sum(wpfinwgt, na.rm = TRUE)] == wagechanges[, sum(balanceweight, na.rm = TRUE)]

# EU and UE balance
wagechanges[UE & switchedOcc, sum(balanceweight, na.rm = TRUE)] == 
 	wagechanges[EU & switchedOcc, sum(balanceweight, na.rm = TRUE)]

# store balanced data
saveRDS(wagechanges, paste0(outputdir,"/balancedwagechanges.RData"))
#rm(list=ls())