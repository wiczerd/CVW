# July 22, 2015
# Functions for occupation switching project
library(dplyr)
library(zoo)

# Guidelines --------------------------------------------------------------

# - refer to incoming dataset as "df"
# - refer to return dataset as "result"
# - explicitly return result
# - make no assumptions about grouping and arrangement of df
# - df is always assumed to contain date, id, and wpfinwgt
# - function name describes an action
# - give brief description of function before each function
# - state pre-condition and post-condition before each function
# - ungroup dataset before returning result
# - "R" in comments refers to survey respondent

# Processing functions ----------------------------------------------------

# Generate recession indicator for individual month and wave
# Precondition:         - df contains wave
# Postcondition:        - new variable recIndic
#                       - new variable waveRec
genRec <- function(df) {
        rec_dates <- as.Date(c("2001-02-01", "2001-12-01","2007-11-01", "2009-07-01"))
        result <- df %>%
                mutate(recIndic = (date > rec_dates[1] & date < rec_dates[2]) |
                               (date > rec_dates[3] & date < rec_dates[4]))
        result <- group_by(result, id, wave) %>%
                mutate(waveRec = as.logical(max(recIndic))) %>%
                ungroup
        return(result)
}

# Generate labor force status variable from employment status recode in the SIPP
# Precondition:         - df contains variable esr
# Postcondition:        - new variable lfStat {1: employed, 2: unemployed, 3: NILF}
genLFStat <- function(df) {
	result <- df %>%
		mutate(lfStat = ifelse(esr >= 1 & esr <= 5, 1, 0)) %>%
		mutate(lfStat = ifelse(esr == 6 | esr == 7, 2, lfStat)) %>%
		mutate(lfStat = ifelse(esr == 8, 3, lfStat))
	return(result)
}

# Impose conditions on occupation, job, and industry codes
# Precondition:         - df contains variables lfStat, occ, job, and ind23
# Postcondition:        - occ is replaced with last occ while employed if R
#                       is unemployed or NILF
#                       - NA job code is replaced with 0
#                       - job is replaced with 0 if R is unemployed or NILF
#                       - ind23 is replaced with last ind23 while employed if
#                       R is unemployed or NILF
fixOccCode <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(id,date) %>%
		mutate(occ = as.integer(ifelse(lfStat == 2 | lfStat == 3, NA, occ))) %>%
		mutate(occ = na.locf(occ, na.rm = FALSE, fromLast = TRUE)) %>%
		mutate(job = as.integer(ifelse(is.na(job), 0, job))) %>%
		mutate(job = as.integer(ifelse(lfStat == 2 | lfStat == 3, 0, job))) %>%
		ungroup
	return(result)
}

# Generate dummy variables to indicate labor flows
# Precondition:		- df contains job, occ, ind23, and lfStat variables
# Postcondition:        - new variable switchedJob 
#                       - new variable switchedOcc
#                       - new variable EE
#                       - new variable UE
genFlowDummies <- function(df) {
	df <- df %>%
		group_by(id) %>%
		arrange(id,date)
	df$switchedJob = (df$job != lead(df$job) & !(lead(df$job) == lag(df$job))  &
				#	!(lead(df$job) == lag(df$job,n=2)) &
				#	!(lead(df$job) == lag(df$job,n=3)) &
				#	!(lead(df$job) == lag(df$job,n=4)) ) & 
					lead(df$id) == df$id )
	df$switchedJob = ifelse(lead(df$id) == df$id, df$switchedJob, NA)
	#all of the unemployment-involving transitions should stay
	df$switchedJob = ifelse( df$job != lead(df$job) & (df$job==0 | lead(df$job==0)), TRUE, df$switchedJob )
	df$switchedOcc = ((df$occ != lead(df$occ)) &
				!(lead(df$occ) == lag(df$occ)) &
				df$switchedJob)
	df$switchedInd = (df$ind23 != lead(df$ind23)) & df$switchedJob
	df$EE = df$lfStat == 1 & lead(df$lfStat) == 1 & df$switchedJob
	df$EE = ifelse(lead(df$id) == df$id, df$EE, NA)
	df$EU = df$lfStat == 1 & lead(df$lfStat) == 2 & df$switchedJob
	df$EU = ifelse(lead(df$id) == df$id, df$EU, NA)
	df$UE = df$lfStat == 2 & lead(df$lfStat) == 1 & df$switchedJob
	df$UE = ifelse(lead(df$id) == df$id, df$UE, NA)
	df <- ungroup(df)
	result <- df
	return(result)
}

# Generate variables to store previous occupation and next occupation
# Precondition:		- df contains switchedJob, job, and occ variables
# Postcondition:	- new variable nextOcc
#			- new variable lastOcc
nextLastOcc <- function(df){
	df <- df %>%
		group_by(id) %>%
		arrange(id,date) %>%
		mutate(leadocc = lead(occ)) %>%
		mutate(leademp = (lead(job) != 0) ) %>%
		mutate(nextOcc = as.integer(ifelse( switchedJob & leademp, leadocc, NA_integer_)) ) %>%
		mutate(nextOcc = na.locf(nextOcc, na.rm = FALSE, fromLast=TRUE)) %>%
		select(-leadocc,-leademp) %>%
		mutate(lagocc = lag(occ)) %>%
		mutate(lagemp = (lag(job)!=0) ) %>%
		mutate(lastOcc = as.integer(ifelse( switchedJob & lagemp, lagocc, NA_integer_))) %>%
		mutate(lastOcc = na.locf(lastOcc, na.rm = FALSE) ) %>%
		select(-lagocc,-lagemp)
	#df$nextOcc = ifelse( df$switchedJob & lead(df$job) != 0 & df$id == lead(df$id), lead(df$occ), NA) 
	#df$nextOcc = na.locf(df$nextOcc, na.rm = FALSE, fromLast=TRUE)
	#df$lastOcc = ifelse( df$switchedJob & df$job != 0  & df$id == lag(df$id), lag(df$occ), NA)
	#df$lastOcc = na.locf(df$lastOcc, na.rm = FALSE) 
	df <- ungroup(df)
	return(df)
}

# Remove observations that meet exclusion criteria
# Precondition:         - df contains age and lfStat variables
# Postcondition:        - age is greater than 16, lfStat is employed or unemployed
#                       with no NA observations
sampleSelect <- function(df) {
        result <- df %>%
                filter(age >= 18 & age <= 65) %>%
                filter(!is.na(esr)) 
        return(result)
}

# drops earnings below 40 hrs per month @min wage = 6.55*40 and topcoded earnings
cleanEarn <- function(df){
	df$nomearn <- df$earnm
	df$earnm <- df$earnm/df$PCEPI*100
	df$badearn <- with(df, abs(log(lead(earnm)/earnm))>2. & abs(log(lead(earnm)/lag(earnm)))<.1 )
	df$badearn <- ifelse(df$id !=lead(df$id) | df$UE | df$EU, F, df$badearn)
	df$earnm <- ifelse(df$badearn, NA_real_ ,df$earnm)
# 	minearn = 6.55*40
# 	maxnomearn = 12500 #come back to this to better take out top coded (follow CEPR lit, should not be more than 12,500 nominal)
# 	df$earnm <- ifelse(df$earnm < minearn, as.numeric(NA), df$earnm)
# 	df$lfStat <- ifelse(df$earnm < minearn, as.integer(NA), df$lfStat)
# 	df$job <- ifelse(df$earnm < minearn, as.integer(NA), df$job)
	#df$earnm <- ifelse(df$nomearnm > maxnomearn, as.numeric(NA), df$earnm)
	#df$lfStat <- ifelse(df$nomearnm > maxnomearn, as.integer(NA), df$lfStat)
	#df$job <- ifelse(df$nomearnm > maxnomearn, as.integer(NA), df$job)
	result <- df %>%
		select(-nomearn,-badearn) #%>%
#		filter(!is.na(earnm))
	return(result)
}
# TEST FUNCTION
# Remove people who are out of the labor force at any point in the panel
removeNILFS <- function(df) {
	result <- df %>%
		group_by(id) %>%
		mutate(everNILF = max(esr, na.rm = TRUE) == 8) %>%
		mutate(everEsrNA = max(is.na(esr)) == 1) %>%
		filter(!everNILF) %>%
		select(-everNILF) %>%
		filter(!everEsrNA) %>%
		select(-everEsrNA)
	ungroup
	return(result)
}

genUnempDuration <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(id,date) 
			# generate unique id for each period respondent enters unemployment
	nspell = sum(result$EU,na.rm=T)
	result$spellID = as.integer(ifelse( lag(result$EU) & result$id == lag(result$id), 1:nspell, NA))
	result$spellID = as.integer(ifelse( result$lfStat==1, 0, result$spellID))
	# carryforward unique id
	result$spellID <- na.locf(result$spellID, na.rm = FALSE)
	# in each spell, calculate cumulative sum of unemployed
	result <- result %>%
		group_by(id, spellID) %>%
		arrange(id,date) %>%
		mutate(poolUnempDur = ifelse(lfStat == 2, as.numeric(cumsum(lfStat == 2)), 0.)) %>%
		mutate(unempDur = ifelse(lfStat == 2 , as.numeric(max(poolUnempDur)), 0.))
		# THIS SEEMS NOT TO WORK: fill it so that this is in the period prior to a spell
	result$unempDur <- ifelse(result$lfStat == 1 & lead(result$lfStat)==2, lead(result$unempDur), result$unempDur)
	result$unempDur <- ifelse(result$lfStat == 1 & !(lead(result$lfStat)==2), 0, result$unempDur)
	result <- result %>%
		select(-spellID) %>%
		ungroup
	return(result)
}


# Wage functions ----------------------------------------------------------

# Generate regressors to use in wage regression
# Precondition:         - df contains PCEPI, earnm, educ, age, and race variables
# Postcondition:        - new variable nomEarnm: nominal earnings
#                       - new variable logEarnm
#                       - new variable yearsSchool
#                       - new variable experience
#                       - new variable black
#                       - new variable hispanic
#                       - new variable year
genRegressors <- function(df) {
	result<- df %>%
		mutate(logEarnm = log(earnm + sqrt(earnm^2 + 1))) %>%
		mutate(yearsSchool = as.integer(ifelse(educ == 1, 9, NA)),
			   yearsSchool = as.integer(ifelse(educ == 2, 12, yearsSchool)),
			   yearsSchool = as.integer(ifelse(educ == 3, 14, yearsSchool)),
			   yearsSchool = as.integer(ifelse(educ == 4, 16, yearsSchool)),
			   yearsSchool = as.integer(ifelse(educ == 5, 18, yearsSchool))) %>%
		mutate(experience = age - yearsSchool) %>%
		mutate(black = (race == 2)) %>%
		mutate(hispanic = (race == 3)) %>%
		mutate(year = as.numeric(format(date, "%Y")))
	return(result)
}

# Calculate wage variable to use for analysis
# Precondition:         - df contains logEarnm, experience, educ, female, black, 
#                       hispanic, and soc2d variables
# Postcondition:        - new variable useWage: residual from wage regression or log earm
calculateUseWage <- function(df, const = 0) {
        if(useRegResid) {
                model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
                                    female + black + hispanic, data = df,
                            na.action = na.exclude, weights = wpfinwgt)
                useWage <- residuals(model) + const
                result <- data.frame(df, useWage)
        } else {
                result <- df %>%
                        mutate(useWage = logEarnm) 
        }
        return(result)
}

# Calculate residual wage within occupation (occ)
# Precondition:         - df contains logEarnm, experience, educ, female, black, 
#                       hispanic, and soc2d variables
# Postcondition:        - new variable occWage: residual from wage regression within
#                       occupation (occ)
calculateOccWage <- function(df, const = 0) {
        df <- group_by(df, soc2d)
        model <- lm(logEarnm ~ experience + I(experience^2) + factor(educ) + 
                            female + black + hispanic, 
                     data = df, na.action = na.exclude, weights = wpfinwgt)
        occWage <- fitted(model) + const
        result <- data.frame(df, occWage)
        result <- result %>% 
                ungroup
        return(result)        
}

# Fill down wage to replace NAs with most recent wage
# Precondition:         - df contains switchedJob, job, useWage, useWageQtr, and occWage
# Postcondition:        - new variable lastWage
#                       - new variable lastWageQtr
#                       - new variable lastOccWage
fillDownWage <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(id, date) %>%
		mutate(lastWage = ifelse(switchedJob & job != 0, lag(useWage), as.numeric(NA))) %>%
		mutate(lastWage = na.locf(lastWage, na.rm = FALSE))
	result <- result %>%
		mutate(lastWageQtr = ifelse(switchedJob & job != 0, lag(useWageQtr, 3), as.numeric(NA))) %>%
		mutate(lastWageQtr = na.locf(lastWageQtr, na.rm = FALSE))
	result <- result %>%
		group_by(id) %>%
		arrange(id, date) %>%
		mutate(lastOccWage = ifelse(switchedJob & job != 0, lag(occWage), as.numeric(NA))) %>%
		mutate(lastOccWage = na.locf(lastOccWage, na.rm = FALSE)) %>%
		ungroup
	return(result)
}

# Fill up wage to replace NAs with next observed wage
# Precondition:         - df contains switchedJob, job, useWage, useWageQtr, and occWage
# Postcondition:        - new variable nextWage
#                       - new variable nextWageQtr
#                       - new variable nextOccWage

fillUpWageDPLYR <- function(df) {
	#this crashes my machine on Linux
	result <- df %>%
		group_by(id) %>%
		arrange(id, date) %>%
	#UEs are the next wage and EUs are NA
		mutate(EmpTmrw = (EE | UE) ) %>%
		mutate(wTmrw = lead(useWage)) %>%
		mutate(occWTmrw = lead(occWage))
	result <- result %>%	
		mutate(nextWage = as.numeric(ifelse(EmpTmrw , wTmrw, NA_real_)) ) %>%
		mutate(nextWage = na.locf(nextWage, na.rm = FALSE, fromLast = TRUE))
	result <- result %>%
		mutate(nextOccWage = as.numeric(ifelse(EmpTmrw, occWTmrw, NA_real_))) %>%
		mutate(nextOccWage = na.locf(nextOccWage, na.rm = FALSE, fromLast = TRUE)) %>%
		select(-EmpTmrw,-wTmrw)
	#result <- result %>%
	#	mutate(EmpNQtr = switchedJob & lead(job,3) != 0) %>%
	#	mutate(wNQtr = lead(useWageQtr,3)) %>%
	#	mutate(nextWageQtr = ifelse( EmpNQtr, wNQtr, NA_real_)) %>%
	#	nextWageQtr = na.locf(nextWageQtr, na.rm = FALSE, fromLast = TRUE) %>%
	#	select(-wNQtr,-EmpNQtr)
	
	ungroup(result)
	return(result)
}
fillUpWage <- function(df) {
	result <- df %>%
		#		group_by(id) %>%
		arrange(id, date) 
	#UEs are the next wage and EUs are NA
	result$EmpTmrw <- (result$EE | result$UE) & lead(result$id) == result$id
	result$nextWage <- ifelse(result$EmpTmrw , lead(result$useWage), NA_real_)
	result$nextWage <- ifelse(result$id != lag(result$id) & result$lfStat ==2
							  , 0., result$nextWage) #force the fill back to respect ID barriers
	result$nextWage = na.locf(result$nextWage, na.rm = FALSE, fromLast = TRUE)
	result$nextWage <- ifelse(result$id != lag(result$id) & result$lfStat ==2
							  , NA_real_, result$nextWage)
	result$nextOccWage = ifelse(result$EmpTmrw, lead(result$occWage), NA_real_)
	result$nextOccWage <- ifelse(result$id != lag(result$id) & result$lfStat ==2
								 , 0., result$nextOccWage) #force the fill back to respect ID barriers
	result$nextOccWage = na.locf(result$nextOccWage, na.rm = FALSE, fromLast = TRUE)
	result$nextOccWage <- ifelse(result$id != lag(result$id) & result$lfStat ==2
							  , NA_real_, result$nextOccWage)
	result$EmpNQtr<- result$switchedJob & lead(result$job,3) != 0 &  lead(result$id,3) == result$id
	result$nextWageQtr = ifelse( result$EmpNQtr, lead(result$useWageQtr, 3), NA_real_)
	result$nextWageQtr = na.locf(result$nextWageQtr, na.rm = FALSE, fromLast = TRUE) 
	
	ungroup(result)
	return(result)
}

# Calculate change in wages between previous occupation and next occupation
# Precondition:         - df contains nextWage, useWage, nextOccWage, occWage, and switchedJob
# Postcondition:        - new variable wageChange
#                       - new variable wageChangeQtr
#                       - new variable occWageChange
#                       - new variable wageChange_stayer
#                       - new variable wageChange_all
# calculateWageChange <- function(df) {
#         result <- df %>%
#                 mutate(wageChange = nextWage - lag(useWage)) %>%
#                 mutate(wageChangeQtr = nextWageQtr - lag(useWageQtr)) %>%
#                 mutate(occWageChange = nextOccWage - lag(occWage)) %>%
#                 mutate(wageChange_stayer = as.numeric(ifelse(!switchedJob, nextWage - lag(useWage), NA) )) %>%
#                 mutate(wageChange_all = as.numeric(ifelse(!switchedJob, wageChange_stayer, wageChange)))
# }


calculateWageChange <- function(df) {
	df <- group_by(df,id)
	df <- arrange(df,id,date)
	tuw <- ifelse(lag(df$id)==df$id,lag(df$useWage), NA)
	tnw <- df$nextWage
	#tnw <- ifelse(!is.finite(tnw) & is.finite(lead(df$useWage)) & lead(df$lfStat)==1 & df$id == lead(df$id),lead(df$useWage),tnw)
	tuw <- ifelse(!is.finite(tuw) & is.finite(df$useWage) & df$lfStat==1,df$useWage,tuw)
	df$wageChange_EUE = ifelse( df$EU | df$EE, tnw - tuw  ,NA)
	df$wageChange = ifelse( df$EE , tnw - tuw  ,NA)
	df$wageChange = ifelse( df$EU , log(1.) - tuw ,df$wageChange)
	df$wageChange = ifelse( df$UE , tnw - log(1.) ,df$wageChange)
	
	df$wageChange_stayer = ifelse(df$job == lead(df$job) & df$job == lag(df$job) & df$id == lead(df$id) & df$id==lag(df$id)
								  , lead(df$useWage) - lag(df$useWage), NA)
	df$wageChange_stayer = ifelse(df$job == 0 & lead(df$job) == 0 & df$id==lag(df$id), 0., df$wageChange_stayer)
	df$wageChange_EUE = ifelse(!df$EU & !df$EE, df$wageChange_stayer, df$wageChange_EUE)
	df$wageChange_all = ifelse(df$job == lead(df$job) & df$job>0  & df$id == lead(df$id)
							   , df$wageChange_stayer, df$wageChange)
	df$wageChange_all = ifelse(df$EE | df$UE | df$EU, df$wageChange,df$wageChange_all)
	
	tuwQ <- lag(df$useWageQtr)
	tnwQ <- df$nextWageQtr
	df$wageChangeQtr = tnwQ - tuwQ
	df$occWageChange = ifelse(df$switchedOcc & !lag(df$switchedOcc),df$nextOccWage - lag(df$occWage),0.)
	df$occWageChange = ifelse(df$switchedOcc & lag(df$switchedOcc),df$nextOccWage - df$occWage,0.)
	ungroup(df)
	return(df)
}    
