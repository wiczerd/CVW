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
                arrange(date) %>%
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
		arrange(date)
	df$switchedJob = (df$job != lead(df$job))
	df$switchedOcc = (df$occ != lead(df$occ)) & df$switchedJob
	df$switchedInd = (df$ind23 != lead(df$ind23)) & df$switchedJob
	df$EE = df$lfStat == 1 & lead(df$lfStat) == 1 & df$switchedJob
	df$EU = df$lfStat == 1 & lead(df$lfStat) == 2 & df$switchedJob
	df$UE = df$lfStat == 2 & lead(df$lfStat) == 1 & df$switchedJob
	ungroup(df)
	return(df)
}

# Generate variables to store previous occupation and next occupation
# Precondition:		- df contains switchedJob, job, and occ variables
# Postcondition:	- new variable nextOcc
#			- new variable lastOcc
nextLastOcc <- function(df){
	df <- df %>%
		group_by(id) %>%
		arrange(date)
	df$nextOcc = ifelse( df$switchedJob & lead(df$job) != 0, lead(df$occ), NA) 
	df$nextOcc = na.locf(df$nextOcc, na.rm = FALSE, fromLast=TRUE)
	df$lastOcc = ifelse( df$switchedJob & df$job != 0, lag(df$occ), NA)
	df$lastOcc = na.locf(df$lastOcc, na.rm = FALSE) 
	ungroup(df)
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
                #filter(esr != 8) %>%
                #filter(!is.na(occ)) %>%
                #filter(!is.na(job))
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

# Wrapper function: runs all processing functions in correct order
processWrapper <- function(df) {
	result <- df %>%
		sampleSelect %>%
		#removeNILFS %>%
		genRec %>%
		genLFStat %>%
		fixOccCode %>%
		genFlowDummies %>%
		nextLastOcc
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
	result <- df %>%
		mutate(nomEarnm = earnm) %>%
		mutate(earnm = earnm/PCEPI*100)
	#mutate(logEarnm = log(earnm)) %>%
	result<- result %>%
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
		mutate(lastWage = as.numeric(ifelse(switchedJob & job != 0, lag(useWage), NA))) %>%
		mutate(lastWage = na.locf(lastWage, na.rm = FALSE))
	result <- result %>%
		mutate(lastWageQtr = as.numeric(ifelse(switchedJob & job != 0, lag(useWageQtr, 3), NA))) %>%
		mutate(lastWageQtr = na.locf(lastWageQtr, na.rm = FALSE))
	result <- result %>%
		group_by(id) %>%
		arrange(id, date) %>%
		mutate(lastOccWage = as.numeric(ifelse(switchedJob & job != 0, lag(occWage), NA))) %>%
		mutate(lastOccWage = na.locf(lastOccWage, na.rm = FALSE)) %>%
		ungroup
	return(result)
}

# Fill up wage to replace NAs with next observed wage
# Precondition:         - df contains switchedJob, job, useWage, useWageQtr, and occWage
# Postcondition:        - new variable nextWage
#                       - new variable nextWageQtr
#                       - new variable nextOccWage
fillUpWage <- function(df) {
	result <- df %>%
		group_by(id) %>%
		arrange(id, date) 
	result$nextWage = as.numeric(ifelse(result$switchedJob & result$job != 0, lead(result$useWage), NA))
	result$nextWage = na.locf(result$nextWage, na.rm = FALSE, fromLast = TRUE)
	result$nextWageQtr = as.numeric(ifelse(result$switchedJob & result$job != 0, lead(result$useWageQtr, 3), NA) )
	result$nextWageQtr = na.locf(result$nextWageQtr, na.rm = FALSE, fromLast = TRUE) 
	result <- result %>%
		group_by(id) %>%
		arrange(date)
	result$nextOccWage = as.numeric(ifelse(result$switchedJob & result$job != 0, lead(result$occWage), NA))
	result$nextOccWage = na.locf(result$nextOccWage, na.rm = FALSE, fromLast = TRUE)
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
	tuw <- lag(df$useWage)
	tnw <- df$nextWage
	df$wageChange_EUE = ifelse( df$EU , tnw - tuw  ,NA)
	df$wageChange = ifelse( df$EE , tnw - tuw  ,NA)
	df$wageChange = ifelse( df$EU , log(1) - tuw ,df$wageChange)
	df$wageChange = ifelse( df$UE , tnw - log(1) ,df$wageChange)
	
	df$wageChange_stayer = ifelse(!df$switchedJob, lead(df$useWage) - lag(df$useWage), NA)
	df$wageChange_all = ifelse(!df$switchedJob, df$wageChange_stayer, df$wageChange)
	
	tuwQ <- lag(df$useWageQtr)
	tnwQ <- df$nextWageQtr
	df$wageChangeQtr = tnwQ - tuwQ
	df$occWageChange = df$nextOccWage - lag(df$occWage)
	return(df)
}    