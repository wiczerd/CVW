# December 21, 2015
# Run wage regressions, calculate residual wages
# 1) combine panels
# 2) create regressor variables
# 3) run regressions (usewage and occwage)
# 4) remove regressor variables
# 5) fill wages upwards to get next observed wage
# 6) save intermediate result, DTall_4.RData
library(data.table)
library(zoo)
library(stats)

#setwd("G:/Research_Analyst/Eubanks/Occupation Switching")
wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
datadir = "~/workspace/CVW/R/Results"
setwd(wd0)


regressors <- c("age", 
		"educ", 
		"female", 
		"race", 
		"yearsschool",
		"experience",
		"black", 
		"hispanic")

DTall <- readRDS(paste0(datadir,"/DTall_3.RData"))

# create regressor variables
DTall[, logearnm := log(earnm + sqrt(earnm^2 + 1))]

DTall[educ == 1, yearsschool := 9]
DTall[educ == 2, yearsschool := 12] 
DTall[educ == 3, yearsschool := 14]
DTall[educ == 4, yearsschool := 16]
DTall[educ == 5, yearsschool := 18]
DTall[ educ<=3, nocol_col:=1  ]
DTall[ educ> 3, nocol_col:=2  ]

DTall[, experience := age - yearsschool]

DTall[, black := race == 2]
DTall[, hispanic := race == 3]

# run regressions
# usewage
const <- DTall[ is.finite(logearnm), sum(logearnm * wpfinwgt, na.rm = TRUE)/sum(wpfinwgt, na.rm = TRUE)]
DTall[lfstat==1, usewage := const +  residuals(lm(logearnm ~ experience + I(experience^2) + factor(nocol_col) + 
				       	female + black + hispanic + factor(month),
				       na.action = na.exclude, weights = wpfinwgt))]

# occwage (NA observations of occ cause issues -- ~10% of employed)
DTall[lfstat==1 & is.finite(occ), occwage := fitted(lm(logearnm ~ experience + I(experience^2) + factor(nocol_col) +
			     	female + black + hispanic + factor(month), 
			     na.action = na.exclude, weights = wpfinwgt)), by = occ]

DTall[ lfstat>=2, usewage:= 0. ]
DTall[ lfstat>=2, occwage:= 0. ]
DTall[ lfstat>=2, earnm  := 0. ]

# remove regressor variables
DTall[, c("race","experience","educ","yearsschool","nocol_col") := NULL]

saveRDS(DTall, paste0(datadir,"/DTall_4.RData"))
rm(list=ls())