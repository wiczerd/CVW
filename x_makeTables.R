# January 26, 2015
# Make occupation switching tables, calculate summary statistics
library(data.table)
library(zoo)
library(Hmisc)
library(reshape2)
library(xtable)

wd0 = "~/workspace/CVW/R"
xwalkdir = "~/workspace/CVW/R/Crosswalks"
setwd(wd0)

wagechanges <- readRDS("./Data/balancedwagechanges.RData")

toKeep <- c("switchedOcc",
			"Young",
			"HSCol",
			"recIndic",
			"wagechange",
			"wagechange_EUE", 
			"wagechange_all", 
			"balanceweight", 
			"EE","EU","UE")

# select toKeep columns only
wagechanges <- wagechanges[, toKeep, with = FALSE]

DTall <- readRDS("./Data/DTall_6.RData")

toKeep <- c(toKeep,"wpfinwgt","switchedJob",)

#drop the EU and UE that are not balanced


# select toKeep columns only
DTall <- DTall[, toKeep, with = FALSE]
DTall <- subset(DTall, is.finite(wpfinwgt) & is.finite(wagechange_all))

DTall[, allwt := wpfinwgt]
DTall[EU==T|UE==T|EE==T, allwt := balanceweight]
DTall[, wagechange_allEUE := ifelse(EU==T, wagechange_EUE,wagechange_all)]
DTall[UE==T, wagechange_allEUE := NA_real_]


# Full sample table-------------------------------------------------------------

tabqtls <- c(.1,.25,.5,.75,.9)

tab_fulldist <- array(0., dim=c(6,length(tabqtls)+1))
tN <- (length(tabqtls)+1)
tab_fulldist[1,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[1,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[2,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[2,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#expansion
tab_fulldist[3,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[3,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[4,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[4,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
#recession
tab_fulldist[5,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[5,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldist[6,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_all,na.rm=T,weights=allwt)]
tab_fulldist[6,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_all,na.rm=T,weights=allwt, probs=tabqtls)]

# EUE wage changes
tab_fulldistEUE <- tab_fulldist
tab_fulldistEUE[1,1]    <- DTall[!(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[1,2:tN] <- DTall[!(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldistEUE[2,1]    <- DTall[  EU==T|UE==T|EE==T, wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[2,2:tN] <- DTall[  EU==T|UE==T|EE==T, wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]
#expansion
tab_fulldistEUE[3,1]     <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[3,2:tN]  <- DTall[recIndic == F & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldistEUE[4,1]     <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[4,2:tN]  <- DTall[recIndic == F &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]
#recession
tab_fulldistEUE[5,1]     <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[5,2:tN]  <- DTall[recIndic == T & !(EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]
tab_fulldistEUE[6,1]     <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.mean(wagechange_allEUE,na.rm=T,weights=allwt)]
tab_fulldistEUE[6,2:tN]  <- DTall[recIndic == T &  (EU==T|UE==T|EE==T), wtd.quantile(wagechange_allEUE,na.rm=T,weights=allwt, probs=tabqtls)]


#output it to tables
tab_fulldist <- data.table(tab_fulldist)
names(tab_fulldist) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
#rownames(tab_fulldist) <- c("Same~Job","Chng~Job","Same~Job,~Exp","Chng~Job,~Exp","Same~Job,~Rec","Chng~Job,~Rec")
rownames(tab_fulldist) <- c("Same\ Job","Chng\ Job","Same\ Job,\ Exp","Chng\ Job,\ Exp","Same\ Job,\ Rec","Chng\ Job,\ Rec")
tab_fulldistEUE <- data.table(tab_fulldist)
names(tab_fulldistEUE) <- c("Mean","0.10","0.25","0.50","0.75","0.90")
rownames(tab_fulldistEUE) <- c("Same\ Job","Chng\ Job","Same\ Job,\ Exp","Chng\ Job,\ Exp","Same\ Job,\ Rec","Chng\ Job,\ Rec")

tab_fulldist <- xtable(tab_fulldist, label="tab:fulldist", digits=2, 
					align="l|l|lllll", caption="Distribution of earnings changes")
print(tab_fulldist,include.rownames=T, hline.after= c(0,nrow(tab_fulldist)), file="fulldist.tex")

tab_fulldistEUE <- xtable(tab_fulldistEUE, label="tab:fulldistEUE", digits=2, 
					   align="l|l|lllll", caption="Distribution of earnings changes, connecting unemployment spells")
print(tab_fulldistEUE,include.rownames=T, hline.after= c(0,nrow(tab_fulldist)), file="fulldistEUE.tex")




# Full sample var decomp -------------------------------------------------------------



dataTable <- c()

for(switchIndic in c(TRUE, FALSE)) {
	
	wagechanges <- subset(wagechangesfull, switchedOcc == switchIndic)
	
	# All job changes
	
	# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
	attach(wagechanges)
	
	qtls <- c(0.1,0.25,0.5,0.75,0.9)
	# Central tendency
	full.Mean <- wtd.mean(wagechange, balanceweight)
	full.Qs <- wtd.quantile(wagechange, balanceweight, probs = qtls)
	full.Med <- full.Qs[3]
	p50 <- full.Qs[3]
	
	# Dispersion
	full.Var <- wtd.var(wagechange, balanceweight)
	full.IQR <- full.Qs[4] - full.Qs[2]
	p90 <- full.Qs[5]
	p10 <- full.Qs[1]
	full.9010 <- p90 - p10
	
	# Skewness
	full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	full.PearsonSkew <- mean( balanceweight*(wagechange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceweight,na.rm=T)
	
	fullRow <- c("full", switchIndic, round(c(full.Mean, full.Med, full.Var, full.IQR, 
						  full.9010, full.Kelly, full.PearsonSkew), 4))
	
	# EUE changes
	
	# Central tendency
	EUE.Mean <- wtd.mean(wagechange_EUE, balanceweight)
	EUE.Qs <- wtd.quantile(wagechange_EUE, balanceweight, probs = qtls)
	EUE.Med <- EUE.Qs[3]
	p50 <- EUE.Qs[3]
	
	# Dispersion
	EUE.Var <- wtd.var(wagechange_EUE, balanceweight)
	EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
	p90 <- EUE.Qs[5]
	p10 <- EUE.Qs[1]
	EUE.9010 <- p90 - p10
	
	# Skewness
	EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	EUE.PearsonSkew <- mean(balanceweight*(wagechange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceweight ,na.rm=T)
	
	EUERow <- c("EUE", switchIndic, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
						EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
	
	# EE distribution
	
	# Central tendency
	EE.Mean <- wtd.mean(wagechange[EE], balanceweight[EE])
	EE.Qs <- wtd.quantile(wagechange[EE], balanceweight[EE], probs = qtls)
	EE.Med <- EE.Qs[3]
	p50 <- EE.Med
	
	# Dispersion
	EE.Var <- wtd.var(wagechange[EE], balanceweight[EE])
	EE.IQR <- EE.Qs[4]-EE.Qs[2]
	p90 <- EE.Qs[5]
	p10 <- EE.Qs[1]
	EE.9010 <- p90 - p10
	
	# Skewness
	EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
	EE.PearsonSkew <- mean(balanceweight[EE]*(wagechange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceweight[EE] ,na.rm=T)
	
	EERow <- c("EE", switchIndic, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
					      EE.9010, EE.Kelly, EE.PearsonSkew), 4))
	
	dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
	
	detach(wagechanges)
}

dataTable <- data.frame(dataTable)
names(dataTable) <- c("Distribution", "Switched", "Mean", "Median", "Variance", "IQR", "90-10", "Kelly", "Pearson")
levels(dataTable$Switched) <- c("No Switch", "Switch")
for(varNo in 3:9) {
	dataTable[,varNo] <- as.numeric(as.character(dataTable[,varNo]))
}

dataTableWide <- melt(dataTable, id = c("Distribution", "Switched"))
dataTableWide <- reshape(dataTableWide, timevar = c("Switched"), idvar = c("Distribution", "variable"), direction = "wide")
dataTableWide <- reshape(dataTableWide, timevar = c("Distribution"), idvar = c("variable"), direction = "wide")
rownames(dataTableWide) <- NULL

new.names = c()
for(name in names(dataTableWide)){
	new.names <- c(new.names, gsub("value.", x = name, replacement = ""))		       
}

# Young/Old ---------------------------------------------------------------

dataTable <- c()

# Young/Old loop
for(youngIndic in c(TRUE, FALSE)){
	
	wagechangesYoung <- subset(wagechangesfull, Young == youngIndic)
	
	for(switchIndic in c(TRUE, FALSE)) {
		
		wagechanges <- subset(wagechangesYoung, switchedOcc == switchIndic)
		
		# All job changes
		
		# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
		attach(wagechanges)
		
		qtls <- c(0.1,0.25,0.5,0.75,0.9)
		# Central tendency
		full.Mean <- wtd.mean(wagechange, balanceweight)
		full.Qs <- wtd.quantile(wagechange, balanceweight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wagechange, balanceweight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( balanceweight*(wagechange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceweight,na.rm=T)
		
		fullRow <- c("full", NA, switchIndic, youngIndic, round(c(full.Mean, full.Med, full.Var, full.IQR, 
									  full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes
		
		# Central tendency
		EUE.Mean <- wtd.mean(wagechange_EUE, balanceweight)
		EUE.Qs <- wtd.quantile(wagechange_EUE, balanceweight, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wagechange_EUE, balanceweight)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(balanceweight*(wagechange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceweight ,na.rm=T)
		
		EUERow <- c("EUE", NA, switchIndic, youngIndic, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
									EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution
		
		# Central tendency
		EE.Mean <- wtd.mean(wagechange[EE], balanceweight[EE])
		EE.Qs <- wtd.quantile(wagechange[EE], balanceweight[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wagechange[EE], balanceweight[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(balanceweight[EE]*(wagechange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceweight[EE] ,na.rm=T)
		
		EERow <- c("EE", NA, switchIndic, youngIndic, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
								      EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wagechanges)
	}
	
}

# HSCol -------------------------------------------------------------------

for(HSColIndic in c(0, 1, 2)){
	
	wagechangesHSCol <- subset(wagechangesfull, HSCol == HSColIndic)
	
	for(switchIndic in c(TRUE, FALSE)) {
		
		wagechanges <- subset(wagechangesHSCol, switchedOcc == switchIndic)
		
		# All job changes
		
		# avoid attaching until all changes to dataframe are finished. attached dataframe does not update!
		attach(wagechanges)
		
		qtls <- c(0.1,0.25,0.5,0.75,0.9)
		# Central tendency
		full.Mean <- wtd.mean(wagechange, balanceweight)
		full.Qs <- wtd.quantile(wagechange, balanceweight, probs = qtls)
		full.Med <- full.Qs[3]
		p50 <- full.Qs[3]
		
		# Dispersion
		full.Var <- wtd.var(wagechange, balanceweight)
		full.IQR <- full.Qs[4] - full.Qs[2]
		p90 <- full.Qs[5]
		p10 <- full.Qs[1]
		full.9010 <- p90 - p10
		
		# Skewness
		full.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		full.PearsonSkew <- mean( balanceweight*(wagechange - full.Mean)^3 ,na.rm=T)/full.Var^(3/2)/mean(balanceweight,na.rm=T)
		
		fullRow <- c("full", HSColIndic, switchIndic, NA, round(c(full.Mean, full.Med, full.Var, full.IQR, 
									  full.9010, full.Kelly, full.PearsonSkew), 4))
		
		# EUE changes
		
		# Central tendency
		EUE.Mean <- wtd.mean(wagechange_EUE, balanceweight)
		EUE.Qs <- wtd.quantile(wagechange_EUE, balanceweight, probs = qtls)
		EUE.Med <- EUE.Qs[3]
		p50 <- EUE.Qs[3]
		
		# Dispersion
		EUE.Var <- wtd.var(wagechange_EUE, balanceweight)
		EUE.IQR <- EUE.Qs[4]-EUE.Qs[2]
		p90 <- EUE.Qs[5]
		p10 <- EUE.Qs[1]
		EUE.9010 <- p90 - p10
		
		# Skewness
		EUE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EUE.PearsonSkew <- mean(balanceweight*(wagechange_EUE-EUE.Mean)^3 ,na.rm=T)/EUE.Var^(3/2)/mean(balanceweight ,na.rm=T)
		
		EUERow <- c("EUE", HSColIndic, switchIndic, NA, round(c(EUE.Mean, EUE.Med, EUE.Var, EUE.IQR, 
									EUE.9010, EUE.Kelly, EUE.PearsonSkew), 4))
		
		# EE distribution
		
		# Central tendency
		EE.Mean <- wtd.mean(wagechange[EE], balanceweight[EE])
		EE.Qs <- wtd.quantile(wagechange[EE], balanceweight[EE], probs = qtls)
		EE.Med <- EE.Qs[3]
		p50 <- EE.Med
		
		# Dispersion
		EE.Var <- wtd.var(wagechange[EE], balanceweight[EE])
		EE.IQR <- EE.Qs[4]-EE.Qs[2]
		p90 <- EE.Qs[5]
		p10 <- EE.Qs[1]
		EE.9010 <- p90 - p10
		
		# Skewness
		EE.Kelly <- ((p90 - p50) - (p50 - p10))/(p90 - p10)
		EE.PearsonSkew <- mean(balanceweight[EE]*(wagechange[EE]-EE.Mean)^3 ,na.rm=T)/EE.Var^(3/2)/mean(balanceweight[EE] ,na.rm=T)
		
		EERow <- c("EE", HSColIndic, switchIndic, NA, round(c(EE.Mean, EE.Med, EE.Var, EE.IQR, 
								      EE.9010, EE.Kelly, EE.PearsonSkew), 4))
		
		dataTable <- rbind(dataTable, fullRow, EUERow, EERow)
		
		detach(wagechanges)
	}
	
}

# Table output ------------------------------------------------------------

dataTable <- data.frame(dataTable)
names(dataTable) <- c("Distribution", "HSCol", "Switched", "Young", "Mean", "Median", "Variance", "IQR", "90-10", "Kelly", "Pearson")
levels(dataTable$Switched) <- c("No Switch", "Switch")
levels(dataTable$Young) <- c("Old", "Young")
levels(dataTable$HSCol) <- c("<HS", "HS", "College")

for(varNo in 5:11) {
	dataTable[,varNo] <- as.numeric(as.character(dataTable[,varNo]))
}

YoungOldLong <- subset(dataTable, !is.na(Young))
YoungOldLong$HSCol <- NULL

HSColLong <- subset(dataTable, !is.na(HSCol))
HSColLong$Young <- NULL

YoungOldWide <- melt(YoungOldLong, id = c("Distribution", "Switched", "Young"))
YoungOldWide <- reshape(YoungOldWide, timevar = c("Switched"), idvar = c("Distribution", "Young", "variable"), direction = "wide")
YoungOldWide <- reshape(YoungOldWide, timevar = c("Young"), idvar = c("Distribution", "variable"), direction = "wide")
YoungOldWide <- reshape(YoungOldWide, timevar = c("Distribution"), idvar = c("variable"), direction = "wide")
rownames(YoungOldWide) <- NULL

HSColWide <- melt(HSColLong, id = c("Distribution", "Switched", "HSCol"))
HSColWide <- reshape(HSColWide, timevar = "Switched", idvar = c("Distribution", "HSCol", "variable"), direction = "wide")
HSColWide <- reshape(HSColWide, timevar = "HSCol", idvar = c("Distribution", "variable"), direction = "wide")
HSColWide <- reshape(HSColWide, timevar = "Distribution", idvar = c("variable"), direction = "wide")
rownames(HSColWide) <- NULL

new.names = c()
for(name in names(YoungOldWide)){
	new.names <- c(new.names, gsub("value.", x = name, replacement = ""))		       
}
names(YoungOldWide) <- new.names

new.names = c()
for(name in names(HSColWide)){
	new.names <- c(new.names, gsub("value.", x = name, replacement = ""))		       
}
names(HSColWide) <- new.names

YoungOldAll <- YoungOldWide[1:6,c(1,2:5)]
YoungOldEUE <- YoungOldWide[1:6,c(1,6:9)]
YoungOldEE <- YoungOldWide[1:6,c(1,10:13)]
HSColAll <- HSColWide[1:6,c(1,2:7)]
HSColEUE <- HSColWide[1:6,c(1,8:13)]
HSColEE <- HSColWide[1:6,c(1,14:19)]

