# Read simulation data from CVW program


library(ggplot2)
library(reshape2)
library(data.table)
library(zoo)

setwd("~/workspace/CVW/c/")

# read grids:
Alev <- as.matrix(read.csv("Alev.csv", header=F, strip.white=T))
alev <- as.matrix(read.csv("zlev.csv", header=F, strip.white=T))
tlev <- as.matrix(read.csv("Thlev.csv", header=F, strip.white=T))
Plev <- as.matrix(read.csv("Plev.csv", header=F, strip.white=T))



whist <- fread(paste0("./whist.csv"), header = F, sep=",")
uhist <- fread(paste0("./uhist.csv"), header = F, sep=",")
jhist <- fread(paste0("./jhist.csv"), header = F, sep=",")
thhist <- fread(paste0("./thhist.csv"), header = F, sep=",")
zhist <- fread(paste0("./zhist.csv"), header = F, sep=",")
xShist <- fread(paste0("./xShist.csv"), header = F, sep=",")
xGhist <- fread(paste0("./xGhist.csv"), header = F, sep=",")
Ahist <- fread(paste0("./Ahist.csv"), header = F, sep=",")

Nsim <- dim(jhist)[1]
Tsim <- dim(jhist)[2]

whist$id <- seq(1:Nsim)
w_melt <- melt(whist, id.var = "id", value.name = "wage")
w_melt <- data.table(w_melt)
for(it in seq(1,Tsim)){ w_melt[ variable==paste0("V",it), t:= it] }
w_melt<- subset(w_melt, select = c("id","wage","t"))
simDat <- w_melt

uhist$id <- seq(1:Nsim)
u_melt <- melt(uhist, id.var = "id", value.name = "unemp")
u_melt <- data.table(u_melt)
for(it in seq(1,Tsim)){ u_melt[ variable==paste0("V",it), t:= it] }
u_melt<- subset(u_melt, select = c("id","unemp","t"))
simDat <- merge(simDat,u_melt,by=c("id","t"))

jhist$id <- seq(1:Nsim)
j_melt <- melt(jhist, id.var = "id", value.name = "occ")
j_melt <- data.table(j_melt)
for(it in seq(1,Tsim)){ j_melt[ variable==paste0("V",it), t:= it] }
j_melt<- subset(j_melt, select = c("id","occ","t"))
simDat <- merge(simDat,j_melt,by=c("id","t"))

thhist$id <- seq(1:Nsim)
th_melt <- melt(thhist, id.var = "id", value.name = "theta")
th_melt <- data.table(th_melt)
for(it in seq(1,Tsim)){ th_melt[ variable==paste0("V",it), t:= it] }
th_melt<- subset(th_melt, select = c("id","theta","t"))
simDat <- merge(simDat,th_melt,by=c("id","t"))

zhist$id <- seq(1:Nsim)
z_melt <- melt(zhist, id.var = "id", value.name = "z")
z_melt <- data.table(z_melt)
for(it in seq(1,Tsim)){ z_melt[ variable==paste0("V",it), t:= it] }
z_melt<- subset(z_melt, select = c("id","z","t"))
simDat <- merge(simDat,z_melt,by=c("id","t"))

xShist$id <- seq(1:Nsim)
xS_melt <- melt(xShist, id.var = "id", value.name = "xS")
xS_melt <- data.table(xS_melt)
for(it in seq(1,Tsim)){ xS_melt[ variable==paste0("V",it), t:= it] }
xS_melt<- subset(xS_melt, select = c("id","xS","t"))
simDat <- merge(simDat,xS_melt,by=c("id","t"))

xGhist$id <- seq(1:Nsim)
xG_melt <- melt(xGhist, id.var = "id", value.name = "xG")
xG_melt <- data.table(xG_melt)
for(it in seq(1,Tsim)){ xG_melt[ variable==paste0("V",it), t:= it] }
xG_melt<- subset(xG_melt, select = c("id","xG","t"))
simDat <- merge(simDat,xG_melt,by=c("id","t"))


Ahist$t <- seq(1:Tsim)
Ahist$Aidx <- Ahist$V1
simDat <- merge( simDat, subset( Ahist, select=c("t","Aidx") ) ,by= "t",all = T )


simDat[ , wagechng := log(wage) - log(shift(wage)), by=id]

simDat[ , occ_chng := occ != shift(occ), by=id]
