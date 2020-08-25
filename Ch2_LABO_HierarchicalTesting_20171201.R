install.packages("ggplot2")
install.packages("gridExtra")
install.packages("lme4")
install.packages("GGally")
install.packages("AICcmodavg")
install.packages("interplot")
library(gridExtra)
library(lme4)
library(GGally)
library(ggplot2)
library(AICcmodavg)
library(interplot)

data<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/2017_AcousticData_ShortenedVersion_20171130.csv")
ggpairs(data, columns=8:16)

#Labo:
#Step 1: Landscape Level
LBstep1a<-glm(LABO~Z_TreeCov250, data=data, family=poisson())
LBstep1b<-glm(LABO~Z_TreeCov500, data=data, family=poisson())
LBstep1c<-glm(LABO~Z_TreeCov1000, data=data, family=poisson())
LBstep1d<-glm(LABO~Z_TreeCov3000, data=data, family=poisson())
LBstep1null<-glm(LABO~1, data=data, family=poisson())

cand.mod.names <- c("LBstep1a", "LBstep1b", "LBstep1c", "LBstep1d", "LBstep1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 2: Patch level- base model is LBstep1null- no interactions as cover is not carried forward
LBstep2null<-glm(LABO~1, data=data, family=poisson())
LBstep2a<-glm(LABO~Z_Area, data=data, family=poisson())
LBstep2b<-glm(LABO~Z_Edge, data=data, family=poisson())
LBstep2c<-glm(LABO~Z_Area+Z_Edge, data=data, family=poisson())

cand.mod.names <- c("LBstep2a", "LBstep2b", "LBstep2c", "LBstep2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.1- Local level- base model is LBstep2b
LBstep3.1null<-glm(LABO~Z_Edge, data=data, family=poisson())
LBstep3.1a<-glm(LABO~Z_Edge+Z_BA_Total, data=data, family=poisson())
LBstep3.1b<-glm(LABO~Z_Edge+Z_DBH_Total, data=data, family=poisson())
LBstep3.1c<-glm(LABO~Z_Edge+Z_Canopy, data=data, family=poisson())
LBstep3.1d<-glm(LABO~Z_Edge+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
LBstep3.1e<-glm(LABO~Z_Edge+Z_BA_Total+Z_Canopy, data=data, family=poisson())
LBstep3.1f<-glm(LABO~Z_Edge+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
LBstep3.1g<-glm(LABO~Z_Edge+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("LBstep3.1a", "LBstep3.1b", "LBstep3.1c", "LBstep3.1d", "LBstep3.1e", "LBstep3.1f", "LBstep3.1g", "LBstep3.1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.2- Local level- base model is LBstep2c
LBstep3.2null<-glm(LABO~Z_Area+Z_Edge, data=data, family=poisson())
LBstep3.2a<-glm(LABO~Z_Area+Z_Edge+Z_BA_Total, data=data, family=poisson())
LBstep3.2b<-glm(LABO~Z_Area+Z_Edge+Z_DBH_Total, data=data, family=poisson())
LBstep3.2c<-glm(LABO~Z_Area+Z_Edge+Z_Canopy, data=data, family=poisson())
LBstep3.2d<-glm(LABO~Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
LBstep3.2e<-glm(LABO~Z_Area+Z_Edge+Z_BA_Total+Z_Canopy, data=data, family=poisson())
LBstep3.2f<-glm(LABO~Z_Area+Z_Edge+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
LBstep3.2g<-glm(LABO~Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("LBstep3.2a", "LBstep3.2b", "LBstep3.2c", "LBstep3.2d", "LBstep3.2e", "LBstep3.2f", "LBstep3.2g", "LBstep3.2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Labo- model averaging
top.mod.names <- c("LBstep3.1d", "LBstep3.2d", "LBstep3.2g")
top.mods <- list( ) 
for(i in 1:length(top.mod.names)) {top.mods[[i]] <- get(top.mod.names[i]) }

modavg(top.mods, parm="Z_Area", modnames = NULL)
modavg(top.mods, parm="Z_Edge", modnames = NULL)
modavg(top.mods, parm="Z_BA_Total", modnames = NULL)
modavg(top.mods, parm="Z_DBH_Total", modnames = NULL)
modavg(top.mods, parm="Z_Canopy", modnames = NULL)

