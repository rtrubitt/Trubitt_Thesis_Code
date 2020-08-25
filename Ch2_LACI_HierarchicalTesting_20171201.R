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

#Laci:
#Step 1: Landscape Level
LCstep1a<-glm(LACI~Z_TreeCov250, data=data, family=poisson())
LCstep1b<-glm(LACI~Z_TreeCov500, data=data, family=poisson())
LCstep1c<-glm(LACI~Z_TreeCov1000, data=data, family=poisson())
LCstep1d<-glm(LACI~Z_TreeCov3000, data=data, family=poisson())
LCstep1null<-glm(LACI~1, data=data, family=poisson())

cand.mod.names <- c("LCstep1a", "LCstep1b", "LCstep1c", "LCstep1d", "LCstep1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 2.1: Patch Level and Landscape/patch interactions- Base model is LCstep1b
LCstep2.1null<-glm(LACI~Z_TreeCov500, data=data, family=poisson())
LCstep2.1a<-glm(LACI~Z_TreeCov500+Z_Area, data=data, family=poisson())
LCstep2.1b<-glm(LACI~Z_TreeCov500+Z_Edge, data=data, family=poisson())
LCstep2.1c<-glm(LACI~Z_TreeCov500+Z_Area+Z_Edge, data=data, family=poisson())
LCstep2.1d<-glm(LACI~Z_TreeCov500*Z_Area, data=data, family=poisson())
LCstep2.1e<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge, data=data, family=poisson())
LCstep2.1f<-glm(LACI~Z_TreeCov500*Z_Edge, data=data, family=poisson())
LCstep2.1g<-glm(LACI~Z_TreeCov500*Z_Edge+Z_Area, data=data, family=poisson())
LCstep2.1h<-glm(LACI~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area, data=data, family=poisson())

cand.mod.names <- c("LCstep2.1a", "LCstep2.1b", "LCstep2.1c", "LCstep2.1d", "LCstep2.1e", "LCstep2.1f", "LCstep2.1g", "LCstep2.1h", "LCstep2.1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 2.2: Patch level- base model is LCstep1null- no interactions as cover is not carried forward
LCstep2.2null<-glm(LACI~1, data=data, family=poisson())
LCstep2.2a<-glm(LACI~Z_Area, data=data, family=poisson())
LCstep2.2b<-glm(LACI~Z_Edge, data=data, family=poisson())
LCstep2.2c<-glm(LACI~Z_Area+Z_Edge, data=data, family=poisson())

cand.mod.names <- c("LCstep2.2a", "LCstep2.2b", "LCstep2.2c", "LCstep2.2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.1- Local level- base model is LCstep2.1e
LCstep3.1null<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge, data=data, family=poisson())
LCstep3.1a<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_BA_Total, data=data, family=poisson())
LCstep3.1b<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_DBH_Total, data=data, family=poisson())
LCstep3.1c<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_Canopy, data=data, family=poisson())
LCstep3.1d<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
LCstep3.1e<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_BA_Total+Z_Canopy, data=data, family=poisson())
LCstep3.1f<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
LCstep3.1g<-glm(LACI~Z_TreeCov500*Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("LCstep3.1a", "LCstep3.1b", "LCstep3.1c", "LCstep3.1d", "LCstep3.1e", "LCstep3.1f", "LCstep3.1g", "LCstep3.1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.2- Local level- base model is LCstep2.2c
LCstep3.2null<-glm(LACI~Z_Area+Z_Edge, data=data, family=poisson())
LCstep3.2a<-glm(LACI~Z_Area+Z_Edge+Z_BA_Total, data=data, family=poisson())
LCstep3.2b<-glm(LACI~Z_Area+Z_Edge+Z_DBH_Total, data=data, family=poisson())
LCstep3.2c<-glm(LACI~Z_Area+Z_Edge+Z_Canopy, data=data, family=poisson())
LCstep3.2d<-glm(LACI~Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
LCstep3.2e<-glm(LACI~Z_Area+Z_Edge+Z_BA_Total+Z_Canopy, data=data, family=poisson())
LCstep3.2f<-glm(LACI~Z_Area+Z_Edge+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
LCstep3.2g<-glm(LACI~Z_Area+Z_Edge+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("LCstep3.2a", "LCstep3.2b", "LCstep3.2c", "LCstep3.2d", "LCstep3.2e", "LCstep3.2f", "LCstep3.2g", "LCstep3.2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Laci- model averaging
top.mod.names <- c("LCstep3.1f", "LCstep3.1b", "LCstep3.2f")
top.mods <- list( ) 
for(i in 1:length(top.mod.names)) {top.mods[[i]] <- get(top.mod.names[i]) }

modavg(top.mods, parm="Z_Area", modnames = NULL, exclude="Z_TreeCov500:Z_Area")
modavg(top.mods, parm="Z_Edge", modnames = NULL)
modavg(top.mods, parm="Z_DBH_Total", modnames = NULL)
modavg(top.mods, parm="Z_Canopy", modnames = NULL)
modavg(top.mods, parm="Z_TreeCov500:Z_Area", modnames = NULL)

#Graphing interaction
LCint<-interplot(m=LCstep3.1f, var1="Z_Area", var2="Z_TreeCov500") +
  xlab("Tree Cover Zscore") + ylab("Estimated Coef for Area") +
  geom_hline(yintercept=0, colour="red")

