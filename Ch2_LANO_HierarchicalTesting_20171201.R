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

#Lano:
#Step 1: Landscape Level
LNstep1a<-glm(LANO~Z_TreeCov250, data=data, family=poisson())
LNstep1b<-glm(LANO~Z_TreeCov500, data=data, family=poisson())
LNstep1c<-glm(LANO~Z_TreeCov1000, data=data, family=poisson())
LNstep1d<-glm(LANO~Z_TreeCov3000, data=data, family=poisson())
LNstep1null<-glm(LANO~1, data=data, family=poisson())

cand.mod.names <- c("LNstep1a", "LNstep1b", "LNstep1c", "LNstep1d", "LNstep1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 2: Patch Level and Landscape/patch interactions
LNstep2null<-glm(LANO~Z_TreeCov500, data=data, family=poisson())
LNstep2a<-glm(LANO~Z_TreeCov500+Z_Area, data=data, family=poisson())
LNstep2b<-glm(LANO~Z_TreeCov500+Z_Edge, data=data, family=poisson())
LNstep2c<-glm(LANO~Z_TreeCov500+Z_Area+Z_Edge, data=data, family=poisson())
LNstep2d<-glm(LANO~Z_TreeCov500*Z_Area, data=data, family=poisson())
LNstep2e<-glm(LANO~Z_TreeCov500*Z_Area+Z_Edge, data=data, family=poisson())
LNstep2f<-glm(LANO~Z_TreeCov500*Z_Edge, data=data, family=poisson())
LNstep2g<-glm(LANO~Z_TreeCov500*Z_Edge+Z_Area, data=data, family=poisson())
LNstep2h<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area, data=data, family=poisson())

cand.mod.names <- c("LNstep2a", "LNstep2b", "LNstep2c", "LNstep2d", "LNstep2e", "LNstep2f", "LNstep2g", "LNstep2h", "LNstep2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3- Local level
LNstep3null<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area, data=data, family=poisson())
LNstep3a<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_BA_Total, data=data, family=poisson())
LNstep3b<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_DBH_Total, data=data, family=poisson())
LNstep3c<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_Canopy, data=data, family=poisson())
LNstep3d<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
LNstep3e<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_BA_Total+Z_Canopy, data=data, family=poisson())
LNstep3f<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
LNstep3g<-glm(LANO~Z_TreeCov500*Z_Edge+Z_TreeCov500*Z_Area+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("LNstep3a", "LNstep3b", "LNstep3c", "LNstep3d", "LNstep3e", "LNstep3f", "LNstep3g", "LNstep3null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#LANO model averaging
top.mod.names <- c("LNstep3a", "LNstep3d", "LNstep3e")
top.mods <- list( ) 
for(i in 1:length(top.mod.names)) {top.mods[[i]] <- get(top.mod.names[i]) }

modavg(top.mods, parm="Z_BA_Total", modnames = NULL)
modavg(top.mods, parm="Z_DBH_Total", modnames = NULL)
modavg(top.mods, parm="Z_Canopy", modnames = NULL)
modavg(top.mods, parm="Z_TreeCov500:Z_Area", modnames = NULL)
modavg(top.mods, parm="Z_TreeCov500:Z_Edge", modnames = NULL)

#Interactions
LNintedge<-interplot(m=LNstep3a, var1="Z_Edge", var2="Z_TreeCov500") +
  xlab("Tree Cover Zscore") + ylab("Estimated Coef for Edge/Area Ratio") +
  geom_hline(yintercept=0, colour="red")

LNintarea<-interplot(m=LNstep3a, var1="Z_Area", var2="Z_TreeCov500") +
  xlab("Tree Cover Zscore") + ylab("Estimated Coef for Area") +
  geom_hline(yintercept=0, colour="red")
