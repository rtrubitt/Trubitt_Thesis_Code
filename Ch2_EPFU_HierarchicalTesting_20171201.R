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

#Epfu:
#Step 1: Landscape Level
EFstep1a<-glm(EPFU~Z_TreeCov250, data=data, family=poisson())
EFstep1b<-glm(EPFU~Z_TreeCov500, data=data, family=poisson())
EFstep1c<-glm(EPFU~Z_TreeCov1000, data=data, family=poisson())
EFstep1d<-glm(EPFU~Z_TreeCov3000, data=data, family=poisson())
EFstep1null<-glm(EPFU~1, data=data, family=poisson())

cand.mod.names <- c("EFstep1a", "EFstep1b", "EFstep1c", "EFstep1d", "EFstep1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 2: Patch Level and Landscape/patch interactions
EFstep2null<-glm(EPFU~Z_TreeCov250, data=data, family=poisson())
EFstep2a<-glm(EPFU~Z_TreeCov250+Z_Area, data=data, family=poisson())
EFstep2b<-glm(EPFU~Z_TreeCov250+Z_Edge, data=data, family=poisson())
EFstep2c<-glm(EPFU~Z_TreeCov250+Z_Area+Z_Edge, data=data, family=poisson())
EFstep2d<-glm(EPFU~Z_TreeCov250*Z_Area, data=data, family=poisson())
EFstep2e<-glm(EPFU~Z_TreeCov250*Z_Area+Z_Edge, data=data, family=poisson())
EFstep2f<-glm(EPFU~Z_TreeCov250*Z_Edge, data=data, family=poisson())
EFstep2g<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area, data=data, family=poisson())
EFstep2h<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_TreeCov250*Z_Area, data=data, family=poisson())

cand.mod.names <- c("EFstep2a", "EFstep2b", "EFstep2c", "EFstep2d", "EFstep2e", "EFstep2f", "EFstep2g", "EFstep2h", "EFstep2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.1- Local level with EFstep2g as base model
EFstep3.1null<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area, data=data, family=poisson())
EFstep3.1a<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_BA_Total, data=data, family=poisson())
EFstep3.1b<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_DBH_Total, data=data, family=poisson())
EFstep3.1c<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_Canopy, data=data, family=poisson())
EFstep3.1d<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
EFstep3.1e<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_BA_Total+Z_Canopy, data=data, family=poisson())
EFstep3.1f<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
EFstep3.1g<-glm(EPFU~Z_TreeCov250*Z_Edge+Z_Area+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("EFstep3.1a", "EFstep3.1b", "EFstep3.1c", "EFstep3.1d", "EFstep3.1e", "EFstep3.1f", "EFstep3.1g","EFstep3.1null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#Step 3.2- Local level with EFstep2g as base model
EFstep3.2null<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area, data=data, family=poisson())
EFstep3.2a<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_BA_Total, data=data, family=poisson())
EFstep3.2b<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_DBH_Total, data=data, family=poisson())
EFstep3.2c<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_Canopy, data=data, family=poisson())
EFstep3.2d<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_BA_Total+Z_DBH_Total, data=data, family=poisson())
EFstep3.2e<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_BA_Total+Z_Canopy, data=data, family=poisson())
EFstep3.2f<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_DBH_Total+Z_Canopy, data=data, family=poisson())
EFstep3.2g<-glm(EPFU~Z_TreeCov250+Z_Edge+Z_Area+Z_BA_Total+Z_DBH_Total+Z_Canopy, data=data, family=poisson())

cand.mod.names <- c("EFstep3.2a", "EFstep3.2b", "EFstep3.2c", "EFstep3.2d", "EFstep3.2e", "EFstep3.2f", "EFstep3.2g","EFstep3.2null")
cand.mods <- list( ) 
for(i in 1:length(cand.mod.names)) {cand.mods[[i]] <- get(cand.mod.names[i]) }
weights<-aictab(cand.set = cand.mods, modnames = cand.mod.names)

#EPFU model averaging
top.mod.names <- c("EFstep3.1a", "EFstep3.1d", "EFstep3.1e", "EFstep3.2a", "EFstep3.2d", "EFstep3.2e")
top.mods <- list( ) 
for(i in 1:length(top.mod.names)) {top.mods[[i]] <- get(top.mod.names[i]) }

modavg(top.mods, parm="Z_TreeCov250", modnames = NULL, exclude="Z_TreeCov250:Z_Edge")
modavg(top.mods, parm="Z_Area", modnames = NULL)
modavg(top.mods, parm="Z_Edge", modnames = NULL, exclude="Z_TreeCov250:Z_Edge")
modavg(top.mods, parm="Z_BA_Total", modnames = NULL)
modavg(top.mods, parm="Z_DBH_Total", modnames = NULL)
modavg(top.mods, parm="Z_Canopy", modnames = NULL)
modavg(top.mods, parm="Z_TreeCov250:Z_Edge", modnames = NULL)

#Take a look at the interaction term (I just plotted it with the top model that included the interaction- all three models
#produce really similar graphs)
EFint<-interplot(m=EFstep3.1a, var1="Z_Edge", var2="Z_TreeCov250") +
  xlab("Tree Cover Zscore") + ylab("Estimated Coef for Edge/Area Ratio") +
  geom_hline(yintercept=0, colour="red")

