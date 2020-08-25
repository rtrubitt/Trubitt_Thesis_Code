#Trubitt Masters Thesis Documentation
#Chapter 1: Habitat Associations of Rangeland Bats

#Fieldwork took place from mid May to mid August, 2016. 
#Response variable is minutes with detection (see Miller, 2001).

#Data file. This file contains Z-scores rather than raw data. Z scores are used so that coefficients from variables with different magnitudes
#can be compared. Field descriptions: "Point"- Survey point number. "Detector" - Bat detector ID. "TimeRec"- minutes the detector was recording.
# "EpfuMin etc"- number of minutes with a detection of that species. Epfu= E. fuscus, Labo=L. borealis, Lano= L. noctivagans, Laci=L.cinereus.
#"TotalMin"- number of minutes with any bat detection. "Long" and "Lat"- longitude and latitude. All of the following variables are Z scores of 
#measured variables. The code to measure these variables is included below. Variables that started with P are proximate variables- "distance to" variables.
#Variables starting with L are landscape variables. The Z in the second character indicates that these are Z scores. The third area is the landscape 
#feature being measured- trees, water, structures etc. For landscape level variables, the radius that the variable is measured at is listed last in the
#name- these can be 250, 500, 1000 or 3000m radii. Note- "Ratio" is the ratio of tree edge length to tree area, and "Edge" is just the length of the
#tree edge. I used "Ratio" for modelling, as tree edge length is highly correlated with tree area.
act<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/2016_Chapter1ModelVariables.csv")
act<-act[1:237,]

#Packages:
install.packages("lme4")
install.packages("AICcmodavg")
install.packages("ggplot2")
install.packages("gridExtra")
library(lme4)
library(AICcmodavg)
library(ggplot2)
library(gridExtra)

install.packages("pacman")
pacman::p_load(maps, ggplot2, ggmap, broom, rgdal, rgeos,plyr, maptools )
install.packages("raster")
library(raster)

#Step 1: Data exploration
##Density histograms to check distrubution of bat activity data- determines which family/link function is used for generalized linear models.
##The shown code is for total number of calls, but you can also break it down by species.
Ch1Distribution<- ggplot(act, aes(x=TotalMin)) + theme_bw() + 
  geom_histogram(aes(y=..density..),      
                 binwidth=1,
                 colour="black", fill="lightgreen") +
  geom_density(alpha=.2, fill="#FF6666")

##Use boxplot/Kruskal-Wallis test to test for differences between bat detectors. Determines whether GLM or GLMM is used.
Ch1Det<-ggplot(act, aes(x=Detector, y=TotalMin)) + geom_boxplot()

kruskal.test(act$TotalMin, act$Detector)

#Step 2: GIS code to gather landscape variable data
##Make layers
proj<-readOGR(dsn= "C:/Users/rebecca.trubitt/Desktop/Thesis GIS", layer= "Landcover_Projected")
proj@data$id <- rownames(proj@data)
proj.points <- fortify(proj, region="id")
proj.df <- join(proj.points, proj@data, by="id")

pt<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/Pts_R.csv")
coordinates(pt) = ~x + y

##Make buffers around survey points at 250, 500, 1000, 3000m radii
buff250<-gBuffer(pt186, byid=TRUE, id=NULL, width=250, quadsegs=8, capStyle="ROUND",
                 joinStyle="ROUND", mitreLimit=1.0)

buff500<-gBuffer(pt186, byid=TRUE, id=NULL, width=500, quadsegs=16, capStyle="ROUND",
                 joinStyle="ROUND", mitreLimit=1.0)

buff1000<-gBuffer(pt186, byid=TRUE, id=NULL, width=1000, quadsegs=32, capStyle="ROUND",
                  joinStyle="ROUND", mitreLimit=1.0)

buff3000<-gBuffer(pt186, byid=TRUE, id=NULL, width=3000, quadsegs=96, capStyle="ROUND",
                  joinStyle="ROUND", mitreLimit=1.0)

proj<- gBuffer(proj, byid=TRUE, width=0) #This command helps fix overlaps in hand-drawn layer

##Clip Landcover shapefile to buffers
int250<-intersect(buff250, proj)
int500<-intersect(buff500, proj)
int1000<-intersect(buff1000, proj)
int3000<-intersect(buff3000, proj)

##Calculate landcover areas and aggregate- gathers numbers for Tree Cover, Water Cover, Crop Cover, and Wetland Cover
area250 <- data.frame(area=sapply(int250@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area250) <- sapply(int250@polygons, FUN=function(x) {slot(x, 'ID')})
attArea <- spCbind(int250, area250)
agg250<-aggregate(area~attArea@data$Class+attArea@data$Point.No, data=attArea, FUN=sum)

area500 <- data.frame(area=sapply(int500@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area500) <- sapply(int500@polygons, FUN=function(x) {slot(x, 'ID')})
attArea500 <- spCbind(int500, area500)
agg500<-aggregate(area~attArea500@data$Class+attArea500@data$Point.No, data=attArea500, FUN=sum)

area1000 <- data.frame(area=sapply(int1000@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area1000) <- sapply(int1000@polygons, FUN=function(x) {slot(x, 'ID')})
attArea1000 <- spCbind(int1000, area1000)
agg1000<-aggregate(area~attArea1000@data$Class+attArea1000@data$Point.No, data=attArea1000, FUN=sum)

area3000 <- data.frame(area=sapply(int3000@polygons, FUN=function(x) {slot(x, 'area')}))
row.names(area3000) <- sapply(int3000@polygons, FUN=function(x) {slot(x, 'ID')})
attArea3000 <- spCbind(int3000, area3000)
agg3000<-aggregate(area~attArea3000@data$Class+attArea3000@data$Point.No, data=attArea3000, FUN=sum)

write.csv(agg250, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/area250_date.csv") #Change parameters to save new file for each buffer

##Find landscape level tree edge: Import new shapefile (only tree cover), re-do clips with only trees, calculate and aggregate shape perimeter
tree<-readOGR(dsn= "C:/Users/rebecca.trubitt/Desktop/Thesis GIS", layer= "TreeCover")
tree<- gBuffer(tree, byid=TRUE, width=0)

int250<-intersect(buff250, tree)
int500<-intersect(buff500, tree)
int1000<-intersect(buff1000, tree)
int3000<-intersect(buff3000, tree)

per250 <- data.frame(gLength(int250, byid=TRUE))
rownames(per250)<-sapply(int250@polygons, FUN=function(x) {slot(x, 'ID')})
attPer250 <- spCbind(int250, per250)
aggper250<-aggregate(attPer250@data$gLength.int250..byid...TRUE.~attPer250@data$Class+attPer250@data$Point.No, data=attPer250, FUN=sum)
write.csv(aggper250, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/edge250_date.csv")


per500 <- data.frame(gLength(int500, byid=TRUE))
rownames(per500)<-sapply(int500@polygons, FUN=function(x) {slot(x, 'ID')})
attPer500 <- spCbind(int500, per500)
aggper500<-aggregate(attPer500@data$gLength.int500..byid...TRUE.~attPer500@data$Class+attPer500@data$Point.No, data=attPer500, FUN=sum)
write.csv(aggper500, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/edge500_date.csv")


per1000 <- data.frame(gLength(int1000, byid=TRUE))
rownames(per1000)<-sapply(int1000@polygons, FUN=function(x) {slot(x, 'ID')})
attPer1000 <- spCbind(int1000, per1000)
aggper1000<-aggregate(attPer1000@data$gLength.int1000..byid...TRUE.~attPer1000@data$Class+attPer1000@data$Point.No, data=attPer1000, FUN=sum)
write.csv(aggper1000, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/edge1000_date.csv")


per3000 <- data.frame(gLength(int3000, byid=TRUE))
rownames(per3000)<-sapply(int3000@polygons, FUN=function(x) {slot(x, 'ID')})
attPer3000 <- spCbind(int3000, per3000)
aggper3000<-aggregate(attPer3000@data$gLength.int3000..byid...TRUE.~attPer3000@data$Class+attPer3000@data$Point.No, data=attPer3000, FUN=sum)
write.csv(aggper3000, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/edge3000_date.csv")

##Find road density with similar protocol to that for tree edge length
road<-readOGR(dsn= "C:/Users/rebecca.trubitt/Desktop/Thesis GIS/Roads", layer= "Road_combo")

int250<-intersect(road, buff250)
int500<-intersect(road, buff500)
int1000<-intersect(road, buff1000)
int3000<-intersect(road, buff3000)

rd250 <- data.frame(gLength(int250, byid=TRUE))
rownames(rd250)<-sapply(int250@lines, FUN=function(x) {slot(x, 'ID')})
attPer250 <- spCbind(int250, rd250)
aggper250<-aggregate(attPer250@data$gLength.int250..byid...TRUE.~attPer250@data$Point.No, data=attPer250, FUN=sum)
write.csv(aggper250, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/rd250_date.csv")

rd500 <- data.frame(gLength(int500, byid=TRUE))
rownames(rd500)<-sapply(int500@lines, FUN=function(x) {slot(x, 'ID')})
attPer500 <- spCbind(int500, rd500)
aggper500<-aggregate(attPer500@data$gLength.int500..byid...TRUE.~attPer500@data$Point.No, data=attPer500, FUN=sum)
write.csv(aggper500, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/rd500_date.csv")

rd1000 <- data.frame(gLength(int1000, byid=TRUE))
rownames(rd1000)<-sapply(int1000@lines, FUN=function(x) {slot(x, 'ID')})
attPer1000 <- spCbind(int1000, rd1000)
aggper1000<-aggregate(attPer1000@data$gLength.int1000..byid...TRUE.~attPer1000@data$Point.No, data=attPer1000, FUN=sum)
write.csv(aggper1000, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/rd1000_date.csv")

rd3000 <- data.frame(gLength(int3000, byid=TRUE))
rownames(rd3000)<-sapply(int3000@lines, FUN=function(x) {slot(x, 'ID')})
attPer3000 <- spCbind(int3000, rd3000)
aggper3000<-aggregate(attPer3000@data$gLength.int3000..byid...TRUE.~attPer3000@data$Point.No, data=attPer3000, FUN=sum)
write.csv(aggper3000, file="C:/Users/rebecca.trubitt/Desktop/ResearchR/data/rd3000_date.csv")

##Additional variables (Distance to tree, distance to open water, distance to human built structure) were measured manually in ArcGIS.
##There is also a 'distance to' protocol in ArcGIS that you can use if you have a shapefile of the landscape feature in question.

#Step 3: Use vifmer to check for collinearity between predictor variables
##vif.mer- This function calculates Variance Inflation Factors for each variable in a model. I required VIF factors to be lower than 3, following
##the protocol suggested by Zuur et al, 2009. The function for vif.mer (vif function for mixed models) is from https://github.com/aufrank/R-hacks/blob/master/mer-utils.R.

vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}

LN.global.250<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LN.global.500<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LN.global.1000<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)

vif.mer(LN.global.1000)

LC.global.250<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LC.global.500<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LC.global.1000<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)

vif.mer(LC.global.1000)


EF.global.250<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
EF.global.500<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
EF.global.1000<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)

vif.mer(EF.global.1000)

#Step 4: Run GLMMs and compare using AICc

##L. noctivagans
LN.global.250<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LN.landscape.250<-glmer(LanoMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LN.landcover.250<-glmer(LanoMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+(1|Detector), data=act, family=poisson)
LN.proximate.250<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LN.roost.250<-glmer(LanoMin ~ PZTree+PZStruct+LZTree250+(1|Detector), data=act, family=poisson)
LN.tree.250<-glmer(LanoMin ~ PZTree+LZTree250+LZRatio250+(1|Detector), data=act, family=poisson)
LN.water.250<-glmer(LanoMin ~ PZWater+LZWater250+LZWetland250+(1|Detector), data=act, family=poisson)
LN.devo.250<-glmer(LanoMin ~ PZStruct+LZCrop250+LZRoad250+(1|Detector), data=act, family=poisson)
LN.null.250<-glmer(LanoMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LN.global.250, LN.landscape.250, LN.landcover.250, LN.proximate.250, LN.roost.250, LN.tree.250, LN.water.250, LN.devo.250, LN.null.250, k=2, REML=NULL)
##chooses global

LN.global.500<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LN.landscape.500<-glmer(LanoMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LN.landcover.500<-glmer(LanoMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+(1|Detector), data=act, family=poisson)
LN.proximate.500<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LN.roost.500<-glmer(LanoMin ~ PZTree+PZStruct+LZTree500+(1|Detector), data=act, family=poisson)
LN.tree.500<-glmer(LanoMin ~ PZTree+LZTree500+LZRatio500+(1|Detector), data=act, family=poisson)
LN.water.500<-glmer(LanoMin ~ PZWater+LZWater500+LZWetland500+(1|Detector), data=act, family=poisson)
LN.devo.500<-glmer(LanoMin ~ PZStruct+LZCrop500+LZRoad500+(1|Detector), data=act, family=poisson)
LN.null.500<-glmer(LanoMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LN.global.500, LN.landscape.500, LN.landcover.500, LN.proximate.500, LN.roost.500, LN.tree.500, LN.water.500, LN.devo.500, LN.null.500, k=2, REML=NULL)
##chooses global

LN.global.1000<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LN.landscape.1000<-glmer(LanoMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LN.landcover.1000<-glmer(LanoMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+(1|Detector), data=act, family=poisson)
LN.proximate.1000<-glmer(LanoMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LN.roost.1000<-glmer(LanoMin ~ PZTree+PZStruct+LZTree1000+(1|Detector), data=act, family=poisson)
LN.tree.1000<-glmer(LanoMin ~ PZTree+LZTree1000+LZRatio1000+(1|Detector), data=act, family=poisson)
LN.water.1000<-glmer(LanoMin ~ PZWater+LZWater1000+LZWetland1000+(1|Detector), data=act, family=poisson)
LN.devo.1000<-glmer(LanoMin ~ PZStruct+LZCrop1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LN.null.1000<-glmer(LanoMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LN.global.1000, LN.landscape.1000, LN.landcover.1000, LN.proximate.1000, LN.roost.1000, LN.tree.1000, LN.water.1000, LN.devo.1000, LN.null.1000, k=2, REML=NULL)
##Chooses global

AICc(LN.global.1000, LN.global.500, LN.global.250, k=2, REML=NULL)
##Chooses 500

##L. cinereus
LC.global.250<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LC.landscape.250<-glmer(LaciMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LC.landcover.250<-glmer(LaciMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+(1|Detector), data=act, family=poisson)
LC.proximate.250<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LC.roost.250<-glmer(LaciMin ~ PZTree+LZTree250+(1|Detector), data=act, family=poisson, offset=log(TimeRec))
LC.tree.250<-glmer(LaciMin ~ PZTree+LZTree250+LZRatio250+(1|Detector), data=act, family=poisson)
LC.water.250<-glmer(LaciMin ~ PZWater+LZWater250+LZWetland250+(1|Detector), data=act, family=poisson)
LC.devo.250<-glmer(LaciMin ~ PZStruct+LZCrop250+LZRoad250+(1|Detector), data=act, family=poisson)
LC.null.250<-glmer(LaciMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LC.global.250, LC.landscape.250, LC.landcover.250, LC.proximate.250, LC.roost.250, LC.tree.250, LC.water.250, LC.devo.250, LC.null.250, k=2, REML=NULL)
##Chooses global

LC.global.500<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LC.landscape.500<-glmer(LaciMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LC.landcover.500<-glmer(LaciMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+(1|Detector), data=act, family=poisson)
LC.proximate.500<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LC.roost.500<-glmer(LaciMin ~ PZTree+LZTree500+(1|Detector), data=act, family=poisson)
LC.tree.500<-glmer(LaciMin ~ PZTree+LZTree500+LZRatio500+(1|Detector), data=act, family=poisson)
LC.water.500<-glmer(LaciMin ~ PZWater+LZWater500+LZWetland500+(1|Detector), data=act, family=poisson)
LC.devo.500<-glmer(LaciMin ~ PZStruct+LZCrop500+LZRoad500+(1|Detector), data=act, family=poisson)
LC.null.500<-glmer(LaciMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LC.global.500, LC.landscape.500, LC.landcover.500, LC.proximate.500, LC.roost.500, LC.tree.500, LC.water.500, LC.devo.500, LC.null.500, k=2, REML=NULL)
##Chooses global

LC.global.1000<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LC.landscape.1000<-glmer(LaciMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LC.landcover.1000<-glmer(LaciMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+(1|Detector), data=act, family=poisson)
LC.proximate.1000<-glmer(LaciMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LC.roost.1000<-glmer(LaciMin ~ PZTree+LZTree1000+(1|Detector), data=act, family=poisson)
LC.tree.1000<-glmer(LaciMin ~ PZTree+LZTree1000+LZRatio1000+(1|Detector), data=act, family=poisson)
LC.water.1000<-glmer(LaciMin ~ PZWater+LZWater1000+LZWetland1000+(1|Detector), data=act, family=poisson)
LC.devo.1000<-glmer(LaciMin ~ PZStruct+LZCrop1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LC.null.1000<-glmer(LaciMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LC.global.1000, LC.landscape.1000, LC.landcover.1000, LC.proximate.1000, LC.roost.1000, LC.tree.1000, LC.water.1000, LC.devo.1000, LC.null.1000, k=2, REML=NULL)
##Chooses global

AICc(LC.global.1000, LC.global.500, LC.global.250, k=2, REML=NULL)
##Chooses 500

##E. fuscus
EF.global.250<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
EF.landscape.250<-glmer(EpfuMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
EF.landcover.250<-glmer(EpfuMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+(1|Detector), data=act, family=poisson)
EF.proximate.250<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
EF.roost.250<-glmer(EpfuMin ~ PZTree+PZStruct+LZTree250+(1|Detector), data=act, family=poisson)
EF.tree.250<-glmer(EpfuMin ~ PZTree+LZTree250+LZRatio250+(1|Detector), data=act, family=poisson)
EF.water.250<-glmer(EpfuMin ~ PZWater+LZWater250+LZWetland250+(1|Detector), data=act, family=poisson)
EF.devo.250<-glmer(EpfuMin ~ PZStruct+LZCrop250+LZRoad250+(1|Detector), data=act, family=poisson)
EF.null.250<-glmer(EpfuMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(EF.global.250, EF.landscape.250, EF.landcover.250, EF.proximate.250, EF.roost.250, EF.tree.250, EF.water.250, EF.devo.250, EF.null.250, k=2, REML=NULL)
#Chooses global

EF.global.500<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
EF.landscape.500<-glmer(EpfuMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
EF.landcover.500<-glmer(EpfuMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+(1|Detector), data=act, family=poisson)
EF.proximate.500<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
EF.roost.500<-glmer(EpfuMin ~ PZTree+PZStruct+LZTree500+(1|Detector), data=act, family=poisson)
EF.tree.500<-glmer(EpfuMin ~ PZTree+LZTree500+LZRatio500+(1|Detector), data=act, family=poisson)
EF.water.500<-glmer(EpfuMin ~ PZWater+LZWater500+LZWetland500+(1|Detector), data=act, family=poisson)
EF.devo.500<-glmer(EpfuMin ~ PZStruct+LZCrop500+LZRoad500+(1|Detector), data=act, family=poisson)
EF.null.500<-glmer(EpfuMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(EF.global.500, EF.landscape.500, EF.landcover.500, EF.proximate.500, EF.roost.500, EF.tree.500, EF.water.500, EF.devo.500, EF.null.500, k=2, REML=NULL)
#Chooses global

EF.global.1000<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
EF.landscape.1000<-glmer(EpfuMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
EF.landcover.1000<-glmer(EpfuMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+(1|Detector), data=act, family=poisson)
EF.proximate.1000<-glmer(EpfuMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
EF.roost.1000<-glmer(EpfuMin ~ PZTree+PZStruct+LZTree1000+(1|Detector), data=act, family=poisson)
EF.tree.1000<-glmer(EpfuMin ~ PZTree+LZTree1000+LZRatio1000+(1|Detector), data=act, family=poisson)
EF.water.1000<-glmer(EpfuMin ~ PZWater+LZWater1000+LZWetland1000+(1|Detector), data=act, family=poisson)
EF.devo.1000<-glmer(EpfuMin ~ PZStruct+LZCrop1000+LZRoad1000+(1|Detector), data=act, family=poisson)
EF.null.1000<-glmer(EpfuMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(EF.global.1000, EF.landscape.1000, EF.landcover.1000, EF.proximate.1000, EF.roost.1000, EF.tree.1000, EF.water.1000, EF.devo.1000, EF.null.1000, k=2, REML=NULL)
#Chooses global

AICc(EF.global.1000, EF.global.500, EF.global.250, k=2, REML=NULL)
#Chooses 1000

##L. borealis
##Does not converge- LB.global.250<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LB.landscape.250<-glmer(LaboMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+LZRatio250+LZRoad250+(1|Detector), data=act, family=poisson)
LB.landcover.250<-glmer(LaboMin ~ LZTree250+LZWater250+LZWetland250+LZCrop250+(1|Detector), data=act, family=poisson)
LB.proximate.250<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LB.roost.250<-glmer(LaboMin ~ PZTree+LZTree250+(1|Detector), data=act, family=poisson)
LB.tree.250<-glmer(LaboMin ~ PZTree+LZTree250+LZRatio250+(1|Detector), data=act, family=poisson)
LB.water.250<-glmer(LaboMin~ PZWater+LZWater250+LZWetland250+(1|Detector), data=act, family=poisson)
LB.devo.250<-glmer(LaboMin ~ PZStruct+LZCrop250+LZRoad250+(1|Detector), data=act, family=poisson)
LB.null.250<-glmer(LaboMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LB.landscape.250, LB.landcover.250, LB.proximate.250, LB.roost.250, LB.tree.250, LB.water.250, LB.devo.250, LB.null.250, k=2, REML=NULL)
##Chooses landcover

##Does not converge LB.global.500<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson, offset=log(TimeRec))
LB.landscape.500<-glmer(LaboMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+LZRatio500+LZRoad500+(1|Detector), data=act, family=poisson)
LB.landcover.500<-glmer(LaboMin ~ LZTree500+LZWater500+LZWetland500+LZCrop500+(1|Detector), data=act, family=poisson)
LB.proximate.500<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LB.roost.500<-glmer(LaboMin ~ PZTree+LZTree500+(1|Detector), data=act, family=poisson)
LB.tree.500<-glmer(LaboMin ~ PZTree+LZTree500+LZRatio500+(1|Detector), data=act, family=poisson)
LB.water.500<-glmer(LaboMin ~ PZWater+LZWater500+LZWetland500+(1|Detector), data=act, family=poisson)
LB.devo.500<-glmer(LaboMin ~ PZStruct+LZCrop500+LZRoad500+(1|Detector), data=act, family=poisson)
LB.null.500<-glmer(LaboMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LB.landscape.500, LB.landcover.500, LB.proximate.500, LB.roost.500, LB.tree.500, LB.water.500, LB.devo.500, LB.null.500, k=2, REML=NULL)
##Chooses landcover

##Does not converge LB.global.1000<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson, offset=log(TimeRec))
LB.landscape.1000<-glmer(LaboMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+LZRatio1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LB.landcover.1000<-glmer(LaboMin ~ LZTree1000+LZWater1000+LZWetland1000+LZCrop1000+(1|Detector), data=act, family=poisson)
LB.proximate.1000<-glmer(LaboMin ~ PZTree+PZWater+PZStruct+(1|Detector), data=act, family=poisson)
LB.roost.1000<-glmer(LaboMin ~ PZTree+LZTree1000+(1|Detector), data=act, family=poisson)
LB.tree.1000<-glmer(LaboMin ~ PZTree+LZTree1000+LZRatio1000+(1|Detector), data=act, family=poisson)
LB.water.1000<-glmer(LaboMin ~ PZWater+LZWater1000+LZWetland1000+(1|Detector), data=act, family=poisson)
LB.devo.1000<-glmer(LaboMin ~ PZStruct+LZCrop1000+LZRoad1000+(1|Detector), data=act, family=poisson)
LB.null.1000<-glmer(LaboMin ~ 1+(1|Detector), data=act, family=poisson)
AICc(LB.landscape.1000, LB.landcover.1000, LB.proximate.1000, LB.roost.1000, LB.tree.1000, LB.water.1000, LB.devo.1000, LB.null.1000, k=2, REML=NULL)
##Chooses landcover

AICc(LB.landcover.1000, LB.landcover.500, LB.landcover.250, k=2, REML=NULL)
##250 and 500 equally explainatory, 250 AIC slightly lower


#Step 5: Determine coefficients and confidence intervals for top models
LN.global.500
confint(LN.global.500)

LC.global.500
confint(LC.global.500)

EF.global.1000
confint(EF.global.1000)

LB.landcover.250
confint(LB.landcover.250)

LB.landcover.500
confint(LB.landcover.500)


#Step 6: Graph
lanocoef<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/LANO_FigureCoefficients_20170914.csv")
lanocoef$Var<- factor(lanocoef$Var, as.character(lanocoef$Var))

LANO <- ggplot(lanocoef, aes(x=Var, y=LANO_500)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL500,ymax=CIU500), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. noctivagans", y="Coefficients")+
  theme(axis.text.y = element_text(size="16", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  expand_limits(y=c(-0.8, 0.8))+
  coord_flip()
LANO


lacicoef<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/LACI_FigureCoefficients_20170914.csv")
lacicoef$Var<- factor(lacicoef$Var, as.character(lacicoef$Var))

LACI <- ggplot(lacicoef, aes(x=Var, y=LACI_500)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL500,ymax=CIU500), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. cinereus", y="Coefficients")+
  theme(axis.text.y = element_text(size="16", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  expand_limits(y=c(-0.8, 0.8))+
  coord_flip()
LACI


epfucoef<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/EPFU_FigureCoefficients_20170914.csv")
epfucoef$Var<- factor(epfucoef$Var, as.character(epfucoef$Var))

EPFU <- ggplot(epfucoef, aes(x=Var, y=EPFU_1000)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL1000,ymax=CIU1000), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="E. fuscus", y="Coefficients")+
  theme(axis.text.y = element_text(size="16", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  expand_limits(y=c(-0.8, 0.8))+
  coord_flip()
EPFU


labocoef<-read.csv("C:/Users/rebecca.trubitt/Desktop/ResearchR/data/LABO_FigureCoefficients_20170914.csv")
labocoef$Var<- factor(labocoef$Var, as.character(labocoef$Var))
#This graph is only for the 250m radius- this code can also be used with 500m coef/confints from step 5.
LABO <- ggplot(labocoef, aes(x=Var, y=LABO250)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL250,ymax=CIU250), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. borealis", y="Coefficients")+
  theme(axis.text.y = element_text(size="14", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  expand_limits(y=c(-0.8, 0.8))+
  coord_flip()
LABO

x11()
grid.arrange(LANO, LACI, EPFU, LABO, ncol=2)

