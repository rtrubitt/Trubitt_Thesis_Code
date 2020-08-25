install.packages("ggplot2")
install.packages("gridExtra")
library(ggplot2)
library(gridExtra)


lanocoef<-read.csv("F:/Everything Else/ResearchR/data/CH2_LANOCoefs.csv")
lanocoef<-lanocoef[1:5,]

lanocoef$Var<- factor(lanocoef$Var, as.character(lanocoef$Var))

LANO <- ggplot(lanocoef, aes(x=Var, y=LANO)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL,ymax=CIU), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. noctivagans", y="Coefficients")+
  theme(axis.text.y = element_text(size="14", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  coord_flip()
LANO


lacicoef<-read.csv("F:/Everything Else/ResearchR/data/CH2_LACICoefs.csv")
lacicoef<-lacicoef[1:5,]

lacicoef$Var<- factor(lacicoef$Var, as.character(lacicoef$Var))

LACI <- ggplot(lacicoef, aes(x=Var, y=LACI)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL,ymax=CIU), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. cinereus", y="Coefficients")+
  theme(axis.text.y = element_text(size="14", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  coord_flip()
LACI


epfucoef<-read.csv("F:/Everything Else/ResearchR/data/CH2_EPFUCoefs.csv")
epfucoef<-epfucoef[1:7,]

epfucoef$Var<- factor(epfucoef$Var, as.character(epfucoef$Var))

EPFU <- ggplot(epfucoef, aes(x=Var, y=EPFU)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL,ymax=CIU), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="E. fuscus", y="Coefficients")+
  theme(axis.text.y = element_text(size="14", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  coord_flip()
EPFU


labocoef<-read.csv("F:/Everything Else/ResearchR/data/CH2_LABOCoefs.csv")
labocoef<-labocoef[1:5,]

labocoef$Var<- factor(labocoef$Var, as.character(labocoef$Var))

LABO <- ggplot(labocoef, aes(x=Var, y=LABO)) + 
  geom_point(stat="identity") + theme_bw()+
  geom_hline(yintercept=0, colour="red")+
  geom_errorbar(aes(ymin=CIL,ymax=CIU), width=0.3)+
  scale_x_discrete(name="Model Term")+
  labs(title="L. borealis", y="Coefficients")+
  theme(axis.text.y = element_text(size="14", vjust=0.6), axis.text.x= element_text(size="14"), plot.title=element_text(size="16"))+
  coord_flip()
LABO