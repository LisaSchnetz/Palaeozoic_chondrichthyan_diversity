
library(nlme)
library(ggplot2)
library(tidyverse)
library(deeptime)
library(MuMIn)

########################################################################
####################GLS analyses between diversity curves#################
########################################################################

raw <-read.csv("./data/raw_data.csv")

sqs <-read.csv("./data/sqs_data.csv")

squares <-read.csv("./data/squares_data.csv")

q0.5 <- subset(sqs,quorum_level=="0.5")
q0.6 <- subset(sqs,quorum_level=="0.6")
q0.7 <- subset(sqs,quorum_level=="0.7")
q0.8 <- subset(sqs,quorum_level=="0.8")


q0.5richness <- q0.5$qD

q0.6richness <- q0.6$qD
q0.7richness <- q0.7$qD
q0.8richness <- q0.8$qD


rawrichness <- raw$count_taxa
rawrichness1 <- rawrichness[-c(1,2,3,4)]
                            
squaresrichness <- squares$squares_list
squaresrichness1 <- squaresrichness[-c(1,2,3,4)]

gls1<-gls(rawrichness ~ squaresrichness,method = "ML", correlation = corARMA(p=1),na.action = na.omit) 
summary(gls1)
r.squaredLR(gls1) 

cor.test(rawrichness, squaresrichness, method ="pearson")

gls2<-gls(rawrichness1 ~ q0.5richness,method = "ML", correlation = corARMA(p=1),na.action = na.omit) 
summary(gls2)
r.squaredLR(gls2) 

cor.test(rawrichness1, q0.5richness, method ="pearson")


gls3<-gls(rawrichness1 ~ q0.6richness,method = "ML", correlation = corARMA(p=1)) 
summary(gls3)
r.squaredLR(gls3) 

cor.test(rawrichness1, q0.6richness, method ="pearson")

gls4<-gls(rawrichness1 ~ q0.7richness,method = "ML", correlation = corARMA(p=1)) 
summary(gls4)
r.squaredLR(gls4) 

cor.test(rawrichness1, q0.7richness, method ="pearson")

gls5<-gls(rawrichness1 ~ q0.8richness,method = "ML", correlation = corARMA(p=1),na.action = na.omit) 
summary(gls5)
r.squaredLR(gls5)

cor.test(rawrichness1, q0.8richness, method ="pearson")

gls7<-gls(squaresrichness1 ~ q0.5richness,method = "ML", correlation = corARMA(p=1)) 
summary(gls7)
r.squaredLR(gls7) 

cor.test(squaresrichness1, q0.5richness, method ="pearson")

gls8<-gls(squaresrichness1 ~ q0.6richness,method = "ML", correlation = corARMA(p=1)) 
summary(gls8)
r.squaredLR(gls8) 

cor.test(squaresrichness1, q0.6richness, method ="pearson")

gls9<-gls(squaresrichness1 ~ q0.7richness,method = "ML", correlation = corARMA(p=1)) 
summary(gls9)
r.squaredLR(gls9) 

cor.test(squaresrichness1, q0.7richness, method ="pearson")


gls10<-gls(squaresrichness1 ~ q0.8richness,method = "ML", correlation = corARMA(p=1),na.action = na.omit) 
summary(gls10)
r.squaredLR(gls10) 

cor.test(squaresrichness1, q0.8richness, method ="pearson")
