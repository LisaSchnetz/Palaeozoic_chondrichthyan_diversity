###########################################################################
#                                                                         #
#                       Spatial diversity analyses                        #  
#                                                                         #
###########################################################################
#                                                                         #
#                        Lisa Schnetz - Oktober 2023                      #
#                                                                         #
###########################################################################

## Packages used in this script:
install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("deeptime")
install.packages("ggpubr")
install.packages("divvy")
install.packages("terra")
install.packages("divDyn")
install.packages("rgplates")
install.packages("chronosphere")

# Libraries used in this script
library(tidyverse)
library(dplyr)
library(ggplot2)
library(deeptime)
library(ggpubr)
library(divvy)
library(terra)
library(divDyn)
library(rgplates)
library(chronosphere)

## First make sure that your environment is clean so that you don't mix up data
rm(list=ls()) 

## Read in the datasets:
shark_data <- read.csv("./data/Total_chondrichthyes_coordinates_new.csv")

intervals <- read.csv2("./data/Intervals.csv")

# Subset data to continents
shark_NA <- subset(shark_data, CONTINENT=="North America")

shark_Europe <- subset(shark_data, CONTINENT=="Europe")

shark_Asia <- subset(shark_data, CONTINENT=="Asia")

shark_Africa <- subset(shark_data, CONTINENT=="Africa")

shark_SA <- subset(shark_data, CONTINENT=="South America")

shark_Ant <- subset(shark_data, CONTINENT=="Antarctica")

shark_Aus <- subset(shark_data, CONTINENT=="Australia and Oceania")

####################################################################
############################## Raw curves ############################
#####################################################################
## First: Calculate raw richness counts

# North America counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_NA %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

NA_data<-as.data.frame(t(genus_data))
NA_data<-cbind(NA_data, intervals)

# Europe counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_Europe %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

Europe_data<-as.data.frame(t(genus_data))
Europe_data<-cbind(Europe_data, intervals)

# Asia counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_Asia %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

Asia_data<-as.data.frame(t(genus_data))
Asia_data<-cbind(Asia_data, intervals)

# Africa counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_Africa %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

Africa_data<-as.data.frame(t(genus_data))
Africa_data<-cbind(Africa_data, intervals)
         
# South America counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_SA %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

SA_data<-as.data.frame(t(genus_data))
SA_data<-cbind(SA_data, intervals)

# Antarctica counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_Ant %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

Ant_data<-as.data.frame(t(genus_data))
Ant_data<-cbind(Ant_data, intervals)

# Australia counts
genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- shark_Aus %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

Aus_data<-as.data.frame(t(genus_data))
Aus_data<-cbind(Aus_data, intervals)

## Second: Plot the data

rawplot_continents <- ggplot()+
  
  #set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+
  
  #manually add grey bars to differentiate stages
  geom_rect(aes(xmax=467.3, xmin = 458.4, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=453, xmin = 445.2, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=443.8, xmin = 440.8, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=438.5, xmin = 433.4, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=430.5, xmin = 427.4, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=425.6, xmin = 423, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=419.2, xmin = 410.8, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=290.1, xmin = 283.5, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=273.01, xmin = 266.9, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=264.28, xmin = 259.51, ymin = 0, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=254.14, xmin = 251.902, ymin = 0, ymax = Inf),fill= "grey90") +
  
  # below adds the squares estimates to the plot
  geom_line(data=NA_data, aes(x = MID.POINT, y = Genus_count, colour="NA"),size = 1)+
  geom_line(data=Europe_data, aes(x = MID.POINT, y = Genus_count, colour="EU"),  size = 1) +
  geom_line(data=Asia_data, aes(x = MID.POINT, y = Genus_count, colour="Asia"),  size = 1) +
  geom_line(data=Africa_data, aes(x = MID.POINT, y = Genus_count, colour="Africa"),  size = 1) +
  geom_line(data=SA_data, aes(x = MID.POINT, y = Genus_count, colour="SA"),  size = 1) +
  geom_line(data=Ant_data, aes(x = MID.POINT, y = Genus_count, colour="Ant"),  size = 1) +
  geom_line(data=Aus_data, aes(x = MID.POINT, y = Genus_count, colour="Aus"),  size = 1) +
  
  # set up axes
  scale_colour_manual(name="Groups",values =c("NA" = "#000000", "EU"="#E69F00", "Asia"="#56B4E9", "Africa"="#009E73",
                                              "SA"="#F0E442", "Ant"="#0072B2", "Aus"="#D55E00"))+
  scale_x_reverse(expand=c(0,0),limits = c(467.3, 251.902),breaks = c(250,300,350,400,450)) + labs(x = "Ma", y = "Raw richness") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25, 50, 75, 100), limits = c(0, 100))

# Add time scale
rawplot_continents <- rawplot_continents+ coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
            dat = list("stages","periods"),
            height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))


#######################################################################
#        Squares analyses continents                                  #
#######################################################################

## This squares part is a modified version of the squares scripts by Allen et al. 2020 - "The latitudinal diversity gradient of tetrapods across the Permo-Triassic mass extinction and recovery interval"
## Diversity here is estimated using the squares method by Alroy (2018) 

# Create a vector giving the chronological order of stages

stages <- c("Darriwilian", "Sandbian", "Katian",  "Hirnantian","Rhuddanian", "Aeronian", "Telychian",
            "Sheinwoodian", "Homerian","Gorstian","Ludfordian","Pridoli","Lochkovian", "Pragian", "Emsian", "Eifelian","Givetian", 
            "Frasnian", "Famennian", "Tournaisian",  "Visean", "Serpukhovian", "Bashkirian",  "Moscovian", "Kasimovian", "Gzhelian", 
            "Asselian", "Sakmarian", "Artinskian", "Kungurian", "Roadian","Wordian","Capitanian", "Wuchiapingian","Changhsingian")

#Create a vector of stage midpoints
midpoints <-c(462.85,455.7,449.1,444.5,442.3,439.65,435.95,431.95,428.95,426.5,424.3,421.1,415,409.2, 400.45,390.5,385.2,377.45,365.55,
              352.8,338.8,327.05,319.2,311.1,305.35,301.3,296.21,291.81,286.8,278.225,270.875,266.95,262.1,256.62,253.021)


## First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

# North America
gen_freq1 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_NA %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq1[[i]] <- gen_list
}

names(gen_freq1) <- stages

# Europe

gen_freq2 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_Europe %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq2[[i]] <- gen_list
}

names(gen_freq2) <- stages

# Asia

gen_freq3 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_Asia%>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq3[[i]] <- gen_list
}

names(gen_freq3) <- stages

# Africa

gen_freq4 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_Africa %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq4[[i]] <- gen_list
}

names(gen_freq4) <- stages

# South America

gen_freq5 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_SA %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq5[[i]] <- gen_list
}

names(gen_freq5) <- stages

# Antarctica

gen_freq6 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_Ant %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq6[[i]] <- gen_list
}

names(gen_freq6) <- stages

# Australia

gen_freq7 <- list()

for (i in 1:length(stages)) {
  gen_list <- shark_Aus %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq7[[i]] <- gen_list
}

names(gen_freq7) <- stages


## Second: Calculate squares estimator from the list created above

# North America

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq1)) {
  freq_list <- gen_freq1[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_NA <- data.frame(midpoints, squares_list)

# Europe

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq2)) {
  freq_list <- gen_freq2[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_EU <- data.frame(midpoints, squares_list)

# Asia

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq3)) {
  freq_list <- gen_freq3[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_Asia <- data.frame(midpoints, squares_list)

# Africa

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq4)) {
  freq_list <- gen_freq4[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_Africa <- data.frame(midpoints, squares_list)

# South America

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq5)) {
  freq_list <- gen_freq5[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_SA <- data.frame(midpoints, squares_list)

# Antarctica

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq6)) {
  freq_list <- gen_freq6[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_Ant <- data.frame(midpoints, squares_list)

# Australia

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq7)) {
  freq_list <- gen_freq7[[i]]
  if(is.na(freq_list[1])){freq_list <- 0}
  if(freq_list[1] == 0){squares <- 0} else {
    sp_count <- length(freq_list)
    sing_count <- sum(freq_list == 1)
    ind_count <- sum(freq_list)
    sum_nsq <- sum(freq_list^2)
    squares <- sp_count + (((sing_count^2)*sum_nsq)/((ind_count^2) - (sing_count*sp_count)))
    if(squares == Inf){squares <- length(freq_list)}
  }
  squares_list <- append(squares_list, squares)
}

squares_Aus <- data.frame(midpoints, squares_list)

## Third: Plot the data using ggplot

squaresplot_cont <- ggplot() +
  #manually add grey bars to differentiate stages
  geom_rect(aes(xmax=467.3, xmin = 458.4, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=453, xmin = 445.2, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=443.8, xmin = 440.8, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=438.5, xmin = 433.4, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=430.5, xmin = 427.4, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=425.6, xmin = 423, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=419.2, xmin = 410.8, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=290.1, xmin = 283.5, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=273.01, xmin = 266.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=264.28, xmin = 259.51, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=254.14, xmin = 251.902, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  
  # below adds the squares estimates to the plot
  geom_line(data=squares_NA, aes(x = midpoints, y = squares_list,colour="NA"), size=1) +
  geom_line(data=squares_EU, aes(x = midpoints, y = squares_list,colour="EU"), size=1) +
  geom_line(data=squares_Asia, aes(x = midpoints, y = squares_list,colour="Asia"), size=1) +
  geom_line(data=squares_Africa, aes(x = midpoints, y = squares_list,colour="Africa"), size=1) +
  geom_line(data=squares_SA, aes(x = midpoints, y = squares_list,colour="SA"), size=1) +
  geom_line(data=squares_Ant, aes(x = midpoints, y = squares_list,colour="Ant"), size=1) +
  geom_line(data=squares_Aus, aes(x = midpoints, y = squares_list,colour="Aus"), size=1) +
  
  # set up axes, themes, margins, etc.:
  scale_colour_manual(name="Groups",values =c("NA" = "#000000", "EU"="#E69F00", "Asia"="#56B4E9", "Africa"="#009E73",
                                              "SA"="#F0E442", "Ant"="#0072B2", "Aus"="#D55E00"))+
  
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  labs(x = "Time (Ma)", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25,50,75,100,125), limits = c(0, 140)) +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))

squaresplot_cont

# Add time scale
squaresplot_cont <- squaresplot_cont + coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
                                 dat = list("stages","periods"),
                                 height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

##########################################################################################
############################## Spatial subsampling using divvy ############################
#########################################################################################

## 1. Prepare data for divvy

#Create a vector giving the chronological order of stages
stages_names <- intervals$interval_name

bins <- c(467.3, 458.4, 453.0,445.2,443.8,440.8,438.5, 433.4,430.5,427.4,425.6,423.0,419.2,410.8,407.6,
          393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,315.2,307.0,303.7,298.9,293.52,290.1,283.5,
          273.01,266.9,264.28,259.51,254.14,251.902) 

binDframe <- data.frame(bin = stages_names,
                        max_ma = as.numeric(bins[1:(length(bins)-1)]), 
                        min_ma = as.numeric(bins[2:(length(bins))]))

binDframe$mid_ma <- (binDframe$max_ma+binDframe$min_ma)/2

shark_data["new_bin"] <- NA
shark_data$new_bin <- as.numeric(as.character(shark_data$new_bin))


#=== Bin in any bin it falls within ===
All_Bin_List <- list()
for(s in 1:nrow(binDframe)){ # for each new bin
  temp_recs <- data.frame()
  for (o in 1:nrow(shark_data)){ 
    if (shark_data$MAX_DATE[o] > binDframe$min_ma[s] && shark_data$MIN_DATE[o] < binDframe$max_ma[s]){ # If occurrence max. age is greater than Bin min. age AND if occurrence min. age is less then Bin max. age (i.e. falls within bin at some point)
      temp_recs <- rbind(temp_recs, shark_data[o,]) # Add that occurrence to binlist
    }
  }
  if (nrow(temp_recs) > 0){
    temp_recs$new_bin <- s
  }
  All_Bin_List[[s]] <- temp_recs
}

shark_data_new <- do.call("rbind", All_Bin_List)

for (i in 1:length(stages_names)) {
  shark_data_new$new_bin[shark_data_new$new_bin== i] <- stages_names[i]}

##Tidy up and select only the columns necessary for the model input

names(shark_data_new)[names(shark_data_new) == 'new_bin'] <- 'stage'
names(shark_data_new)[names(shark_data_new) == 'MAX_DATE'] <- 'max_ma'
names(shark_data_new)[names(shark_data_new) == 'MIN_DATE'] <- 'min_ma'
names(shark_data_new)[names(shark_data_new) == 'MID.POINT'] <- 'mid_ma'

## 2. Calculate palaeocoordinates for chondrichthyans using palaeoverse

## Read in the dataset which includes the palaeo-rotated coordinates:
## For details regarding rotations of coordinates into palaeocoordinates, visit https://adamtkocsis.com/rgplates/index.html. 
sharks <- read.csv("./data/Palaeocoord.csv")

## 3. Divvy protocol as per Antell et al. 2023: https://gawainantell.github.io/divvy/ 

# Initialise Equal Earth projected coordinates
rWorld <- rast()
prj <- 'EPSG:8857'
rPrj <- project(rWorld, prj, res = 200000) # 200,000m is approximately 2 degrees
values(rPrj) <- 1:ncell(rPrj)

# Coordinate column names for the current and target coordinate reference system
xyCartes <- c('plng','plat')
xyCell   <- c('cellX','cellY')

# Extract cell number and centroid coordinates associated with each occurrence
llOccs <- vect(sharks, geom = xyCartes, crs = 'epsg:4326')
prjOccs <- project(llOccs, prj)
sharks$cell <- cells(rPrj, prjOccs)[,'cell']
sharks[, xyCell] <- xyFromCell(rPrj, sharks$cell)

## 4. Calculate metrics for global dataset

#Create a vector giving the chronological order of stages
stages_names <- intervals$interval_name

Stage_List <- list()
for(s in 1:length(stages_names)){
  temp <- data.frame()
  for(o in 1:nrow(sharks)){
    if(sharks$stage[o] == stages_names[s]){
      temp <- rbind(temp, sharks[o,]) 
    }
  }
  if (nrow(temp) > 0){
    temp$new_bin <- s
  }
  Stage_List[[s]] <- temp
}

names(Stage_List) <- stages_names

Stage_List <-within(Stage_List, rm(Katian))
Stage_List <-within(Stage_List, rm(Hirnantian))

summary_list <- list()

for(i in 1:length(Stage_List)){
  summary_list[[i]] <- sdSumry(Stage_List[i], taxVar = 'GENUS', xy = xyCell, crs = prj)
}

spatial_measures <- bind_rows(summary_list) 
names(spatial_measures)[names(spatial_measures) == "id"] <- "interval_name"

spatial_measures <- merge(spatial_measures, intervals, by="interval_name")
str(spatial_measures)

spatial_plot <- ggplot(spatial_measures) +
  geom_line(aes(x = MID.POINT, y = minSpanTree, colour = "#FDE725FF"), size = 1) +
  geom_line(aes(x = MID.POINT, y = greatCircDist,colour = "#29AF7FFF"), size = 1) +
  geom_line(aes(x = MID.POINT, y = meanPairDist, colour = "#39568CFF"), size = 1) +
  
  geom_point(aes(x = MID.POINT, y = minSpanTree, colour = "#FDE725FF"), pch=19, size = 2) +
  geom_point(aes(x = MID.POINT, y = greatCircDist,colour = "#29AF7FFF"), pch=19, size = 2) +
  geom_point(aes(x = MID.POINT, y = meanPairDist, colour = "#39568CFF"), pch=19, size = 2) +
  
  labs(x = "Time (Ma)", y = "Distance (km)") +
  scale_colour_manual(values= c("#FDE725FF","#29AF7FFF","#39568CFF"),
                      labels=c("summed MST length","GCD", "mean PD")) +
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  #scale_y_continuous(expand=c(0,0), breaks = c(0,25,50,75,100,125), limits = c(0, 125)) +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5),
        legend.title = element_blank(),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.02, vjust = -5))+
  #theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+
  ggtitle("Global")

spatial_plot

## 5. Calculate metrics for continent datasets

#Example: North America
shark_NA <- subset(sharks, CONTINENT=="North America")

#Create a vector giving the chronological order of stages
stages_names <- intervals$interval_name

Stage_List_NA <- list()
for(s in 1:length(stages_names)){
  temp <- data.frame()
  for(o in 1:nrow(shark_NA)){
    if(shark_NA$stage[o] == stages_names[s]){
      temp <- rbind(temp, shark_NA[o,]) 
    }
  }
  if (nrow(temp) > 0){
    temp$new_bin <- s
  }
  Stage_List_NA[[s]] <- temp
}

names(Stage_List_NA) <- stages_names

# Remove stages without data for divvy 
Stage_List_NA <-within(Stage_List_NA, rm(Darriwilian))
Stage_List_NA <-within(Stage_List_NA, rm(Katian))
Stage_List_NA <-within(Stage_List_NA, rm(Hirnantian))
Stage_List_NA <-within(Stage_List_NA, rm(Rhuddanian))
Stage_List_NA <-within(Stage_List_NA, rm(Aeronian))

summary_list_NA <- list()

for(i in 1:length(Stage_List_NA)){
  summary_list_NA[[i]] <- sdSumry(Stage_List_NA[i], taxVar = 'GENUS', xy = xyCell, crs = prj)
}

spatial_measures_NA <- bind_rows(summary_list_NA) 
names(spatial_measures_NA)[names(spatial_measures_NA) == "id"] <- "interval_name"

spatial_measures_NA <- merge(spatial_measures_NA, intervals, by="interval_name")

spatial_plot_NA <- ggplot(spatial_measures_NA) +
  geom_line(aes(x = MID.POINT, y = minSpanTree, colour = "#FDE725FF"), size = 1) +
  geom_line(aes(x = MID.POINT, y = greatCircDist,colour = "#29AF7FFF"), size = 1) +
  geom_line(aes(x = MID.POINT, y = meanPairDist, colour = "#39568CFF"), size = 1) +
  
  geom_point(aes(x = MID.POINT, y = minSpanTree, colour = "#FDE725FF"), pch=19, size = 2) +
  geom_point(aes(x = MID.POINT, y = greatCircDist,colour = "#29AF7FFF"), pch=19, size = 2) +
  geom_point(aes(x = MID.POINT, y = meanPairDist, colour = "#39568CFF"), pch=19, size = 2) +
  
  labs(x = "Time (Ma)", y = "Distance (km)") +
  scale_colour_manual(values= c("#FDE725FF","#29AF7FFF","#39568CFF"),
                      labels=c("summed MST length","GCD", "mean PD")) +
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  #scale_y_continuous(expand=c(0,0), breaks = c(0,25,50,75,100,125), limits = c(0, 125)) +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5),
        legend.title = element_blank(),
        plot.title = element_text(size = 10, face = "bold", hjust = 0.02, vjust = -5))+
  #theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+
  ggtitle("North America")

spatial_plot_NA


## Repeat for the other continents
## Plot together to compare patterns

######################################################End_of_script###################################################