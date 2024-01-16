###########################################################################
#                                                                         #
#                       Spatial diversity analyses                        #  
#                                                                         #
###########################################################################
#                                                                         #
#                        Lisa Schnetz - May 2023                          #
#                                                                         #
###########################################################################

# Libraries used in this script
library(tidyverse)
library(dplyr)
library(ggplot2)
library(deeptime)
library(ggpubr)
library(viridis)

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

#add time scale to your plot
rawplot_continents <- rawplot_continents+ coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
            dat = list("stages","periods"),
            height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

ggsave(plot = rawplot_continents,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/spatial_biases_raw.pdf", useDingbats=FALSE)

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
#add time scale to your plot
squaresplot_cont <- squaresplot_cont + coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
                                 dat = list("stages","periods"),
                                 height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

ggsave(plot = squaresplot_cont,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/spatial_biases_squares.pdf", useDingbats=FALSE)

#######################
#
#
###### Plot Palaeocoordinates ##
#
#
#################
library(rgplates)
library(chronosphere)

# Get offline palaeomodel data

#Meredith et al. model
merdith2021<- chronosphere::fetch(src="EarthByte", ser="MERDITH2021")

shark_data <- read.csv("./data/Total_chondrichthyes_coordinates_new.csv")

collections <- unique(shark_data[, c("COLLECTION", "Longitude", "Latitude",  "MID.POINT")])
collections <- na.omit(collections) 

# Merdith et al one
paleoCoords2 <- rgplates::reconstruct(collections[, c("Longitude","Latitude")], age=collections$MID.POINT, model=merdith2021, enumerate=FALSE)
colnames(paleoCoords2) <- c("plng", "plat")

colls2 <- cbind(collections, paleoCoords2)

shark_data_merdith <- merge(shark_data, colls2, by=c("COLLECTION", "Longitude", "Latitude", "MID.POINT"))

## Save data as csv file
#write_csv(shark_data_merdith, "./data/merdith_coord_021023.csv")

coord_data <- read.csv("./data/merdith_coord_021023.csv")

## Plot coordinates on simple plot
# Count the number of taxa per collection (i.e. their frequency):
taxa_freqs_2 <- count(coord_data, COLLECTION)

taxa_freqs <- as.data.frame(table(coord_data[, c('COLLECTION','GENUS')]$COLLECTION))
names(taxa_freqs)[1] <- "COLLECTION" 
head(taxa_freqs) 




## Subset lat_data to only the columns we need:
coord_data <- coord_data%>% 
  select(COLLECTION, plat, plng, MID.POINT) %>% 
  distinct() %>% na.omit()

## Add add the frequency information:
coord_data <- left_join(taxa_freqs, coord_data, by = "COLLECTION")

## Before we plot, let's order the frequencies and remove any NAs that have crept in:
coord_data <- coord_data %>% arrange(Freq) %>% na.omit()


## Set up our ggplot layers
lat_plot <- ggplot(data = coord_data, aes(x = MID.POINT, y = plat, colour = Freq)) +
 # geom_vline(xintercept = int_boundaries, lty = 2, col = "grey90") +
  geom_hline(yintercept = 0, colour = "grey10") +
  scale_color_viridis(trans = "log", breaks = c(1, 2, 10, 95), direction = -1, option = "D") + # set the break= to match your richness data
  scale_y_continuous(labels = function(x) format(x, width = 5), limits = c(-90, 90), breaks = seq(from = -90, to = 90, by = 20)) +
  scale_x_reverse() + 
  theme_minimal() + 
  theme(legend.direction = "vertical", 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
        axis.title = element_text(size = 12)) +
  labs(x = "", y = "Paleolatitude (º)") +
  geom_point(size = 4, alpha = 0.5) # (alpha sets point transparency)
lat_plot # call to plot window


lat_plot_strat <- lat_plot + coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
                                                    dat = list("stages","periods"),
                                                    height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

lat_plot_strat




## Set dimensions and save plot (as pdf)
ggsave(plot = lat_plot,
       width = 20, height = 10, dpi = 500, units = "cm", 
       filename = "./plots/lat_alpha_div.pdf", useDingbats=FALSE)

