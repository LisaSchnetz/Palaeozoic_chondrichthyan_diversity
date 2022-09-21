###########################################################################
#                                                                         #
#           Chondrichthyan group squares diversity analyses               #  
#                                                                         #
###########################################################################
#                                                                         #
#                      Lisa Schnetz - September 2022                      #
#                                                                         #
###########################################################################

## Packages used in this script:
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("deeptime")

library(tidyverse)
library(geoscale)
library(deeptime)

## Parts of this script are a modified version of the squares scripts by Allen et al. 2020 - "The latitudinal diversity gradient of tetrapods across the Permo-Triassic mass extinction and recovery interval"
## Diversity here is estimated using the squares method by Alroy (2018) 

#Create a vector giving the chronological order of stages

stages <- c("Darriwilian", "Sandbian", "Katian",  "Hirnantian","Rhuddanian", "Aeronian", "Telychian",
            "Sheinwoodian", "Homerian","Gorstian","Ludfordian","Pridoli","Lochkovian", "Pragian", "Emsian", "Eifelian","Givetian", 
            "Frasnian", "Famennian", "Tournaisian",  "Visean", "Serpukhovian", "Bashkirian",  "Moscovian", "Kasimovian", "Gzhelian", 
            "Asselian", "Sakmarian", "Artinskian", "Kungurian", "Roadian","Wordian","Capitanian", "Wuchiapingian","Changhsingian")

#Create a vector of stage midpoints
midpoints <-c(462.85,455.7,449.1,444.5,442.3,439.65,435.95,431.95,428.95,426.5,424.3,421.1,415,409.2, 400.45,390.5,385.2,377.45,365.55,
              352.8,338.8,327.05,319.2,311.1,305.35,301.3,296.21,291.81,286.8,278.225,270.875,266.95,262.1,256.62,253.021)


#Read in dataset
shark_data <- read.csv2("./data/Total_chondrichthyes_input_groups.csv")


intervals <- read.csv2("./data/iNEXTintervals.csv")


######################################################
#                Squares analyses                    #
######################################################

## Elasmobranchii
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

Elasmobranchii <- subset(shark_data, SUBCLASS=="Elasmobranchii")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Elasmobranchii %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_plot_elasmo <- data.frame(midpoints, squares_list)

## Holocephali
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

Holocephali <- subset(shark_data, SUBCLASS=="Holocephali")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Holocephali %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_plot_holo <- data.frame(midpoints, squares_list)


## Toothplates only
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr


Toothplates <- subset(shark_data, SUBCLASS=="Holocephali"& TOOTH_PLATE=="Y")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Toothplates %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_plot_plates <- data.frame(midpoints, squares_list)


## Stem
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

Stem <- subset(shark_data, SUBCLASS=="Stem")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Stem %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_plot_stem <- data.frame(midpoints, squares_list)


## Third: Plot the data using ggplot

squaresplot_subclass <- ggplot()+
  
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
  geom_line(data=to_plot_elasmo, aes(x = midpoints, y = squares_list, colour="Elasmo"),size = 1)+
  geom_line(data=to_plot_holo, aes(x = midpoints, y = squares_list, colour="Holo"),  size = 1) +
  geom_line(data=to_plot_plates, aes(x = midpoints, y = squares_list, colour="Plate"),  size = 1) +
  geom_line(data=to_plot_stem, aes(x = midpoints, y = squares_list, colour="Stem"),  size = 1) +
  
  # set up axes
  scale_colour_manual(name="Groups",values =c("Elasmo"="#FFA500","Holo"="#000000", "Plate"="grey", "Stem"="red"))+
  scale_x_reverse(expand=c(0,0),limits = c(467.3, 251.902),breaks = c(250,300,350,400,450)) + labs(x = "Ma", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25, 50, 75,100), limits = c(0, 100))

#add time scale to your plot
squaresplot_subclass <- gggeo_scale(squaresplot_subclass, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
squaresplot_subclass <- gggeo_scale(squaresplot_subclass, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)


ggsave(plot = squaresplot_subclass,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/squares_plot_subclass_all.pdf", useDingbats=FALSE)



############################################################################################
#                Squares analyses using the two major groups by Ginter et al. 2010         #
############################################################################################

## Euchondrocephali
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

Euchondrocephali2 <- subset(shark_data, GROUP_GINTER2=="Euchondrocephali")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Euchondrocephali2 %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_ploteuch2 <- data.frame(midpoints, squares_list)

## Elasmobranchii
# First: Generate a list of frequencies by stage by using and trimming the 'count' function in dplyr

Elasmobranchii2 <- subset(shark_data, GROUP_GINTER2=="Elasmobranchii")

gen_freq <- list()

for (i in 1:length(stages)) {
  gen_list <- Elasmobranchii2 %>% filter(MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE) %>%
    count(., GENUS) %>% arrange(desc(n)) %>% select(n)
  gen_list <- unlist(gen_list, use.names = F)
  gen_freq[[i]] <- gen_list
}

names(gen_freq) <- stages

## Second: Calculate squares estimator from the list created above

squares_list <- vector("numeric", length = 0)

for(i in 1:length(gen_freq)) {
  freq_list <- gen_freq[[i]]
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

to_plotelasmo2 <- data.frame(midpoints, squares_list)


## Third: Plot the data using ggplot

squaresplot_subclass_ginter <- ggplot()+
  
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
  geom_line(data=to_plotelasmo2, aes(x = midpoints, y = squares_list, colour="Elasmo"),size = 1)+
  geom_line(data=to_ploteuch2, aes(x = midpoints, y = squares_list, colour="Holo"),  size = 1) +
  
  # set up axes
  scale_colour_manual(name="Groups",values =c("Elasmo"="#FFA500","Holo"="#000000"))+
  scale_x_reverse(expand=c(0,0),limits = c(467.3, 251.902),breaks = c(250,300,350,400,450)) + labs(x = "Ma", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25, 50, 75), limits = c(0, 75))


#add time scale to your plot
squaresplot_subclass_ginter <- gggeo_scale(squaresplot_subclass_ginter, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
squaresplot_subclass_ginter <- gggeo_scale(squaresplot_subclass_ginter, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)


ggsave(plot = squaresplot_subclass_ginter,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/squares_plot_Ginter.pdf", useDingbats=FALSE)

