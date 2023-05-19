###########################################################################
#                                                                         #
#               Chondrichthyan group raw diversity analyses               #  
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


#Read in dataset
shark_data <- read.csv2("./data/Total_chondrichthyes_input_groups.csv")


intervals <- read.csv2("./data/Intervals.csv")


#######################################################################
#           Group divisions based on recent phylogenies               #
#   (Coates et al. 2018, Dearden et al. 2019 and Frey et al. 2020)    # 
#######################################################################

## First: count the number of genera for each group

## Elasmobranchii
Elasmobranchii <- subset(shark_data, SUBCLASS=="Elasmobranchii")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Elasmobranchii %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

elasmo_data<-as.data.frame(t(genus_data))
elasmo_data<-cbind(elasmo_data, intervals)


## Holocephali

Holocephali <- subset(shark_data, SUBCLASS=="Holocephali")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Holocephali %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

holo_data<-as.data.frame(t(genus_data))
holo_data<-cbind(holo_data, intervals)


## Tooth plates only

Toothplates <- subset(shark_data, SUBCLASS=="Holocephali"& TOOTH_PLATE=="Y")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Toothplates %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

plate_data<-as.data.frame(t(genus_data))
plate_data<-cbind(plate_data, intervals)


## Stem only

Stem <- subset(shark_data, SUBCLASS=="Stem")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Stem %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

stem_data<-as.data.frame(t(genus_data))
stem_data<-cbind(stem_data, intervals)


## Second: Plot the data

rawplot_subclass <- ggplot()+
  
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
  geom_line(data=elasmo_data, aes(x = MID.POINT, y = Genus_count, colour="Elasmo"),size = 1)+
  geom_line(data=holo_data, aes(x = MID.POINT, y = Genus_count, colour="Holo"),  size = 1) +
  geom_line(data=plate_data, aes(x = MID.POINT, y = Genus_count, colour="Plate"),  size = 1) +
  geom_line(data=stem_data, aes(x = MID.POINT, y = Genus_count, colour="Stem"),  size = 1) +
  
  # set up axes
  scale_colour_manual(name="Groups",values =c("Elasmo"="#FFA500","Holo"="#000000", "Plate"="grey", "Stem"="red"))+
  scale_x_reverse(expand=c(0,0),limits = c(467.3, 251.902),breaks = c(250,300,350,400,450)) + labs(x = "Ma", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25, 50, 75), limits = c(0, 75))

#add time scale to your plot
rawplot_subclass <- gggeo_scale(rawplot_subclass, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
rawplot_subclass <- gggeo_scale(rawplot_subclass, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)


ggsave(plot = rawplot_subclass,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/squares_plot_subclass_raw.pdf", useDingbats=FALSE)


############################################################################################
#                 Raw analyses using the two major groups by Ginter et al. 2010            #
############################################################################################

## First: count the number of genera for each group

## Elasmobranchii

Elasmobranchii2 <- subset(shark_data, GROUP_GINTER2=="Elasmobranchii")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Elasmobranchii2 %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

elasmo_data2<-as.data.frame(t(genus_data))
elasmo_data2<-cbind(elasmo_data2, intervals)

## Euchondrocephali

Euchondrocephali2 <- subset(shark_data, GROUP_GINTER2=="Euchondrocephali")

genus_interval <- lapply(1:nrow(intervals), function(i) {
  chosen_interval <- Euchondrocephali2 %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # filter the data to a single interval
  no_colls <- (length(unique(chosen_interval$GENUS))) # count the total no. of collections in that interval
  no_colls
})
names(genus_interval) <- intervals$interval_name # give each list element its correct interval name

genus_data <- bind_rows(genus_interval)
rownames(genus_data) = "Genus_count"

holo_data2<-as.data.frame(t(genus_data))
holo_data2<-cbind(holo_data2, intervals)


## Second: Plot the data

rawplot_subclass_ginter <- ggplot()+
  
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
  geom_line(data=elasmo_data2, aes(x = MID.POINT, y = Genus_count, colour="Elasmo"),size = 1)+
  geom_line(data=holo_data2, aes(x = MID.POINT, y = Genus_count, colour="Holo"),  size = 1) +
  
  # set up axes
  scale_colour_manual(name="Groups",values =c("Elasmo"="#FFA500","Holo"="#000000"))+
  scale_x_reverse(expand=c(0,0),limits = c(467.3, 251.902),breaks = c(250,300,350,400,450)) + labs(x = "Ma", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25, 50, 75), limits = c(0, 75))

#add time scale to your plot
rawplot_subclass_ginter <- gggeo_scale(rawplot_subclass_ginter, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
rawplot_subclass_ginter <- gggeo_scale(rawplot_subclass_ginter, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)


ggsave(plot = rawplot_subclass_ginter,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/plot_Ginter_raw.pdf", useDingbats=FALSE)

