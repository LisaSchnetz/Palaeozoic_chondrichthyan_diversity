
library(ggplot2)
library(tidyverse)
library(deeptime)
library(viridis)
#devtools::install_github("willgearty/deeptime")
#Read in the datasets
shark_data <- read.csv2("./data/Chondrichthyes_input_groups.csv")
intervals <- read.csv2("./data/iNEXTintervals.csv")



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
elasmo_data$group <- "Elasmobranchii"

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
holo_data$group <- "Holocephali"

total <- rbind(elasmo_data, holo_data)

areaplot <- ggplot(total) +
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
  
  
  geom_area(aes(x=MID.POINT, y=Genus_count, fill = group),colour = "black") +
  
  # set up axes, themes, margins, etc.:
  #scale_fill_viridis(option="D", discrete=T) +
  scale_fill_brewer(name = "Subclass", palette = "YlOrRd") +
  #scale_fill_discrete(values=c("#FFA500","#000000")) +
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  labs(x = "Time (Ma)", y = "Sampled-in-bin richness") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,20,40), limits = c(0, 40)) +
  theme_classic()+
  theme(legend.position=c(0.12, 0.9),legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))


areaplot <-areaplot+ coord_geo(xlim = c(467.3, 251.902), ylim = c(0, 40), pos = as.list(rep("bottom",2)),
            dat = list("stages","periods"),
            height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))
areaplot

#Save a copy of the plot to your plots folder:
ggsave(plot = areaplot,
       width = 7.01, height = 5.96, dpi = 600, units = "in", 
       filename = "./stacked_area_subclasses.png")