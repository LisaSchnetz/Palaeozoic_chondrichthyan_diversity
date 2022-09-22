####################################################
#                                                  #
#    Estimate sampling proxies from raw data        #  
#                                                  #
####################################################
#                                                  #
#         Lisa Schnetz - September 2022            #
#                                                  #
####################################################

## Packages used in this script:

install.packages("tidyverse")
install.packages("deeptime")
install.packages("ggplot2")
install.packages("dplyr")

library(tidyverse)
library(deeptime)
library(ggplot2)
library(dplyr)

# First make sure that your environment is clean so that you don't mix up data
rm(list=ls()) 

## Read in the datasets:
shark_data <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")

Acanthodians <-read.csv2("./data/Acanthodian_input_R_sampling.csv")
allsharks <- rbind(shark_data, Acanthodians)

intervals <- read.csv2("./data/iNEXTintervals.csv")

glimpse(allsharks)


#Count your data:

#Here, we are going to count raw numbers for the different units we want to consider. For sampling biases, 
#this should contain number of taxa, formations and localities (= collections)

# 1. Taxa per interval 

count_taxa <- vector("numeric") 
for (i in 1:nrow(intervals)) {
  out <- subset(allsharks, MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE)
  count_taxa[i] <- (length(unique(out$GENUS)))
  print(count_taxa[i])
}


# 2. Collections per interval
allsharks2 <- distinct(allsharks, GENUS, COLLECTION, .keep_all = TRUE)

count_colls <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(allsharks2, MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE)
  count_colls[i] <- (length(unique(out$COLLECTION)))
  print(count_colls[i])
}

# Formations per interval
allsharks3 <- distinct(allsharks, GENUS, FORMATION, .keep_all = TRUE)

count_formations <- vector("numeric")
for (i in 1:nrow(intervals)) {
  out <- subset(allsharks3, MAX_DATE > intervals[i,]$MIN_DATE & MIN_DATE < intervals[i,]$MAX_DATE)
  count_formations[i] <- (length(unique(out$FORMATION)))
  print(count_formations[i])
}


#Add all the information into a dataframe to make plotting easier####

proxy_counts <- data.frame(intervals$interval_name, intervals$MID.POINT, count_taxa, count_colls, count_formations)

# Rename the columns for ease:
proxy_counts <- rename(proxy_counts, 
                       "interval_name" = "intervals.interval_name", 
                       "mid_ma" = "intervals.MID.POINT")


## Finally, convert all zero's to NAs for plotting 
proxy_counts[proxy_counts == 0] <- NA 


#Plot the data 

colors <-c("Collections" ="orange", "Formations"="brown1","Genera"="black")


## Set your ggplot theme:
raw_plot <- ggplot(data=proxy_counts, aes(x=mid_ma)) + 
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
  
  # the lines below add the different proxies
  geom_line(data=proxy_counts, aes(y=count_colls, colour="Collections"), size=0.5) +
  geom_point(data=proxy_counts, aes(y=count_colls, colour="Collections"), size=2)+
  
  geom_line(data=proxy_counts, aes(y=count_formations, colour="Formations"), size=0.5) +
  geom_point(data=proxy_counts, aes(y=count_formations, colour="Formations"),size=2) +
  
  geom_line(data=proxy_counts, aes(y=count_taxa, colour="Genera"),size=0.5) +
  geom_point(data=proxy_counts, aes(y=count_taxa, colour="Genera"),size=2) +
 
   # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+
   scale_colour_manual(name="Sampling",values=c("orange","brown1","black"),labels=c("Collections","Formations","Genera"))+
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,50,100,150,200,250,300), limits = c(0, 300)) +
  labs(x = "Time (Ma)", y = "Total count")

raw_plot # check plot

Proxyplot <- gggeo_scale(raw_plot, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
Proxyplot <- gggeo_scale(Proxyplot, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = Proxyplot,
       width = 16, height = 13, dpi = 600, units = "cm", 
       filename = "./plots/sampling_proxies.pdf", useDingbats=FALSE)


# Lastly, perform a simple regression analysis to check correlation

## Raw diversity vs. collections
ggplot(proxy_counts, aes(x=count_taxa, y=count_colls)) + 
  geom_point(shape=17, size = 6, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1")  +
  theme_minimal()

## Raw diversity vs. formations
ggplot(proxy_counts, aes(x=count_taxa, y=count_formations)) + 
  geom_point(shape=16, size = 5, colour = "orange")+
  geom_smooth(method=lm, colour = "orange4", fill = "orange1") +
  theme_minimal()


## Let's quantify these relationships:
lm_colls = lm(count_colls ~ count_taxa, proxy_counts) # linear model
summary(lm_colls) # summary of results

lm_forms = lm(count_formations ~ count_taxa, proxy_counts) # linear model
summary(lm_forms) # summary of results


