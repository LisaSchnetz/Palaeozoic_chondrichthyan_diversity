###########################################################################
#                                                                         #
#           Chondrichthyan group SQS diversity analyses               #  
#                                                                         #
###########################################################################
#                                                                         #
#                      Lisa Schnetz - October 2023                      #
#                                                                         #
###########################################################################

## Packages used in this script:
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("deeptime")
devtools::install_version("iNEXT", version = "2.0.20")

library(tidyverse)
library(geoscale)
library(deeptime)
library(iNEXT)


intervals <- read.csv2("./data/iNEXTintervals2.csv")

#Load acanthodian and non-acanthodian chondrichthyans
shark_data <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")
shark_data <- subset(shark_data, select=c(GENUS, SPECIES, EARLIEST, LATEST, MAX_DATE, MIN_DATE, COLLECTION))

Acanthodians <-read.csv2("./data/Acanthodian_input_R_sampling.csv")
Acanthodians <- subset(Acanthodians, select=c(GENUS, SPECIES, EARLIEST, LATEST, MAX_DATE, MIN_DATE, COLLECTION))

# First: Create an incidence matrix (presence/absence matrix)
#Acanthodians
incidence_data_acan <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- Acanthodians %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_acan) <- intervals$interval_name 
incidence_data_acan <- incidence_data_acan[!sapply(incidence_data_acan, is.null)] 

### Need to remove stages with only one taxa
incidence_data_acan <-within(incidence_data_acan, rm(Wordian))
incidence_data_acan <-within(incidence_data_acan, rm(Capitanian))
incidence_data_acan <-within(incidence_data_acan, rm(Wuchiapingian))
incidence_data_acan <-within(incidence_data_acan, rm(Asselian))
incidence_data_acan <-within(incidence_data_acan, rm(Roadian))
incidence_data_acan <-within(incidence_data_acan, rm(Aeronian))
incidence_data_acan <-within(incidence_data_acan, rm(Telychian))
incidence_data_acan <-within(incidence_data_acan, rm(Kasimovian))
incidence_data_acan <-within(incidence_data_acan, rm(Gzhelian))

#Non-acanthodian chondrichthyans

incidence_data_chon <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- shark_data%>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_chon) <- intervals$interval_name 
incidence_data_chon <- incidence_data_chon[!sapply(incidence_data_chon, is.null)] 

### Need to remove stages with only one taxa
incidence_data_chon <-within(incidence_data_chon, rm(Gorstian))
incidence_data_chon <-within(incidence_data_chon, rm(Ludfordian))

# Second: Compute diversity estimates for coverage-based subsampling

## Create vector of quorum levels to run analysis at 
quorum_levels <- round(seq(from = 0.4, to = 0.7, by = 0.1), 1)

## Compute diversity estimate with estimateD()
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_acan, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_acan, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

# The output from the estimateD() analysis is a list object, so we'll turn it into a dataframe for plotting
plotting_data_acan <- bind_rows(estD_output) 

## Now join this up with data from the intervals dataset
plotting_data_acan <- plotting_data_acan %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 

## Filter down the data to only the quorum level(s) you want to plot:
plotting_data_acan <- filter(plotting_data_acan, quorum_level %in% quorum_levels[1:4])

## Ensure the interval_name and quorum_level columns are being treated as factors to avoid errors while plotting:
plotting_data_acan$interval_name <- as.factor(plotting_data_acan$interval_name) 
plotting_data_acan$quorum_level <- as.factor(plotting_data_acan$quorum_level)

## Compute diversity estimate with estimateD()
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_chon, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_chon, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

# The output from the estimateD() analysis is a list object, so we'll turn it into a dataframe for plotting

plotting_data_chon <- bind_rows(estD_output) 

## Now join this up with data from the intervals dataset
plotting_data_chon <- plotting_data_chon %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 

## Filter down the data to only the quorum level(s) you want to plot:
plotting_data_chon <- filter(plotting_data_chon, quorum_level %in% quorum_levels[1:4])

## Ensure the interval_name and quorum_level columns are being treated as factors to avoid errors while plotting:
plotting_data_chon$interval_name <- as.factor(plotting_data_chon$interval_name) 
plotting_data_chon$quorum_level <- as.factor(plotting_data_chon$quorum_level)


## Create a colour gradient if plotting more than one quorum level 
blue_gradient <- scales::seq_gradient_pal("cyan2","yellow", "Lab")(seq(0, 1, length.out = 4))

# Start plot

cov_rare_plot <- ggplot(plotting_data_acan, aes(x = MID.POINT, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  
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
  
  # the lines below add each quroum level to the plot and a colour from the colour gradient
  geom_ribbon(data=subset(plotting_data_acan, quorum_level == 0.4), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[1], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_acan, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[2], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_acan, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[3], alpha = 0.2) +
  
  # this sets line and point sizes/shapes/colours:
  geom_line(size = 1) +
  geom_point(aes(pch = Method), size = 3) +
  scale_shape_manual(values=c(15, 16, 17)) +
  scale_colour_manual(values = blue_gradient) +
  
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(450, 250), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,30,50,70), limits = c(0, 70)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot # check your plot

blue_gradient <- scales::seq_gradient_pal("darkblue","black", "Lab")(seq(0, 1, length.out = 5))

cov_rare_plot2 <- ggplot(plotting_data_chon, aes(x = MID.POINT, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
  
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
  
  # the lines below add each quroum level to the plot and a colour from the colour gradient
  geom_ribbon(data=subset(plotting_data_chon, quorum_level == 0.4), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[1], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_chon, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[2], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_chon, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[3], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_chon, quorum_level == 0.7), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[4], alpha = 0.2) +
  
  # this sets line and point sizes/shapes/colours:
  geom_line(size = 1) +
  geom_point(aes(pch = Method), size = 3) +
  scale_shape_manual(values=c(15, 16, 17)) +
  scale_colour_manual(values = blue_gradient) +
  
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(450, 250), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,30,50,70), limits = c(0, 70)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot2 # check your plot



#create plot using one quorum
plotting_data_chon <- filter(plotting_data_chon, quorum_level %in% quorum_levels[3])
plotting_data_acan <- filter(plotting_data_acan, quorum_level %in% quorum_levels[3])

new_plot <- ggplot()+
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
  geom_line(data=plotting_data_chon, aes(x = MID.POINT, y = qD, colour = "shark")) +
  geom_point(data=plotting_data_chon, aes(x = MID.POINT, y = qD, colour = "shark"),size=2) +
  geom_line(data=plotting_data_acan, aes(x = MID.POINT, y = qD,colour = "acan")) +
  geom_point(data=plotting_data_acan, aes(x = MID.POINT, y = qD,colour = "acan"),size=2) +
  
  scale_colour_manual(name=NULL, labels=c("Non-acanthodian chondrichthyans","Acanthodian chondrichthyans"),values =c("shark"="darkblue", "acan"="cyan2")) +
  
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(450, 250), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,30,50,70), limits = c(0, 70)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")
new_plot


#Subset to Elasmobranchii and Holocephali
#Read in dataset
shark_data <- read.csv2("./data/Total_chondrichthyes_input_groups.csv")

Elasmobranchii <- subset(shark_data, SUBCLASS=="Elasmobranchii")
Holocephali <- subset(shark_data, SUBCLASS=="Holocephali")

# First: Create an incidence matrix (presence/absence matrix)

incidence_data <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- genus_data %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data) <- intervals$interval_name 
incidence_data <- incidence_data[!sapply(incidence_data, is.null)] 

incidence_data[[1]] # check that it has worked correctly 


# Second: Compute diversity estimates for coverage-based subsampling

## Create vector of quorum levels to run analysis at 
quorum_levels <- round(seq(from = 0.3, to = 0.6, by = 0.1), 1)
quorum_levels <- round(seq(from = 0.5, to = 0.9, by = 0.1), 1)


## Compute diversity estimate with estimateD()
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 
