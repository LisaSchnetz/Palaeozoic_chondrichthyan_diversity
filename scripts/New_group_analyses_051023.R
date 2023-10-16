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
#devtools::install_version("iNEXT", version = "2.0.20")

library(tidyverse)
library(geoscale)
library(deeptime)
library(iNEXT)


intervals <- read.csv2("./data/iNEXTintervals.csv")

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


# Start plot

cov_rare_plot <- ggplot() + 
  
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
  geom_ribbon(data=subset(plotting_data_acan, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "cyan2", alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data_chon, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "darkblue", alpha = 0.2) +
 
  # this sets line and point sizes/shapes/colours:
  geom_line(data=subset(plotting_data_acan, quorum_level == 0.6), aes(x = MID.POINT, y=qD, colour ="cyan2"), size = 1) +
  geom_line(data=subset(plotting_data_chon, quorum_level == 0.6), aes(x = MID.POINT, y=qD, colour ="darkblue"), size = 1) +
  
  geom_point(data=subset(plotting_data_acan, quorum_level == 0.6), aes(x = MID.POINT, y=qD, colour ="cyan2",pch = Method), size = 3) +
  geom_point(data=subset(plotting_data_chon, quorum_level == 0.6), aes(x = MID.POINT, y=qD, colour ="darkblue", pch = Method), size = 3) +
  scale_shape_manual(name = 'method',values=c(15, 16, 17)) +
  scale_colour_manual(name=NULL, labels=c("Acanthodian chondrichthyans","Non-acanthodian chondrichthyans"),values = c("cyan2","darkblue")) +
  
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20,30,40), limits = c(0, 40)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot # check your plot

cov_rare_plot <- cov_rare_plot + coord_geo(xlim = c(467.3, 251.902), ylim = c(0, 40), pos = as.list(rep("bottom",2)),
            dat = list("stages","periods"),
            height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

ggsave(plot = cov_rare_plot,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/chond_acan_SQS_plot.pdf", useDingbats=FALSE)

#####################
#####
####          Chondrichthyan groups
###
####################

shark_data <- read.csv2("./data/Total_chondrichthyes_input_groups.csv")

## Elasmobranchii
Elasmobranchii <- subset(shark_data, SUBCLASS=="Elasmobranchii")

incidence_data_elasmo <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- Elasmobranchii %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_elasmo) <- intervals$interval_name 
incidence_data_elasmo <- incidence_data_elasmo[!sapply(incidence_data_elasmo, is.null)] 

### Need to remove stages with only one taxa
incidence_data_elasmo <-within(incidence_data_elasmo, rm(Lochkovian))
incidence_data_elasmo <-within(incidence_data_elasmo, rm(Pragian))


## Holocephali
Holocephali <- subset(shark_data, SUBCLASS=="Holocephali")

incidence_data_holo <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- Holocephali %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_holo) <- intervals$interval_name 
incidence_data_holo <- incidence_data_holo[!sapply(incidence_data_holo, is.null)] 

### Need to remove stages with only one taxa
incidence_data_holo <-within(incidence_data_holo, rm(Eifelian))
incidence_data_holo <-within(incidence_data_holo, rm(Givetian))
incidence_data_holo <-within(incidence_data_holo, rm(Capitanian))

## Tooth plates only
Toothplates <- subset(shark_data, SUBCLASS=="Holocephali"& TOOTH_PLATE=="Y")

incidence_data_tooth <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- Toothplates %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_tooth) <- intervals$interval_name 
incidence_data_tooth <- incidence_data_tooth[!sapply(incidence_data_tooth, is.null)] 

### Need to remove stages with only one taxa

incidence_data_tooth <-within(incidence_data_tooth, rm(Frasnian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Famennian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Givetian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Asselian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Sakmarian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Artinskian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Kungurian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Roadian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Wuchiapingian))
incidence_data_tooth <-within(incidence_data_tooth, rm(Changhsingian))


## Stem only
Stem <- subset(shark_data, SUBCLASS=="Stem")

incidence_data_stem <- lapply(1:nrow(intervals), function(i) { 
  tmp_data <- Stem %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data_stem) <- intervals$interval_name 
incidence_data_stem <- incidence_data_stem[!sapply(incidence_data_stem, is.null)] 

### Need to remove stages with only one taxa
incidence_data_stem <-within(incidence_data_stem, rm(Sandbian))
incidence_data_stem <-within(incidence_data_stem, rm(Rhuddanian))
incidence_data_stem <-within(incidence_data_stem, rm(Kasimovian))
incidence_data_stem <-within(incidence_data_stem, rm(Gzhelian))
incidence_data_stem <-within(incidence_data_stem, rm(Asselian))
incidence_data_stem <-within(incidence_data_stem, rm(Artinskian))
incidence_data_stem <-within(incidence_data_stem, rm(Roadian))
incidence_data_stem <-within(incidence_data_stem, rm(Wordian))
incidence_data_stem <-within(incidence_data_stem, rm(Capitanian))
incidence_data_stem <-within(incidence_data_stem, rm(Wuchiapingian))

##########
## Compute diversity estimate with estimateD()
###########

#Elasmobranchii
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_elasmo, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_elasmo, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

plotting_elasmo <- bind_rows(estD_output) 
plotting_elasmo <- plotting_elasmo %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 
plotting_elasmo <- filter(plotting_elasmo, quorum_level %in% quorum_levels[1:4])
plotting_elasmo$interval_name <- as.factor(plotting_elasmo$interval_name) 
plotting_elasmo$quorum_level <- as.factor(plotting_elasmo$quorum_level)

#Holocephali
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_holo, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_holo, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

plotting_holo <- bind_rows(estD_output) 
plotting_holo <- plotting_holo %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 
plotting_holo <- filter(plotting_holo, quorum_level %in% quorum_levels[1:4])
plotting_holo$interval_name <- as.factor(plotting_holo$interval_name) 
plotting_holo$quorum_level <- as.factor(plotting_holo$quorum_level)

# Tooth plates
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_tooth, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_tooth, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

plotting_tooth <- bind_rows(estD_output) 
plotting_tooth <- plotting_tooth %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 
plotting_tooth <- filter(plotting_tooth, quorum_level %in% quorum_levels[1:4])
plotting_tooth$interval_name <- as.factor(plotting_tooth$interval_name) 
plotting_tooth$quorum_level <- as.factor(plotting_tooth$quorum_level)

# Stem
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data_stem, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data_stem, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

plotting_stem <- bind_rows(estD_output) 
plotting_stem <- plotting_stem %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 
plotting_stem <- filter(plotting_stem, quorum_level %in% quorum_levels[1:4])
plotting_stem$interval_name <- as.factor(plotting_stem$interval_name) 
plotting_stem$quorum_level <- as.factor(plotting_stem$quorum_level)

###########
#
########### Plot the data
#
###########

cov_rare_plot_groups <- ggplot() + 
  
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
  geom_ribbon(data=subset(plotting_elasmo, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "#FFA500", alpha = 0.2) +
  geom_ribbon(data=subset(plotting_holo, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "#000000", alpha = 0.2) +
  geom_ribbon(data=subset(plotting_tooth, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "grey", alpha = 0.2) +
  geom_ribbon(data=subset(plotting_stem, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = "red", alpha = 0.2) +
  
  # this sets line and point sizes/shapes/colours:
  geom_line(data=subset(plotting_elasmo, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="elasmo"), size = 1) +
  geom_line(data=subset(plotting_holo, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="holo"), size = 1) +
  geom_line(data=subset(plotting_tooth, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="tooth"), size = 1) +
  geom_line(data=subset(plotting_stem, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="stem"), size = 1) +
  
  geom_point(data=subset(plotting_elasmo, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="elasmo",pch = Method), size = 3) +
  geom_point(data=subset(plotting_holo, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="holo", pch = Method), size = 3) +
  geom_point(data=subset(plotting_tooth, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="tooth", pch = Method), size = 3) +
  geom_point(data=subset(plotting_stem, quorum_level == 0.5), aes(x = MID.POINT, y=qD, colour ="stem", pch = Method), size = 3) +
  
  scale_shape_manual(name = 'Method',values=c(15, 16, 17)) +
  scale_colour_manual(name="Groups",values =c("elasmo"="#FFA500","holo"="#000000", "tooth"="grey", "stem"="red"), labels=c("Elasmobranchii","Holocephali", "Tooth plates", "Stem-group")) +
  scale_fill_manual(values =c("elasmo"="#FFA500","holo"="#000000", "tooth"="grey", "stem"="red"))+
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,5,10,15,20,25), limits = c(0,25)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot_groups # check your plot

cov_rare_plot_groups <- cov_rare_plot_groups + coord_geo(xlim = c(467.3, 251.902), pos = as.list(rep("bottom",2)),
                                           dat = list("stages","periods"),
                                           height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))

ggsave(plot = cov_rare_plot_groups,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/group_SQS_plot.pdf", useDingbats=FALSE)
