###########################################################################
#                                                                         #
#                           iNext analyses                                #  
#                                                                         #
###########################################################################
#                                                                         #
#                      Lisa Schnetz - Oktober 2023                        #
#                                                                         #
###########################################################################

##This script primarily uses the package iNEXT implemented by Hsieh et al. 2016 | https://cran.r-project.org/package=iNEXT

## Packages used in this script:
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("deeptime")
install.packages("iNEXT")
install.packages("directlabels")
install.packages('remotes')

## Load packages:
library(iNEXT)
library(ggplot2)
library(tidyverse)
library(deeptime)
library(directlabels)
library(remotes)

# First make sure that your environment is clean so that you don't mix up data
rm(list=ls()) 

## Read in the datasets:
shark_data <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")

#Acanthodian data from Schnetz et al. 2022 - Palaeontology: https://doi.org/10.1111/pala.12616
Acanthodians <-read.csv2("./data/Acanthodian_input_R_sampling.csv")

##This interval data does not contain the Darriwilian, Sandbian, Katian and Hirnantian stages of the Ordovician 
# as they would lead to error messages in any of the iNEXT analyses. Due to the limited number of occurrences in those stages, 
# they do not add significant information to the diversity analyses and can be removed.

intervals <- read.csv2("./data/iNEXTintervals.csv")


### Combine the data for total chondrichthyan diversity and remove unnecessary columns
allsharks <- rbind(shark_data, Acanthodians)
genus_data <- subset(allsharks, select=c(GENUS, SPECIES, EARLIEST, LATEST, MAX_DATE, MIN_DATE, COLLECTION))

## The code below has been modified but is primarily taken from Dunne et al. 2018 - "Diversity change during the rise of tetrapods and the impact of the ‘Carboniferous rainforest collapse"
## doi:10.1098/rspb.2017.2730

## 1. Coverage based subsampling/SQS of total chondrichtyans using iNEXT

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
quorum_levels <- round(seq(from = 0.5, to = 0.8, by = 0.1), 1)

## Compute diversity estimate with estimateD()
## Might show an error when trying for the first time, check that individual quorums are working and try again.
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$Order.q == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 

# Third: Plot your data
# The output from the estimateD() analysis is a list object, so we'll turn it into a dataframe for plotting
plotting_data <- bind_rows(estD_output) 

## Now join this up with data from the intervals dataset
plotting_data <- plotting_data %>% rename(interval_name = Assemblage) %>% full_join(.,intervals, by = "interval_name") 

## Ensure the interval_name and quorum_level columns are being treated as factors to avoid errors while plotting:
plotting_data$interval_name <- as.factor(plotting_data$interval_name) 
plotting_data$quorum_level <- as.factor(plotting_data$quorum_level)

## Create a colour gradient if plotting more than one quorum level 
blue_gradient <- scales::seq_gradient_pal("cyan2", "darkslategrey", "Lab")(seq(0, 1, length.out = 4))

# Start plot
cov_rare_plot <- ggplot(plotting_data, aes(x = MID.POINT, y = qD, ymin = qD.LCL, ymax = qD.UCL, colour = quorum_level)) + 
 
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
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[1], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[2], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.7), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[3], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.8), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[4], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.9), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[4], alpha = 0.2) +
  
  # this sets line and point sizes/shapes/colours:
  geom_line(size = 1) +
  geom_point(aes(pch = Method), size = 3) +
  scale_shape_manual(values=c(15, 16, 17)) +
  scale_colour_manual(values = blue_gradient) +
  
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,30,50,70,90), limits = c(0, 90)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot # check your plot

#add time scale to your plot
cov_rare_plot <- cov_rare_plot + coord_geo(xlim = c(470, 250), pos = as.list(rep("bottom",2)),
                                  dat = list("stages","periods"),
                                  height = list(unit(1.5, "lines"),unit(1.5,"lines")), rot = list(0,0), size = list(2.5, 2.5), abbrv = list(TRUE, FALSE))
cov_rare_plot

## 2. Coverage-based rarefaction for each interval by period     
# Assemble data again but with abbreviated stages (to make plotting later easier)

incidence_data <- lapply(1:nrow(intervals), function(i) { # begin loop
  tmp_data <- genus_data %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) 
  if (nrow(tmp_data) > 1) {
    incidence_raw <- tmp_data %>% .[,c("COLLECTION","GENUS")] %>% distinct %>% table %>% t 
    incidence_raw <- as.data.frame.matrix(incidence_raw) 
    incidence_raw[incidence_raw > 1] <- rep(1, length(incidence_raw[incidence_raw > 1])) 
    incidence_raw 
  }
})
names(incidence_data) <- intervals$interval_name_abb # name each list item from the intervals data
incidence_data <- incidence_data[!sapply(incidence_data, is.null)] 

#Run analysis by applying the iNEXT function of the package
#You might have to run it more than once as it sometimes gives you an error. 

inc.data <- iNEXT(incidence_data, q = 0, datatype = "incidence_raw") 

cov_rare <- inc.data$iNextEst

for(i in 1:length(cov_rare)) {
  
  cov_rare[[i]]$stage_int <- names(cov_rare)[i]
  
}

cov_rare_size <- cov_rare$size_based %>% as_tibble() #convert to tibble for ease of plotting

#Add a geological period column to make plotting easier!
cov_rare_size[which(cov_rare_size$Assemblage %in% intervals$interval_name_abb[1:8]), "Period"] <- "Silurian"
cov_rare_size[which(cov_rare_size$Assemblage %in% intervals$interval_name_abb[9:15]), "Period"] <- "Devonian"
cov_rare_size[which(cov_rare_size$Assemblage %in% intervals$interval_name_abb[16:22]), "Period"] <- "Carboniferous"
cov_rare_size[which(cov_rare_size$Assemblage %in% intervals$interval_name_abb[23:31]), "Period"] <- "Permian"

cov_rare_cov <- cov_rare$coverage_based %>% as_tibble() #convert to tibble for ease of plotting

#Add a geological period column to make plotting easier!
cov_rare_cov[which(cov_rare_cov$Assemblage %in% intervals$interval_name_abb[1:8]), "Period"] <- "Silurian"
cov_rare_cov[which(cov_rare_cov$Assemblage %in% intervals$interval_name_abb[9:15]), "Period"] <- "Devonian"
cov_rare_cov[which(cov_rare_cov$Assemblage %in% intervals$interval_name_abb[16:22]), "Period"] <- "Carboniferous"
cov_rare_cov[which(cov_rare_cov$Assemblage %in% intervals$interval_name_abb[23:31]), "Period"] <- "Permian"


## Plot data. We will divide the data into two separate plots to make visualization easier:
## One for Silurian-Devonian, one for Carboniferous-Permian

colour_scheme <- c("#F04028","#F04028","#67A599","#F04028", "#F04028","#67A599","#67A599","#F04028",
                   "#67A599", "#F04028","#F04028","#67A599","#67A599","#67A599","#F04028","#F04028")

names(colour_scheme) <- c("Ar", "As", "Ba", "Ca", "Cha", "Gz","Ka","Ku","Mo","Ro", "Sa","Se","To", "Vi","Wo","Wu")
        
cov_rare1 <-subset(cov_rare_cov, Period=="Carboniferous" | Period=="Permian")

cov_rare_plot1 <- ggplot(data = cov_rare1, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage, colour = Period, lty = Method)) + 
  geom_line(size = 1) + 
  geom_smooth(aes(fill=Assemblage))+
  scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  scale_colour_manual(values = c("#67A599","#F04028")) +
  scale_fill_manual(values = colour_scheme) +
  theme(panel.background = element_blank(),
        legend.position="none",
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=14)) + 
  labs(x = "Sample coverage", y = "Genus diversity") +
  scale_x_continuous(limits = c(0, 1.1), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 130), expand=c(0,0), breaks = seq(0, 130, 30))
cov_rare_plot1
cov_rare_plot1 <- cov_rare_plot1 +geom_dl(data=cov_rare1, aes(label=Assemblage),method=list("last.points",rot=30))


cov_rare2 <-subset(cov_rare_cov, Period=="Silurian" | Period=="Devonian")

colour_scheme2 <- c("#B3E1B6","#B3E1B6","#B3E1B6","#B3E1B6","#B3E1B6","#B3E1B6","#B3E1B6","#B3E1B6","#CB8C37",
                    "#CB8C37","#CB8C37","#CB8C37","#CB8C37","#CB8C37","#CB8C37")

names(colour_scheme2) <- c("Rh",  "Ae",  "Te",  "Sh",  "Ho",  "Go",  "Lu",  "Pr",  "Lo",  "Pra", "Em",  "Ei",  "Gi",  "Fra", "Fam")

  
cov_rare_plot2 <- ggplot(data = cov_rare2, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage, colour = Period, lty = Method)) + 
  geom_line(size = 1) + 
  geom_smooth(aes(fill=Assemblage))+
  scale_fill_manual(values = colour_scheme2)+
  scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  scale_colour_manual(values = c("#CB8C37","#B3E1B6")) +
  # geom_point(data = cov_rare1, aes(x = SC, y = qD, pch = Period, colour = Period), size = 3, inherit.aes = F) + 
  theme(panel.background = element_blank(),
        legend.position="none",
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=14)) + 
  labs(x = "Sample coverage", y = "Genus diversity") +
  scale_x_continuous(limits = c(0, 1.1), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 130), expand=c(0,0), breaks = seq(0, 130, 30))
cov_rare_plot2
cov_rare_plot2 <- cov_rare_plot2 +geom_dl(data=cov_rare2, aes(label=Assemblage),method=list("last.points",rot=30))

#############
#########################Plot standard rarefaction curves###############
#############
cov_rare3 <-subset(cov_rare_size, Period=="Carboniferous" | Period=="Permian")

cov_rare_plot3 <- ggplot(data = cov_rare3, aes(x = t, y = qD, fill = Assemblage, colour = Period, lty = Method)) + 
  geom_line(size = 1) + 
  geom_smooth(aes(fill=Assemblage))+
  scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  scale_colour_manual(values = c("#67A599","#F04028")) +
  scale_fill_manual(values = colour_scheme) +
  theme(panel.background = element_blank(),
        legend.position="none",
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=14)) + 
  labs(x = "Number of occurrences", y = "Genus diversity") +
scale_x_continuous(limits = c(0, 650), expand=c(0,0), breaks = seq(0, 600, 100)) +
  scale_y_continuous(limits = c(0, 130), expand=c(0,0), breaks = seq(0, 130, 30))
cov_rare_plot3
cov_rare_plot3 <- cov_rare_plot3 +geom_dl(data=cov_rare3, aes(label=Assemblage),method=list("last.points",rot=30))

#Second plot
cov_rare4 <-subset(cov_rare_cov, Period=="Silurian" | Period=="Devonian")

cov_rare_plot4 <- ggplot(data = cov_rare4, aes(x = t, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = Assemblage, colour = Period, lty = Method)) + 
  geom_line(size = 1) + 
  geom_smooth(aes(fill=Assemblage))+
  scale_linetype_manual(values=c("dotted", "dotdash", "solid")) +
  scale_colour_manual(values = c("#CB8C37","#B3E1B6")) +
  scale_fill_manual(values = colour_scheme2) +
  # geom_point(data = cov_rare1, aes(x = SC, y = qD, pch = Period, colour = Period), size = 3, inherit.aes = F) + 
  theme(panel.background = element_blank(),
        legend.position="none",
        panel.grid.minor.y = element_line(colour = "grey90"),
        panel.grid.minor.x = element_line(colour = "grey90"),
        panel.grid.major.y = element_line(colour = "grey90"),
        panel.grid.major.x = element_line(colour = "grey90"),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(size=12, angle=0, hjust=0.5),
        axis.text.y = element_text(size=16),
        axis.title = element_text(size=14)) + 
  labs(x = "Number of occurrences", y = "Genus diversity") + 
  scale_x_continuous(limits = c(0, 650), expand=c(0,0), breaks = seq(0, 600, 100)) +
 scale_y_continuous(limits = c(0, 130), expand=c(0,0), breaks = seq(0, 130, 30))
cov_rare_plot4
cov_rare_plot4 <- cov_rare_plot4 +geom_dl(data=cov_rare4, aes(label=Assemblage),method=list("last.points",rot=30))

all<- ggarrange(cov_rare_plot1,cov_rare_plot3,cov_rare_plot2, cov_rare_plot4, nrow=2, ncol=2, labels = c('A', 'C','B','D','E'))

######################################################End_of_script###################################################