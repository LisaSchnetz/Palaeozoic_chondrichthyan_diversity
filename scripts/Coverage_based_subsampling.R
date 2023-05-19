####################################################
#                                                  #
#               iNext analyses                     #  
#                                                  #
####################################################
#                                                  #
#         Lisa Schnetz - September 2022            #
#                                                  #
####################################################


##This script primarily uses the package iNEXT implemented by Hsieh et al. 2016 | https://cran.r-project.org/package=iNEXT

## Packages used in this script:

install.packages("ggplot2")
install.packages("tidyverse")
install.packages("deeptime")
install.packages("iNEXT")
install.packages("directlabels")
install.packages('remotes')

library(iNEXT)
library(ggplot2)
library(tidyverse)
library(deeptime)
library(directlabels)
library(remotes)

#Package version of iNEXT that was used for the below analyses
remotes::install_version("iNEXT", version = "2.0.20")

# First make sure that your environment is clean so that you don't mix up data
rm(list=ls()) 

## Read in the datasets:
shark_data <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")

Acanthodians <-read.csv2("./data/Acanthodian_input_R_sampling.csv")


##This interval data does not contain the Darriwilian, Sandbian, Katian and Hirnantian stages of the Ordovician 
# as they would lead to error messages in any of the iNEXT analyses. Due to the limited number of occurrences in those stages, 
# they do not add significant information to the diversity analyses and can be removed.

intervals <- read.csv2("./data/iNEXTintervals.csv")


### Combine the data for total chondrichthyan diversity and remove unnecessary columns

allsharks <- rbind(shark_data, Acanthodians)

genus_data <- subset(allsharks, select=c(GENUS, SPECIES, EARLIEST, LATEST, MAX_DATE, MIN_DATE, COLLECTION))


glimpse(genus_data)

#############################################################################
##                      Abundance analyses                                  ##
##############################################################################

## The code below has been modified but is primarily taken from Dunne et al. 2018 - "Diversity change during the rise of tetrapods and the impact of the ‘Carboniferous rainforest collapse"


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

## Compute diversity estimate with estimateD()
estD_output <- lapply(1:length(quorum_levels), function(i) { # loop will run over each quorum level set above
  estD_output <- estimateD(incidence_data, datatype="incidence_raw", base="coverage", level=quorum_levels[i]) # main estimateD code - see vignette for more details
  estD_output <- estD_output[estD_output$order == 0, ] # filter to richness (order == 0)
  estD_output$reference_t <- sapply(incidence_data, sum) # tally total occurrences in each bin
  estD_output[which(estD_output$t >= 2 * estD_output$reference_t), c("qD","qD.LCL","qD.UCL")] <- rep(NA, 3) # no more than twice ref sample size
  estD_output$quorum_level <- quorum_levels[i] # create new column in output
  estD_output # returns object
}) 


# Third: Plot your data


# The output from the estimateD() analysis is a list object, so we'll turn it into a dataframe for plotting

plotting_data <- bind_rows(estD_output) 

## Now join this up with data from the intervals dataset
plotting_data <- plotting_data %>% rename(interval_name = site) %>% full_join(.,intervals, by = "interval_name") 

## Filter down the data to only the quorum level(s) you want to plot:
plotting_data <- filter(plotting_data, quorum_level %in% quorum_levels[1:4])

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
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.4), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[1], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.5), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[2], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.6), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[3], alpha = 0.2) +
  geom_ribbon(data=subset(plotting_data, quorum_level == 0.7), aes(x = MID.POINT, ymin = qD.LCL, ymax = qD.UCL), inherit.aes = FALSE, fill = blue_gradient[4], alpha = 0.2) +
  
  # this sets line and point sizes/shapes/colours:
  geom_line(size = 1) +
  geom_point(aes(pch = method), size = 4.5) +
  scale_shape_manual(values=c(15, 16, 17)) +
  scale_colour_manual(values = blue_gradient) +
  
  # set up axes, themes, margins, etc.:
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+ 
  scale_x_reverse(expand=c(0,0), limits = c(450, 250), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,10,20,30,40), limits = c(0, 40)) +
  labs(x = "Time (Ma)", y = "Coverage rarified genus richness")

cov_rare_plot # check your plot

# Add a time scale to your plot
cov_rare_plot <- gggeo_scale(cov_rare_plot, dat = "stages", height = unit(4, "lines"), rot = 90, size = 2.5, abbrv = FALSE)
cov_rare_plot

# Save a copy of the plot to your plots folder:
ggsave(plot = cov_rare_plot,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./Coverage_subsampling_plot.pdf", useDingbats=FALSE)



########################################################################
##          Coverage-based rarefaction for each interval by period     ##
########################################################################

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

inc.data <- iNEXT(incidence_data, q = 0, datatype = "incidence_raw") 

cov_rare <- inc.data$iNextEst

for(i in 1:length(cov_rare)) {
  
  cov_rare[[i]]$stage_int <- names(cov_rare)[i]
  
}

cov_rare <- do.call(rbind, cov_rare) %>% as_tibble() #convert to tibble for ease of plotting

cov_rare[which(cov_rare$stage_int %in% intervals$interval_name_abb[1:8]), "Period"] <- "Silurian"
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name_abb[9:15]), "Period"] <- "Devonian"
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name_abb[16:22]), "Period"] <- "Carboniferous"
cov_rare[which(cov_rare$stage_int %in% intervals$interval_name_abb[23:31]), "Period"] <- "Permian"


# Plot data. We will divide data into two separate plots to make visualisation easier:
# One for Silurian-Devonian, one for Carboniferous-Permian


cov_rare1 <-subset(cov_rare, Period=="Carboniferous" | Period=="Permian")

cov_rare_plot1 <- ggplot(data = cov_rare1, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = stage_int, colour = Period, lty = method)) + 
  geom_line(size = 1) + 
  scale_linetype_manual(values=c("dotted", "solid", "dotdash")) +
  scale_colour_manual(values = c("#67A599","#F04028")) +
  #geom_point(data = cov_rare, aes(x = SC, y = qD, pch = Period, colour = Period), size = 3, inherit.aes = F) + 
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
  labs(x = "Coverage", y = "Species richness") +
  scale_x_continuous(limits = c(0, 1.05), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 150), expand=c(0,0), breaks = seq(0, 210, 30))
cov_rare_plot1
cov_rare_plot1 +geom_dl(data=cov_rare1, aes(label=stage_int),method=list("last.points",rot=30))

#Save a copy of the plot to your plots folder:
  ggsave(plot = cov_rare_plot1,
         width = 20, height = 15, dpi = 600, units = "cm", 
         filename = "./Coverage_rarefaction_plot_Carb_Perm.pdf", useDingbats=FALSE)

#Second plot
  
cov_rare2 <-subset(cov_rare, Period=="Silurian" | Period=="Devonian")

cov_rare_plot2 <- ggplot(data = cov_rare2, aes(x = SC, y = qD, ymin = qD.LCL, ymax = qD.UCL, fill = stage_int, colour = Period, lty = method)) + 
  geom_line(size = 1) + 
  scale_linetype_manual(values=c("dotted", "solid", "dotdash")) +
  scale_colour_manual(values = c("#CB8C37","#B3E1B6")) +
  #geom_point(data = cov_rare, aes(x = SC, y = qD, pch = Period, colour = Period), size = 3, inherit.aes = F) + 
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
  labs(x = "Coverage", y = "Species richness") +
  scale_x_continuous(limits = c(0, 1.05), expand=c(0,0), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 150), expand=c(0,0), breaks = seq(0, 210, 30))
cov_rare_plot2
cov_rare_plot2 +geom_dl(data=cov_rare2, aes(label=stage_int),method=list("last.points",rot=30))

#Save a copy of the plot to your plots folder:
ggsave(plot = cov_rare_plot2,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./Coverage_rarefaction_plot_Sil_Devon.pdf", useDingbats=FALSE)


########################################################################################
##            Create relative abundances for each time interval from lists above       ##
########################################################################################

#To check if stages have similar relative abundances or if there are big changes between stages

##Generate frequency data

genus_freq1 <- lapply(1:nrow(intervals), function(i) {
  tmp <- genus_data %>% filter(MAX_DATE > intervals[i,"MIN_DATE"] & MIN_DATE < intervals[i,"MAX_DATE"]) # # filter the data to a single interval
  gen_counts <- tmp %>% count(GENUS) %>% distinct
  freq_raw <- as.numeric(gen_counts$n)
  freq_raw
})
names(genus_freq1) <- intervals$interval_name # give each list element its correct interval name


#Silurian

df <- genus_freq1$Rhuddanian[order(genus_freq1$Rhuddanian,decreasing = TRUE)]
df1 <- genus_freq1$Aeronian[order(genus_freq1$Aeronian,decreasing = TRUE)]
df2 <- genus_freq1$Telychian[order(genus_freq1$Telychian,decreasing = TRUE)]
df3 <- genus_freq1$Sheinwoodian[order(genus_freq1$Sheinwoodian,decreasing = TRUE)]
df4 <- genus_freq1$Homerian[order(genus_freq1$Homerian,decreasing = TRUE)]
df5 <- genus_freq1$Gorstian[order(genus_freq1$Gorstian,decreasing = TRUE)]
df6 <- genus_freq1$Ludfordian[order(genus_freq1$Ludfordian,decreasing = TRUE)]
df7 <- genus_freq1$Pridoli[order(genus_freq1$Pridoli,decreasing = TRUE)]


par(mfrow=c(2,4))
barplot(df,ylab="Frequency", xlab="Number of genera", main="Rhuddanian",col="#CB8C37")
barplot(df1,ylab="Frequency", xlab="Number of genera", main="Aeronian",col="#CB8C37")
barplot(df2,ylab="Frequency", xlab="Number of genera", main="Telychian",col="#CB8C37")
barplot(df3,ylab="Frequency", xlab="Number of genera", main="Sheinwoodian",col="#CB8C37")
barplot(df4,ylab="Frequency", xlab="Number of genera", main="Homerian",col="#CB8C37")
barplot(df5,ylab="Frequency", xlab="Number of genera", main="Gorstian",col="#CB8C37")
barplot(df6,ylab="Frequency", xlab="Number of genera", main="Ludfordian",col="#CB8C37")
barplot(df7,ylab="Frequency", xlab="Number of genera", main="Pridoli",col="#CB8C37")

#Devonian

df8 <- genus_freq1$Lochkovian[order(genus_freq1$Lochkovian,decreasing = TRUE)]
df9 <- genus_freq1$Pragian[order(genus_freq1$Pragian,decreasing = TRUE)]
df10 <- genus_freq1$Emsian[order(genus_freq1$Emsian,decreasing = TRUE)]
df11 <- genus_freq1$Eifelian[order(genus_freq1$Eifelian,decreasing = TRUE)]
df12 <- genus_freq1$Givetian[order(genus_freq1$Givetian,decreasing = TRUE)]
df13 <- genus_freq1$Frasnian[order(genus_freq1$Frasnian,decreasing = TRUE)]
df14 <- genus_freq1$Famennian[order(genus_freq1$Famennian,decreasing = TRUE)]
#Carboniferous
df15 <- genus_freq1$Tournaisian[order(genus_freq1$Tournaisian,decreasing = TRUE)]


par(mfrow=c(2,4))
barplot(df8,ylab="Frequency", xlab="Number of genera", main="Lochkovian",col="#B3E1B6")
barplot(df9,ylab="Frequency", xlab="Number of genera", main="Pragian",col="#B3E1B6")
barplot(df10,ylab="Frequency", xlab="Number of genera", main="Emsian",col="#B3E1B6")
barplot(df11,ylab="Frequency", xlab="Number of genera", main="Eifelian",col="#B3E1B6")
barplot(df12,ylab="Frequency", xlab="Number of genera", main="Givetian",col="#B3E1B6")
barplot(df13,ylab="Frequency", xlab="Number of genera", main="Frasnian",col="#B3E1B6")
barplot(df14,ylab="Frequency", xlab="Number of genera", main="Famennian",col="#B3E1B6")
barplot(df15,ylab="Frequency", xlab="Number of genera", main="Tournaisian",col="#67A599")

#Carboniferous

df16 <- genus_freq1$Visean[order(genus_freq1$Visean,decreasing = TRUE)]
df17 <- genus_freq1$Serpukhovian[order(genus_freq1$Serpukhovian,decreasing = TRUE)]
df18 <- genus_freq1$Bashkirian[order(genus_freq1$Bashkirian,decreasing = TRUE)]
df19 <- genus_freq1$Moscovian[order(genus_freq1$Moscovian,decreasing = TRUE)]
df20 <- genus_freq1$Kasimovian[order(genus_freq1$Kasimovian,decreasing = TRUE)]
df21 <- genus_freq1$Gzhelian[order(genus_freq1$Gzhelian,decreasing = TRUE)]
#Permian
df22 <- genus_freq1$Asselian[order(genus_freq1$Asselian,decreasing = TRUE)]
df23 <- genus_freq1$Sakmarian[order(genus_freq1$Sakmarian,decreasing = TRUE)]


par(mfrow=c(2,4))
barplot(df16,ylab="Frequency", xlab="Number of genera", main="Visean",col="#67A599")
barplot(df17,ylab="Frequency", xlab="Number of genera", main="Serpukhovian",col="#67A599")
barplot(df18,ylab="Frequency", xlab="Number of genera", main="Bashkirian",col="#67A599")
barplot(df19,ylab="Frequency", xlab="Number of genera", main="Moscovian",col="#67A599")
barplot(df20,ylab="Frequency", xlab="Number of genera", main="Kasimovian",col="#67A599")
barplot(df21,ylab="Frequency", xlab="Number of genera", main="Gzhelian",col="#67A599")
barplot(df22,ylab="Frequency", xlab="Number of genera", main="Asselian",col="#F04028")
barplot(df23,ylab="Frequency", xlab="Number of genera", main="Sakmarian",col="#F04028")

#Permian

df24 <- genus_freq1$Artinskian[order(genus_freq1$Artinskian,decreasing = TRUE)]
df25 <- genus_freq1$Kungurian[order(genus_freq1$Kungurian,decreasing = TRUE)]
df26 <- genus_freq1$Roadian[order(genus_freq1$Roadian,decreasing = TRUE)]
df27 <- genus_freq1$Wordian[order(genus_freq1$Wordian,decreasing = TRUE)]
df28 <- genus_freq1$Capitanian[order(genus_freq1$Capitanian,decreasing = TRUE)]
df29 <- genus_freq1$Wuchiapingian[order(genus_freq1$Wuchiapingian,decreasing = TRUE)]
df30 <- genus_freq1$Changhsingian[order(genus_freq1$Changhsingian,decreasing = TRUE)]


par(mfrow=c(2,4))
barplot(df24,ylab="Frequency", xlab="Number of genera", main="Artinskian",col="#F04028")
barplot(df25,ylab="Frequency", xlab="Number of genera", main="Kungurian",col="#F04028")
barplot(df26,ylab="Frequency", xlab="Number of genera", main="Roadian",col="#F04028")
barplot(df27,ylab="Frequency", xlab="Number of genera", main="Wordian",col="#F04028")
barplot(df28,ylab="Frequency", xlab="Number of genera", main="Capitanian",col="#F04028")
barplot(df29,ylab="Frequency", xlab="Number of genera", main="Wuchiapingian",col="#F04028")
barplot(df30,ylab="Frequency", xlab="Number of genera", main="Changhsingian",col="#F04028")

