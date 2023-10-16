###########################################################################
#                                                                         #
#                        Alpha diversity analyses                         #  
#                                                                         #
###########################################################################
#                                                                         #
#                      Lisa Schnetz - September 2022                      #
#                                                                         #
###########################################################################


## Packages used in this script:

install.packages("magrittr")
install.packages("dplyr") 
install.packages("ggplot2") 
install.packages("deeptime") 
install.packages("divDyn")
install.packages("matrixStats")
install.packages("colorspace")

library(magrittr)
library(dplyr)
library(ggplot2)
library(deeptime)
library(divDyn)
library(colorspace)

## Read in the datasets:

Palshark_data <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")
acantho_data <- read.csv2("./data/Acanthodian_input_R_sampling.csv")

shark_data <- rbind(Palshark_data, acantho_data)

## Remove any occurrences that have no locality information:

shark_data <- filter(shark_data, COLLECTION != "")


#######################################################################
#             Alpha diversity analysis total chondrichthyans         #
#######################################################################


## The script used from here on is mainly taken from local_richness depository from Dr. Emma Dunne. https://github.com/emmadunne/local_richness.git
## The method is based on the paper by Close et al. (2019).

## First count the number of taxa per locality using the table() function

freq_table <- as.data.frame(table(shark_data[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
head(freq_table) 

## Then we select some additional information from the original dataset
locality_info <- select(shark_data, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()


## Lastly, we'll add it to the frequency table  before plotting
alpha_data <- left_join(freq_table, locality_info, by = "COLLECTION")



## Plot the data
mytheme_alpha <- theme(panel.background = element_blank(),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       panel.grid.major.y = element_blank(),
                       panel.grid.major.x = element_blank(),
                       legend.position="right",
                       panel.border = element_rect(colour = "black", fill = NA),
                       axis.text.x = element_text(size=10, angle=0, hjust=0.5),
                       axis.text.y = element_text(size=10),
                       axis.title = element_text(size=12))

alpha_plot <- ggplot() +
  
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
  
  geom_point(data = alpha_data, aes(MID.POINT, Freq), colour = adjustcolor("turquoise4", alpha.f = 0.4), size = 6) + 
  mytheme_alpha +
  labs(x = "Time (Ma)", y = "Number of genera") +
  scale_x_reverse(expand = c(0, 0), limits = c(470, 252), breaks = c(251.902, 298.9, 358.9, 419.2,443.8,467.3)) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25,50,75,100), expand = c(0, 0))


#add time scale to your plot
alpha_plot <- gggeo_scale(alpha_plot, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
alpha_plot <- gggeo_scale(alpha_plot , dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = alpha_plot,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/alpha_plot.pdf", useDingbats=FALSE)



#####################################################################
#             Alpha diversity analysis by environment               #
#####################################################################

## Subset the data to only contain occurrences from freshwater versus marine environments

shark_datamarine <- subset(shark_data, BA_1=="1"|BA_2=="1"|BA_3=="1"|BA_4=="1"|BA_5=="1"|BA_6=="1" )
shark_data0 <- subset(shark_data, BA_0=="1")

## Create frequency table for each BA, add additional information and combine data
freq_table <- as.data.frame(table(shark_datamarine[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION"

locality_info <- select(shark_datamarine, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_datamarine <- left_join(freq_table, locality_info, by = "COLLECTION") 


freq_table <- as.data.frame(table(shark_data0[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION"
locality_info <- select(shark_data0, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data0 <- left_join(freq_table, locality_info, by = "COLLECTION") # warning is ok, R is just cranky again ;)


## Plot the data

alpha_plotenv <- ggplot() +
  
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
  
  geom_point(data = alpha_datamarine, aes(MID.POINT, Freq), colour = adjustcolor("blue", alpha.f = 0.4), size = 6) + 
  geom_point(data = alpha_data0, aes(MID.POINT, Freq), colour = adjustcolor("red", alpha.f = 0.4), size = 6)+
  mytheme_alpha +
  labs(x = "Time (Ma)", y = "Number of genera") +
  scale_x_reverse(expand = c(0, 0), limits = c(470, 250), breaks = c(251.902, 298.9, 358.9, 419.2,443.8,467.3)) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25,50,75,100), expand = c(0, 0))

#add time scale to your plot
alpha_plotenv <- gggeo_scale(alpha_plotenv, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
alpha_plotenv <- gggeo_scale(alpha_plotenv, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = alpha_plotenv,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/alpha_plot_environment.pdf", useDingbats=FALSE)




#####################################################################
#                   Alpha diversity analysis by BAs                #
#####################################################################

## Subset the data to only contain occurrences from each BA

shark_data0 <- subset(shark_data, BA_0=="1")
shark_data1 <- subset(shark_data, BA_1=="1" )
shark_data2 <- subset(shark_data, BA_2=="1" )
shark_data3 <- subset(shark_data, BA_3=="1" )
shark_data4 <- subset(shark_data, BA_4=="1" )
shark_data5 <- subset(shark_data, BA_5=="1" )
shark_data6 <- subset(shark_data, BA_6=="1" )


## Create frequency table for each BA, add additional information and combine data
# BA 0
freq_table <- as.data.frame(table(shark_data0[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION"
locality_info <- select(shark_data0, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data0 <- left_join(freq_table, locality_info, by = "COLLECTION")

# BA 1
freq_table <- as.data.frame(table(shark_data1[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data1, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data1 <- left_join(freq_table, locality_info, by = "COLLECTION") 

# BA 2
freq_table <- as.data.frame(table(shark_data2[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data2, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data2 <- left_join(freq_table, locality_info, by = "COLLECTION") 

# BA 3
freq_table <- as.data.frame(table(shark_data3[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data3, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data3 <- left_join(freq_table, locality_info, by = "COLLECTION") 

# BA 4
freq_table <- as.data.frame(table(shark_data4[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data4, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data4 <- left_join(freq_table, locality_info, by = "COLLECTION") 

# BA 5
freq_table <- as.data.frame(table(shark_data5[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data5, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data5 <- left_join(freq_table, locality_info, by = "COLLECTION") 

# BA 6
freq_table <- as.data.frame(table(shark_data6[, c('COLLECTION','GENUS')]$COLLECTION))
names(freq_table)[1] <- "COLLECTION" 
locality_info <- select(shark_data6, COLLECTION,  
                        EARLIEST, MID.POINT) %>% distinct(COLLECTION, .keep_all = TRUE) %>% na.omit()
alpha_data6 <- left_join(freq_table, locality_info, by = "COLLECTION") 

## Plot the data

## BA0 - BA2

alpha_plotBA0to2 <- ggplot() +
  
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
  
  geom_point(data = alpha_data0, aes(MID.POINT, Freq), colour = adjustcolor("#fc8d59", alpha.f = 0.4), show.legend = T,size = 6)+
  geom_point(data = alpha_data1, aes(MID.POINT, Freq), colour = adjustcolor("#ffffbf", alpha.f = 0.4),show.legend = T, size = 6) + 
  geom_point(data = alpha_data2, aes(MID.POINT, Freq), colour = adjustcolor("#91bfdb", alpha.f = 0.4), show.legend = T,size = 6)+
  mytheme_alpha +
  labs(x = "Time (Ma)", y = "Number of genera") +
  scale_x_reverse(expand = c(0, 0), limits = c(470, 250), breaks = c(251.902, 298.9, 358.9, 419.2,443.8,467.3)) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25,50,75,100), expand = c(0, 0))

#add time scale to your plot
alpha_plotBA0to2 <- gggeo_scale(alpha_plotBA0to2, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
alpha_plotBA0to2 <- gggeo_scale(alpha_plotBA0to2, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = alpha_plotBA0to2,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/alpha_plot_BA0BA2.pdf", useDingbats=FALSE)

## BA3 - BA6

alpha_plotBA3to6 <- ggplot() +
  
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
  
  geom_point(data = alpha_data3, aes(MID.POINT, Freq), colour = adjustcolor("#ca0020", alpha.f = 0.4), size = 6)+
  geom_point(data = alpha_data4, aes(MID.POINT, Freq), colour = adjustcolor("green", alpha.f = 0.4), size = 6) + 
  geom_point(data = alpha_data5, aes(MID.POINT, Freq), colour = adjustcolor("purple", alpha.f = 0.4), size = 6)+
  geom_point(data = alpha_data6, aes(MID.POINT, Freq), colour = adjustcolor("#0571b0", alpha.f = 0.4), size = 6)+
  mytheme_alpha +
  labs(x = "Time (Ma)", y = "Number of genera") +
  scale_x_reverse(expand = c(0, 0), limits = c(470, 250), breaks = c(251.902, 298.9, 358.9, 419.2,443.8,467.3)) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25,50,75,100), expand = c(0, 0))

#add time scale to your plot
alpha_plotBA3to6 <- gggeo_scale(alpha_plotBA3to6, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
alpha_plotBA3to6 <- gggeo_scale(alpha_plotBA3to6, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = alpha_plotBA3to6,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/alpha_plot_BA3BA6.pdf", useDingbats=FALSE)

## all BAs

alpha_plotallBAs <- ggplot() +
  
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
  
  geom_point(data = alpha_data0, aes(MID.POINT, Freq), colour = adjustcolor("orange", alpha.f = 0.4), show.legend = T,size = 6)+
  geom_point(data = alpha_data1, aes(MID.POINT, Freq), colour = adjustcolor("yellow", alpha.f = 0.4),show.legend = T, size = 6) + 
  geom_point(data = alpha_data2, aes(MID.POINT, Freq), colour = adjustcolor("cyan", alpha.f = 0.4), show.legend = T,size = 6)+
  geom_point(data = alpha_data3, aes(MID.POINT, Freq), colour = adjustcolor("#ca0020", alpha.f = 0.4), size = 6)+
  geom_point(data = alpha_data4, aes(MID.POINT, Freq), colour = adjustcolor("green", alpha.f = 0.4), size = 6) + 
  geom_point(data = alpha_data5, aes(MID.POINT, Freq), colour = adjustcolor("purple", alpha.f = 0.4), size = 6)+
  geom_point(data = alpha_data6, aes(MID.POINT, Freq), colour = adjustcolor("#0571b0", alpha.f = 0.4), size = 6)+
  mytheme_alpha +
  labs(x = "Time (Ma)", y = "Number of genera") +
  scale_x_reverse(expand = c(0, 0), limits = c(470, 250), breaks = c(251.902, 298.9, 358.9, 419.2,443.8,467.3)) +
  scale_y_continuous(limits = c(0, 110), breaks = c(0, 25,50,75,100), expand = c(0, 0))


#add time scale to your plot
alpha_plotallBAs <- gggeo_scale(alpha_plotallBAs, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
alpha_plotallBAs <- gggeo_scale(alpha_plotallBAs, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

ggsave(plot = alpha_plotallBAs,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/alpha_plot_allBAs.pdf", useDingbats=FALSE)

