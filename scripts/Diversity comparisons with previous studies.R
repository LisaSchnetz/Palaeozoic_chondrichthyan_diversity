
############
install.packages("tidyverse")
install.packages("deeptime")
install.packages("ggplot2")

library(tidyverse)
#library(geoscale)
library(deeptime)
library(ggplot2)
###########

#######Feichtinger et al. 2021 data###############

feichtinger <- data.frame(SIB  = c(29,37,40,30, 32,23,21), EMSD = c(20.5,27.5, 28, 25.5, 24, 19.5, 16),
                          max_age = c(358.9,346.7,330.9,323.2,315.2,307,303.7),
                          min_age = c(346.7,330.9,323.2,315.2,307,303.7,298.9),
                          mid_age = c(352.80, 338.80, 327.05, 319.20, 311.10, 305.35, 301.30))


###########Sallan & Coates 2010 data#########

sallan <- data.frame(acan_div = c(29, 21, 11, 8, 7, 6),
           chond_div = c(11, 20, 36, 66, 94, 77),
           total = c(40, 41, 47, 74,101,83),
max_age = c(387.7,382.7,372.2,358.9,346.7,330.9),
min_age = c(382.7,372.2,358.9,346.7,330.9,323.2),
mid_age = c(385.20, 377.45, 365.55, 352.80, 338.80, 327.05))


######################################

###################################################
#######################Plot curves##################
##################################################

###Feichtinger et al. 2021

curves <- ggplot() +
   geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  
  geom_line(data=feichtinger, aes(x = mid_age, y = SIB,colour="cyan4"), size=1) +
  geom_line(data=feichtinger, aes(x = mid_age, y = EMSD,colour="cyan4"), size=1, linetype=2) +
  
  scale_colour_manual(values ="cyan4") +
  scale_x_reverse(expand=c(0,0), limits = c(393.3, 293.52), breaks = c(250,300,350)) +
  labs(x = "Time (Ma)", y = "Diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25,50), limits = c(0, 50)) +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))


curves <- gggeo_scale(curves, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
curves <- gggeo_scale(curves, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

windows(30,20)  
curves

ggsave(plot = curves,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/curve_feichtinger.pdf", useDingbats=FALSE)

###########Sallan & Coates 2010#########

curves2 <- ggplot() +
  geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = 0, ymax = Inf),fill= "grey90",linetype="blank") +
  
  geom_line(data=sallan, aes(x = mid_age, y = chond_div,colour= "Chondrichthyes"), size=1) +
  geom_line(data=sallan, aes(x = mid_age, y = acan_div,colour="Acanthodii"), size=1) +
  geom_line(data=sallan, aes(x = mid_age, y = total,colour="Total"), size=1) +
  
  scale_colour_manual(values =c("Chondrichthyes" ="blue", "Acanthodii" = "red", "Total" ="black")) +
  scale_x_reverse(expand=c(0,0), limits = c(393.3, 293.52), breaks = c(250,300,350)) +
  labs(x = "Time (Ma)", y = "Diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,25,50,75,100), limits = c(0, 110)) +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))


curves2 <- gggeo_scale(curves2, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
curves2 <- gggeo_scale(curves2, dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

windows(30,20)  
curves2

ggsave(plot = curves2,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/curve_sallan.pdf", useDingbats=FALSE)