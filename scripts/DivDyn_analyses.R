######Diversity analyses Palaeozoic chondrichthyans####

install.packages("divDyn")
install.packages("matrixStats")
install.packages("colorspace")

library(matrixStats)
library(divDyn)
library(colorspace)


###Dataset with sp. species!####

sharks <- read.csv2("./data/Chondrichthyes_input_R_sampling.csv")

Acanthodians <-read.csv2("./data/Acanthodian_input_R_sampling.csv")

allsharks <- rbind(sharks, Acanthodians)



allsharks["new_bin"] <- NA
allsharks$new_bin <- as.numeric(as.character(allsharks$new_bin))


data(stages)
data(tens)

bins <- c(467.3,458.4,453.0,445.2, 443.8,440.8,438.5,433.4,430.5,427.4,
          425.6,423.0,419.2,410.8,407.6,393.3,387.7,382.7,372.2,358.9,346.7,330.9,323.2,315.2,307.0,303.7,
          298.9,293.52,290.1,283.5,272.95,268.8,265.1,259.1,254.14,251.902)

binDframe <- data.frame(bin = c("Darriwilian", 
                                "Sandbian",
                                "Katian",
                                "Hirnantian",
                                "Rhuddanian",
                                "Aeronian",
                                "Telychian",
                                "Sheinwoodian",
                                "Homerian",
                                "Gorstian",
                                "Ludfordian",
                                "Pridoli",
                                "Lochkovian",
                                "Pragian",
                                "Emsian",
                                "Eifelian",
                                "Givetian",
                                "Frasnian",
                                "Famennian",
                                "Tournaisian",
                                "Visean",
                                "Serpukhovian",
                                "Bashkirian",
                                "Moscovian",
                                "Kasimovian",
                                "Gzhelian",
                                "Asselian",
                                "Sakmarian",
                                "Artinskian",
                                "Kungurian",
                                "Roadian",
                                 "Wordian",
                                "Capitanian",
                                "Wuchiapingian",
                                "Changhsingian"), 
                        FAD = as.numeric(bins[1:(length(bins)-1)]), 
                        LAD = as.numeric(bins[2:(length(bins))]))

binDframe$mid <- (binDframe$FAD+binDframe$LAD)/2




#=== Bin by midpoint ===
#for(s in 1:nrow(binDframe)){
  #for(o in 1:nrow(allsharks)){
   # if(allsharks$MID.POINT[o] <= binDframe$FAD[s] && allsharks$MID.POINT[o] >= binDframe$LAD[s]){
     # allsharks$new_bin[o] <- stages$stg[s]}  }}
#df1 <- allsharks




#=== Bin in any bin it falls within ===
All_Bin_List <- list()
for(s in 1:nrow(binDframe)){ # for each new bin
  temp_recs <- data.frame()
  for (o in 1:nrow(allsharks)){ 
    if (allsharks$MAX_DATE[o] > binDframe$LAD[s] && allsharks$MIN_DATE[o] < binDframe$FAD[s]){ # If occurrence max. age is greater than Bin min. age AND if occurrence min. age is less then Bin max. age (i.e. falls within bin at some point)
      temp_recs <- rbind(temp_recs, allsharks[o,]) # Add that occurrence to binlist
    }
  }
  if (nrow(temp_recs) > 0){
    temp_recs$new_bin <- s
  }
  All_Bin_List[[s]] <- temp_recs
}
df2 <- do.call("rbind", All_Bin_List)

bin_info <- binstat(df2, tax="GENUS", bin="new_bin")
bin_info[bin_info==0] <- NA

ddsharks<-divDyn(df2, tax="GENUS", bin="new_bin")

####ggplot instead of tplot for Second-for-third origination/extinction rates####

ggplotorigination <-cbind(ddsharks,binDframe)

write.csv(ggplotorigination,"D:/Third chapter data 051021/divdyn_results.csv", row.names = FALSE)

#theme_iNEXT <- theme(panel.background = element_blank(),
                    # panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
                    # panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(),
                    # panel.border = element_rect(colour = "black", fill = NA),
                    # axis.text.x = element_text(size=14, angle=0, hjust=0.5),
                     #axis.text.y = element_text(size=14),
                     #axis.title = element_text(size=12))

origination <- ggplot(ggplotorigination, aes(x=mid)) +
  geom_rect(aes(xmax=467.3, xmin = 458.4, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=453, xmin = 445.2, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=443.8, xmin = 440.8, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=438.5, xmin = 433.4, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=430.5, xmin = 427.4, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=425.6, xmin = 423, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=419.2, xmin = 410.8, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=290.1, xmin = 283.5, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=273.01, xmin = 266.9, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=264.28, xmin = 259.51, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  geom_rect(aes(xmax=254.14, xmin = 251.902, ymin = -0.15, ymax = Inf),fill= "grey90",linetype="blank") +
  theme_classic()+
  theme(legend.position="top",legend.background = element_rect(size=0.5))+
  theme(plot.margin=margin(0.5,0.75,0.5,0.5,"cm"))+
  geom_line(data=ggplotorigination, aes(x = mid, y = ori2f3, colour = "blue"), size = 1) +
  geom_line(data=ggplotorigination, aes(x = mid, y = ext2f3, colour="red"),size = 1)+
  scale_colour_manual(values =c("blue","red"), labels=c("origination","extinction")) +
  scale_x_reverse(expand=c(0,0), limits = c(467.3, 251.902), breaks = c(250,300,350,400,450)) +
  scale_y_continuous(expand=c(0,0), breaks = c(0,0.5, 1.0, 1.5, 2.0), limits = c(-0.15, 2.0))+
  labs(x = "Ma", y = "Rate")
  #geom_vline(aes(xintercept = 358.9), linetype = "longdash", colour = "red", size=1)  #Hangenberg extinction

origination_scale <- gggeo_scale(origination, dat = "periods", height = unit(1.5, "lines"),  size = 4, abbrv = FALSE)
origination_scale <- gggeo_scale(origination_scale , dat = "stages", height = unit(1.5, "lines"),  size = 3, abbrv = TRUE)

windows(30,20)  
origination_scale

ggsave(plot = origination_scale,
       width = 20, height = 15, dpi = 600, units = "cm", 
       filename = "./plots/origination_plot.pdf", useDingbats=FALSE)




squaresplot2 <- ggplot(ggplotorigination, aes(x=mid)) + theme_iNEXT+
  geom_rect(aes(xmax=467.3, xmin = 458.4, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=453, xmin = 445.2, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=443.8, xmin = 440.8, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=438.5, xmin = 433.4, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=430.5, xmin = 427.4, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=425.6, xmin = 423, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=419.2, xmin = 410.8, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=407.6, xmin = 393.3, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=387.7, xmin = 382.7, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=372.2, xmin = 358.9, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=346.7, xmin = 330.9, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=323.2, xmin = 315.2, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=307, xmin = 303.7, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=298.9, xmin = 293.52, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=290.1, xmin = 283.5, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=273.01, xmin = 266.9, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=264.28, xmin = 259.51, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_rect(aes(xmax=254.14, xmin = 251.902, ymin = -0.15, ymax = Inf),fill= "grey90") +
  geom_line(data=ggplotorigination, aes(x = mid, y = ori2f3, colour = "blue"), size = 2) +
  geom_line(data=ggplotorigination, aes(x = mid, y = ext2f3, colour="red"),size = 2)+
  scale_colour_manual(name="Rate",values =c("blue","red"), labels=c("origination","extinction")) +
  scale_x_reverse(limits = c(470, 250)) + labs(x = "Ma", y = "Squares diversity") +
  scale_y_continuous(expand=c(0,0), breaks = c(0,0.5, 1.0, 1.5, 2.0), limits = c(-0.15, 2.0))+
  geom_vline(aes(xintercept = 358.9), linetype = "longdash", colour = "red", size=1)  #Hangenberg extinction

squaresplot2

scale <- gggeo_scale(squaresplot2, dat = "stages",height = unit(4, "lines"), rot = 90, size = 2.5, abbrv = FALSE)
scale
