######Diversity analyses Palaeozoic chondrichthyans####

install.packages("divDyn")
install.packages("matrixStats")
install.packages("colorspace")

library(matrixStats)
library(divDyn)
library(colorspace)

setwd("C:/Users/Install/Desktop/Lisa/Desktop/Third chapter")

###Dataset with sp. species!####

Acanthodians <-read.csv2("D:/Third chapter data 051021/Acanthodian_completeness_input_R.csv")
Acanthodians <-subset(Acanthodians, DUPLICATE=="N")

sharks <-read.csv("D:/Third chapter data 051021/Chondrichthyes_input_R.csv")
sharks <-subset(sharks, DUPLICATE=="N")

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
     # allsharks$new_bin[o] <- stages$stg[s]
    }
  }
}
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


SQS1 <- subsample(df2,iter=100, q=0.4, tax="GENUS", bin="new_bin", 
                 output="dist", type="sqs", 
                 duplicates = TRUE, useFailed = TRUE)

SQS2 <- subsample(df2,iter=100, q=0.6, tax="GENUS", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)

SQS3 <- subsample(df2,iter=100, q=0.8, tax="GENUS", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)


tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))
shades(binDframe$mid, SQS1$divRT, col="red")
shades(binDframe$mid, SQS2$divRT, col="blue")
shades(binDframe$mid, SQS3$divRT, col="green")

lines(binDframe$mid, bin_info$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("raw", "SQS (0.4)", "SQS (0.6)", "SQS (0.8)"),
       col=c("black","red", "blue", "green"), lty=c(1,1,1,1), lwd=c(2,2,2,2))


ddsharks<-divDyn(df2, tax="GENUS", bin="new_bin")


# origination plus extinction rate plot
tsplot(stages, boxes=c("short", "system"), shading="stage", xlim=16:52,ylim=c(0,60),
       ylab="Number of taxa originating/ getting extinct",boxes.col=c("col", "systemCol"))
lines(binDframe$mid, ddsharks$tOri, lwd=2, lty=1, col="blue")
lines(binDframe$mid, ddsharks$tExt, lwd=2, lty=1, col="red")

legend("topleft", inset=c(0.1,0.1), legend=c("origination", "extinction"),
       lwd=2, lty=c(1,1), col=c("blue", "red"), bg="white")



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

####Species level#########


bin_info_spec <- binstat(df2, tax="SPECIES", bin="new_bin")
bin_info_spec[bin_info_spec==0] <- NA

SQSS1 <- subsample(df2,iter=100, q=0.4, tax="SPECIES", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)

SQSS2 <- subsample(df2,iter=100, q=0.6, tax="SPECIES", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)

SQSS3 <- subsample(df2,iter=100, q=0.8, tax="SPECIES", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)


tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (species)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info_spec$SIBs, na.rm = TRUE)+(max(bin_info_spec$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))
shades(binDframe$mid, SQSS1$divSIB, col="red")
shades(binDframe$mid, SQSS2$divSIB, col="blue")
shades(binDframe$mid, SQSS3$divSIB, col="green")

lines(binDframe$mid, bin_info_spec$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("SQS (0.4)", "SQS (0.6)", "SQS (0.8)"),
       col=c("red", "blue", "green"), lty=c(1,1,1), lwd=c(2,2,2))


####Seperate major chondrichthyan groups and look at their diversity compared to each other######

Cladodontomorphi <- subset(df2, ORDER=="Symmoriiformes"| ORDER=="Cladoselachiiformes"| ORDER=="Ctenacanthiformes"| ORDER=="Squatinactiformes")

bin_info_clad <- binstat(Cladodontomorphi, tax="GENUS", bin="new_bin")
bin_info_clad[bin_info_clad==0] <- NA


add_row1 <- c(34, NA,NA,NA,NA,NA,NA)
add_row2 <- c(35, NA,NA,NA,NA,NA,NA) 
bin_info_clad <- rbind(bin_info_clad,add_row1, add_row2)

Xenacanthimorpha <- subset(df2, ORDER=="Bransonelliformes"| ORDER=="Xenacanthiformes")

bin_info_xena <- binstat(Xenacanthimorpha, tax="GENUS", bin="new_bin")
bin_info_xena[bin_info_xena==0] <- NA

#add_row1 <- c(34, NA,NA,NA,NA,NA,NA)
#add_row2 <- c(35, NA,NA,NA,NA,NA,NA) 
#bin_info_xena <- rbind(bin_info_xena,add_row1, add_row2)

Euselachii <- subset(df2, ORDER=="Hybodontiformes"| ORDER=="Synechodontiformes")

bin_info_eusel <- binstat(Euselachii, tax="GENUS", bin="new_bin")
bin_info_eusel[bin_info_eusel==0] <- NA

#add_row1 <- c(35, NA,NA,NA,NA,NA,NA) 
#bin_info_eusel <- rbind(bin_info_eusel,add_row1)

Euchondrocephali <- subset(df2, ORDER=="Orodontiformes"| ORDER=="Eugeneodontiformes"| ORDER=="Petalodontiformes"
                           | ORDER=="Chimaeriformes"| ORDER=="Chondrenchelyiformes"| ORDER=="Cochliodontiformes"
                           | ORDER=="Copodontiformes"| ORDER=="Helodontiformes"| ORDER=="Iniopterygia"| ORDER=="Menaspiformes"| ORDER=="Psammodontiformes")

bin_info_euchondro <- binstat(Euchondrocephali, tax="GENUS", bin="new_bin")
bin_info_euchondro[bin_info_euchondro==0] <- NA


tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))

lines(binDframe$mid, bin_info_clad$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)
lines(binDframe$mid, bin_info_xena$SIBs, type = "o", pch = 21, col = "orange", bg = "orange", lwd = 2, cex = 1)
lines(binDframe$mid, bin_info_eusel$SIBs, type = "o", pch = 21, col = "red", bg = "red", lwd = 2, cex = 1)
lines(binDframe$mid, bin_info_euchondro$SIBs, type = "o", pch = 21, col = "blue", bg = "blue", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("Cladodontomorphi", "Xenacanthimorpha", "Euselachii","Euchondrocephali"),
       col=c("black", "orange", "red","blue"), lty=c(1,1,1,1), lwd=c(2,2,2,2))


###Seperate geographical regions and compare!!!!!!

North <- subset(df2, HEMISPHERE=="N"|HEMISPHERE=="N/S")

bin_info_north <- binstat(North, tax="GENUS", bin="new_bin")
bin_info_north[bin_info_north==0] <- NA


South <- subset(df2, HEMISPHERE=="S"|HEMISPHERE=="N/S") 

bin_info_south <- binstat(South, tax="GENUS", bin="new_bin")
bin_info_south[bin_info_south==0] <- NA

add_row1 <- c(35, NA,NA,NA,NA,NA,NA) 
bin_info_south <- rbind(bin_info_south,add_row1)



tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))

lines(binDframe$mid, bin_info_north$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)
lines(binDframe$mid, bin_info_south$SIBs, type = "o", pch = 21, col = "orange", bg = "orange", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("North", "South"),
       col=c("black", "orange"), lty=c(1,1), lwd=c(2,2))


###Compare freshwater versus marine diversity#####

fresh <-subset(df2, BA_0=="1")

bin_info_fresh <- binstat(fresh, tax="GENUS", bin="new_bin")
bin_info_fresh[bin_info_fresh==0] <- NA


SQSfre1 <- subsample(fresh,iter=100, q=0.4, tax="GENUS", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)

SQSfre2 <- subsample(fresh,iter=100, q=0.6, tax="GENUS", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)

SQSfre3 <- subsample(fresh,iter=100, q=0.8, tax="GENUS", bin="new_bin", 
                  output="dist", type="sqs", 
                  duplicates = TRUE, useFailed = TRUE)


marine <-subset(df2, BA_1=="1"| BA_2=="1"| BA_3=="1"|BA_4=="1"|BA_5=="1"|BA_6=="1")

bin_info_marine <- binstat(marine, tax="GENUS", bin="new_bin")
bin_info_marine[bin_info_marine==0] <- NA

SQSmar1 <- subsample(marine,iter=100, q=0.4, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 
                     duplicates = TRUE, useFailed = TRUE)

SQSmar2 <- subsample(marine,iter=100, q=0.6, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 
                     duplicates = TRUE, useFailed = TRUE)

SQSmar3 <- subsample(marine,iter=100, q=0.8, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 
                     duplicates = TRUE, useFailed = TRUE)



 
tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info$SIBs, na.rm = TRUE)+(max(bin_info$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))

lines(binDframe$mid, bin_info_fresh$SIBs, type = "o", pch = 21, col = "green", bg = "green", lwd = 2, cex = 1)
lines(binDframe$mid, bin_info_marine$SIBs, type = "o", pch = 21, col = "blue", bg = "blue", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("freshwater", "marine"),
       col=c("green", "blue"), lty=c(1,1), lwd=c(2,2))

##Plot freshwater with SQS####
tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info_fresh$SIBs, na.rm = TRUE)+(max(bin_info_fresh$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))
shades(binDframe$mid, SQSfre1$divSIB, col="greenyellow")
shades(binDframe$mid, SQSfre2$divSIB, col="green")
shades(binDframe$mid, SQSfre3$divSIB, col="green4")

lines(binDframe$mid, bin_info_fresh$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("raw", "SQS (0.4)", "SQS (0.6)", "SQS (0.8)"),
      col=c("black", "greenyellow", "green", "green4"), lty=c(1,1,1,1), lwd=c(2,2,2,2), pch=c(19,NA,NA,NA))

###Plot marine with SQS####
tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info_marine$SIBs, na.rm = TRUE)+(max(bin_info_marine$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))

shades(binDframe$mid, SQSmar1$divSIB, col="cyan")
shades(binDframe$mid, SQSmar2$divSIB, col="skyblue")
shades(binDframe$mid, SQSmar3$divSIB, col="royalblue4")

lines(binDframe$mid, bin_info_marine$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("raw", "SQS (0.4)", "SQS (0.6)", "SQS (0.8)"),
       col=c("black", "cyan", "skyblue", "royalblue4"), lty=c(1,1,1,1), lwd=c(2,2,2,2),pch=c(19,NA,NA,NA))


####Plot anadromous diversity####

anad <-subset(df2, BA_0=="1"| BA_1=="1")

bin_info_anad <- binstat(anad, tax="GENUS", bin="new_bin")
bin_info_anad[bin_info_anad==0] <- NA


SQSanad1 <- subsample(anad,iter=100, q=0.4, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 
                     duplicates = TRUE, useFailed = TRUE)

SQSanad2 <- subsample(anad,iter=100, q=0.6, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 
                     duplicates = TRUE, useFailed = TRUE)

SQSanad3 <- subsample(anad,iter=100, q=0.8, tax="GENUS", bin="new_bin", 
                     output="dist", type="sqs", 

                                          duplicates = TRUE, useFailed = TRUE)

tsplot(stages, boxes=c("short", "system"), shading="stage",ylab = "Sampled-in-bin diversity (genera)",
       boxes.col=c("col", "systemCol"), ylim=c(0,(max(bin_info_anad$SIBs, na.rm = TRUE)+(max(bin_info_anad$SIBs, na.rm = TRUE)*0.1))),xlim=c(16:52),plot.args=list(axes=TRUE, main="Palaeozoic"))

shades(binDframe$mid, SQSanad1$divSIB, col="yellow")
shades(binDframe$mid, SQSanad2$divSIB, col="orange")
shades(binDframe$mid, SQSanad3$divSIB, col="red")

lines(binDframe$mid, bin_info_anad$SIBs, type = "o", pch = 21, col = "black", bg = "black", lwd = 2, cex = 1)

legend("topleft", bg="white", legend=c("raw", "SQS (0.4)", "SQS (0.6)", "SQS (0.8)"),
       col=c("black", "yellow", "orange", "red"), lty=c(1,1,1,1), lwd=c(2,2,2,2),pch=c(19,NA,NA,NA))



############Sampling biases#############

sharks_sampling <-read.csv2("C:/Users/Install/Desktop/Lisa/Desktop/Third chapter/Chondrichthyes_input_R_sampling.csv")


