setwd("C:/Users/clee41/OneDrive - UBC/Desktop/GradWork/Data/2023/growth_curves/Jan 9 2023-sorbitol glucose mannitol")

library(ggplot2)
library(tidyverse)
library(readr)

data = read.csv("growthdata.csv")

sumdata = data %>%
    group_by(strain, timepoint, treatment, concentration)  %>%
    summarize(meanOD = mean(OD600),sd = sd(OD600))

medias = unique(data$treatment)

for (m in medias){
  media_data = sumdata %>%
    filter(treatment ==m)
  print(m)
  conc = unique(media_data$concentration)
  for (c in conc){
    
    conc_data = media_data %>%
      filter(concentration==c)
    print(c)
    plot1 = ggplot(conc_data, aes(timepoint,meanOD,color = strain))+
        geom_line(size = 2, linetype = "longdash")+
        geom_errorbar(aes(ymin = meanOD-sd,ymax = meanOD+sd),size =1.5)+
        scale_color_manual(values = c("#FF808C","#FFC285","#8080FC"))+
        theme_bw()+
        labs(x="Time (hours)", y = "Average OD600")+
        theme(axis.text = element_text(size = 15, face = "bold"),
              axis.title = element_text(size = 15, face = "bold"),
              legend.text = element_text(size = 15, face = "bold"),
              legend.title =element_text(size = 17, face = "bold") )+
        ylim(0,0.9)
    ggsave(paste(c,m,"curve.png",sep = "_"),plot = plot1) 
  }
  
  
}



