setwd("~/Google Drive/proj/A13-YuanDai/Data")
rawdata=read.csv('Iowa Priavte Well As-SHL data to 2020.csv',na.strings = "")
library(dplyr)
library(ggplot2)
library(tidyverse)
library(CompRandFld)
library(RandomFields)
library(plotly)
library(maps)
library(sf) 
library(leaflet)
data=rawdata%>%filter(!is.na(Latitude) & !is.na(Longitude))%>%select(Collected.Date,Location.City,Sample.Description,Arsenic,MCL,Latitude,Longitude)%>%
    mutate(Arsenic2=as.numeric(replace(as.character(Arsenic), as.character(Arsenic)=='<0.001', '0.0005')))%>%mutate(Y=as.numeric(Arsenic2>=0.01))%>%filter(!is.na(Y))
summary(data)
 


## Get the states and county map, turn into sf object
US_state <- st_as_sf(maps::map("state", plot = FALSE, fill = TRUE))
## Test the function using points in Wisconsin and Oregon
testPoints <- data%>%select(Latitude,Longitude)
# Make it a spatial dataframe, using the same coordinate system as the US spatial dataframe
testPoints <- st_as_sf(testPoints, coords = c("Longitude", "Latitude"), crs = st_crs(US_state))
data$State=st_join(testPoints, US_state)$ID
US_county <- st_as_sf(maps::map("county", plot = FALSE, fill = TRUE)) 
data$County=st_join(testPoints, US_county)$ID
 
data=data%>%filter(State=='iowa')%>%separate(County,c('State Name','County'),sep=',')
  
  

IW_latlonzoom=c(-93.22,42.2280,7)

 
library(ggmap)
library(ggplot2)
register_google("AIzaSyB5l0SAlK3a4VOcvzg_a8-qT7xcHYeCe-s") ## this is my key. You need to register to get your own google API  
 
iowamap=get_map('iowa',zoom=7,maptype='roadmap')
iowadata=data%>%select(Latitude,Longitude,Y)

ggmap(iowamap)+geom_point(data=iowadata,aes(x=Longitude,y=Latitude,col=as.factor(Y),shape=as.factor(Y)))+
  labs(title="Water Quality (Arsenic)", x ="longitude", y = "latitude")+
  scale_colour_discrete(name  ="Arsenic",breaks=c(0,1),labels=c("< MCL", "Exceed MCL")) +
  scale_shape_discrete(name  ="Arsenic", breaks=c(0,1), labels=c("< MCL", "Exceed MCL"))+theme(legend.text = element_text(size=12))

write.csv(data,file='PreProcessedArsenic.Rdata')

#dim(data%>%filter(!is.na(Sample.Description)))%>%mutate(Address=replace(Sample,)
#data[c(which(str_detect(data$Sample.Description,c(' st'))), which(str_detect(data$Sample.Description,c(' ct'))),
#       which(str_detect(data$Sample.Description,c(' blvd'))),which(str_detect(data$Sample.Description,c(' rd')))
#       ) ,]

#library(R.matlab)
#writeMat()
# 
# leaflet(data) %>%
#   addProviderTiles("CartoDB") %>%
#   addMapPane("blocks", 410) %>% 
#   addMapPane("districts", 420) %>% 
#   setView(IW_latlonzoom[1],IW_latlonzoom[2], IW_latlonzoom[3])%>%
#   addMarkers(~Longitude, ~Latitude, popup = ~as.character(Arsenic2), label = ~as.character(Arsenic2))
# 


 
