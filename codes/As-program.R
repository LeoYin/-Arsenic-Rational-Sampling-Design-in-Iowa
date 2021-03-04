rawdata=read.csv('/Users/yinlihao/Desktop/projects/Project-As/Iowa Priavte Well As-SHL data to 2020.csv'
                 ,na.strings = '')

library(dplyr)
library(ggplot2)
library(tidyverse)
library(CompRandFld)
library(RandomFields)
library(plotly)
library(maps)
library(sf) 
library(leaflet)
library(stringr)
library(ggvoronoi)
library(tibble)
library(tidycensus)

## Convert date
Get_Year <- function(Date){
  as.numeric(strsplit(strsplit(as.character(Date),' ')[[1]][1],'/')[[1]][3])
}

Get_Season <- function(Date){
  ceiling(as.numeric(strsplit(strsplit(as.character(Date),' ')[[1]][1],'/')[[1]][1])/3)
}


data=rawdata%>%filter(!is.na(Latitude) & !is.na(Longitude))%>%#select(Collected.Date,Location.City,Sample.Description,Arsenic,MCL,Latitude,Longitude)%>%
  mutate(Arsenic2=as.numeric(replace(as.character(Arsenic), as.character(Arsenic)=='<0.001', '0.0005')))%>%
  mutate(Y=as.numeric(Arsenic2>=0.01))%>%filter(!is.na(Y))

data$Year <- apply(matrix(data$Collected.Date,,1), MARGIN=1, FUN=Get_Year)
data$Season <- apply(matrix(data$Collected.Date,,1), MARGIN=1, FUN=Get_Season)

for(i in 1:nrow(data)){
  if(data$County[i]=="obrien"){
    print(i)
    data$County[i]="o'brien"
  }
}
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

iowamap=get_map('iowa',zoom=6.5,maptype='roadmap')
iowamap=get_map(c(left=-97,bottom=40,right=-90,top=44))
iowadata=data%>%select(Latitude,Longitude,Y)

ggmap(iowamap)+geom_point(data=iowadata,aes(x=Longitude,y=Latitude,col=as.factor(Y),shape=as.factor(Y)))+
  labs(title="Water Quality (Arsenic)", x ="longitude", y = "latitude")+
  scale_colour_discrete(name  ="Arsenic",breaks=c(0,1),labels=c("< MCL", "Exceed MCL")) +
  scale_shape_discrete(name  ="Arsenic", breaks=c(0,1), labels=c("< MCL", "Exceed MCL"))+theme(legend.text = element_text(size=12))

################################################################
library(igraph)
library(Matrix)
library(deldir)
library(MASS)
library(spatstat)
library(RandomFields)
library(fields)
library(cccd)


## Overview
coords=cbind(data$Latitude,data$Longitude)
#id=which((coords[,1]>40)&(coords[,2]< -85))
y=data$Y
g=nng(coords,k=5)
H=Get_H(g)

# g.dataframe <- igraph::as_data_frame(g)
# from <- g.dataframe$from
# to <- g.dataframe$to
# w <- rep(0,length(from))
# for(i in 1:length(from)){
#   w[i]=fields::rdist(coords[c(from[i],to[i]),])[1,2]
# }
# weight=exp(-w^(1/2))

nlambda=25;#lambda.list=10^seq(-5,5,length=nlambda)
lambda.list=10^seq(-2,2,length=nlambda)
# fit=admm.Bern(y, H, lambda.list=lambda.list, w=weight, maxiter=500)
fit=admm.Bern(y, H, lambda.list=lambda.list, maxiter=1000)

round.digits=2 #or 3
BIC=Get_BIC(fit,y,round.digits=round.digits) 
plot(BIC$BIC)

it=which(BIC$BIC==min(BIC$BIC))
beta.hat=round(fit$x[,it],digits=round.digits)
beta.hat=1/(1+exp(-beta.hat))
data_real <- data.frame("beta" = beta.hat, "lon" = coords[,1], "lat" = coords[,2])
ggplot(data_real, aes(lat, lon)) + geom_point(aes(colour = beta)) +
  scale_color_gradientn(colours = rainbow(4),limits=c(0,1))





##################County-wise analysis####################

IA_As <- get_acs(state = "IA",  geography = "county", 
                 variables = "B19013_001", geometry = TRUE)
Convert_Cty<-function(x){
  tolower(substr(x,start=1,stop=nchar(x)-13))
}
IA_As$County <- apply(matrix(IA_As$NAME,,1), MARGIN=1, FUN=Convert_Cty)


data1=data.frame(County=IA_As$County)
data1$Y=rep(0,nrow(data1));data1$n=rep(0,nrow(data1))
ncounty=nrow(data1)
for(i in 1:nrow(data1)){
  tmp=data%>%filter(County==data1$County[i])
  data1$n[i]=nrow(tmp)
  if(nrow(tmp)!=0){
    data1$Y[i]=sum(tmp$Y)
  }
}
rownames(data1)=data1$County
data1$id=1:ncounty

H=matrix(0,0,ncounty)
from=c();to=c();lon=c();lat=c()
for(i in 1:ncounty){
  v1=data1$County[i]
  for(v2 in NB[[as.character(v1)]]){
    j=data1[v2,]$id
    if(j>i){
      tmp=rep(0,ncounty)
      tmp[i]=1;tmp[j]=-1
      H=rbind(H,tmp)
      from=c(from,i)
      to=c(to,j)
    }
  }
  coord=attributes(IA_As$geometry[i])$bbox
  lon=c(lon,mean(coord[c(1,3)]))
  lat=c(lat,mean(coord[c(2,4)]))
}
data1$lon=lon;data1$lat=lat

es=data1$Y/data1$n; 
es[which(is.na(es))]=0
es=log((es+0.01)/(1.01-es))
w=exp(-abs(as.vector(H%*%es)))

nlambda=40;#lambda.list=10^seq(-5,5,length=nlambda)
lambda.list=10^seq(-4,4,length=nlambda)
fit=admm.MN(N=data1$n,y=data1$Y,H=H,lambda.list,w=w,maxiter=500)


weight=1/(as.vector(H%*%fit$x[,16]))^2
fit_ad=admm.MN(N=data1$n,y=data1$Y,H=H,lambda.list,w=weight,maxiter=500)







#################Delete Duplicate Points#####################
Delete_Duplicate<-function(data){
  lon=lat=Y=county=c()
  for(i in 1:nrow(data)){
    s1=data$Longitude[i];s2=data$Latitude[i]
    id=which((lon==s1)&(lat==s2))
    if(setequal(id,integer(0))){
      lon=c(lon,s1);lat=c(lat,s2);
      Y=c(Y,data$Y[i]);county=c(county,data$County[i])
    }else{
      Y[id]=min(Y[id]+data$Y[i],1)
    }
  }
  return(data.frame(lon=lon,lat=lat,Y=Y,County=county))
}

data2=Delete_Duplicate(data)

g=nng(cbind(data2$lon,data2$lat),k=5)
H=Get_H(g)
nlambda=50;#lambda.list=10^seq(-5,5,length=nlambda)
lambda.list=10^seq(-4,4,length=nlambda)
fit2=admm.Bern(data2$Y, H, lambda.list=lambda.list, maxiter=500)


weight=1/(as.vector(H%*%fit2$x[,18]))^2
fit2_ad=admm.Bern(data2$Y, H, lambda.list=lambda.list,w=weight, maxiter=1000)


data_real <- data.frame("beta" = 1/(1+exp(-fit2$x[,18])), "lon" = data2$lon, "lat" = data2$lat)
ggplot(data=IA_As) +
  geom_sf() +
  #coord_sf(crs = 26915) +
  geom_point(data=data_real,aes(lon, lat,colour = beta),alpha=1/7)  +
  coord_sf(xlim=c(-97,-90), ylim = c(40, 44)) +
  scale_color_gradientn(colours = rainbow(3),limits=c(0,0.5))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.key.height=unit(2.2,"cm"),
        legend.title=element_text(size=0))



scale_x_longitude <- function(xmin=-180, xmax=180, step=1, ...) {
  ewbrks <- seq(xmin,xmax,step)
  ewlbls <- unlist(lapply(ewbrks, function(x) ifelse(x < 0, paste(x, "W"), ifelse(x > 0, paste(x, "E"),x))))
  return(scale_x_continuous("Longitude", breaks = ewbrks, labels = ewlbls, expand = c(0, 0), ...))
}
scale_y_latitude <- function(ymin=-90, ymax=90, step=0.5, ...) {
  nsbrks <- seq(ymin,ymax,step)
  nslbls <- unlist(lapply(nsbrks, function(x) ifelse(x < 0, paste(x, "S"), ifelse(x > 0, paste(x, "N"),x))))
  return(scale_y_continuous("Latitude", breaks = nsbrks, labels = nslbls, expand = c(0, 0), ...))
}


ia <- map_data("county", "iowa")
mid_range <- function(x) mean(range(x))
seats <- do.call(rbind, lapply(split(ia, ia$subregion), function(d) {
  data.frame(lat = mid_range(d$lat), long = mid_range(d$long), subregion = unique(d$subregion))
}))

ia$subregion[which(ia$subregion=='obrien')]="o'brien"

ggplot(data=data2) +
  #geom_point(data=data_real,aes(lon, lat,colour = beta),alpha=1/7)+
  geom_voronoi(aes(lon,lat,fill=p2),outline=ia)+
  #scale_color_gradientn(colours = rainbow(3),limits=c(0,0.5))+
  scale_fill_viridis_c(option= "inferno",end=0.8,limits=c(0,0.4))+
  #coord_sf(xlim=c(-96.5,-90), ylim = c(40.5, 43.5)) +
  scale_x_longitude(-97,-89.5)+ scale_y_latitude(40, 44) +
  xlab('Logitide') + ylab('Latitude') +
  labs(fill=expression(bold(hat(p)))) +
  geom_polygon(aes(x=long,y=lat,group = group),data=ia ,fill = NA, colour = "grey60") +
  geom_text(aes(x=long,y=lat,label = subregion), data = seats, size = 3.5, angle = 45,colour='white')+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),
        aspect.ratio=0.4695096,
        legend.text=element_text(size=14),legend.key.height=unit(1.6,"cm"),
        legend.title=element_text(size=18),
        panel.background=element_rect(fill="white", colour = "black"))

p.hat=1/(1+exp(-fit2_ad$x[,18]))
id1=which(p.hat<0.1)
id2=which((p.hat>0.1)&(data2$lon>-95.78983))
id3=which((p.hat>0.1)&(data2$lon< -95.78983))

data2$p=p.hat
data2$p2=rep(0,nrow(data2))
data2$p2[id1]=mean(data2$Y[id1]);
data2$p2[id2]=mean(data2$Y[id2]);
data2$p2[id3]=mean(data2$Y[id3])


data2$label=rep(0,nrow(data2))
data2$label[id1]='Cluster 1';
data2$label[id2]='Cluster 2';
data2$label[id3]='Cluster 3'

ggplot(data=data2) +
  geom_voronoi(aes(lon,lat,fill=label,col=label),outline=ia)+
  scale_x_longitude(-97,-89.5)+ scale_y_latitude(40, 44, 1) +
  scale_fill_viridis_d(option= "inferno",end=0.8) +
  scale_color_viridis_d(option= "inferno",end=0.8) +
  xlab('Logitide') + ylab('Latitude') +
  geom_polygon(aes(x=long,y=lat,group = group),data=ia ,fill = NA, colour = "grey60") +
  geom_text(aes(x=long,y=lat,label = subregion), data = seats, size = 3.5, angle = 45,colour='white')+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        aspect.ratio=0.4844984,
        legend.text=element_text(size=12),#legend.key.height=unit(1.6,"cm"),
        legend.title=element_text(size=0),
        panel.background=element_rect(fill="white", colour = "black"))

## Data exploartion
IA_As <- get_acs(state = "IA",  geography = "county", 
                 variables = "B19013_001", geometry = TRUE)
Convert_Cty<-function(x){
  tolower(substr(x,start=1,stop=nchar(x)-13))
}
IA_As$County <- apply(matrix(IA_As$NAME,,1), MARGIN=1, FUN=Convert_Cty)


data3=data.frame(County=IA_As$County)
ncounty=nrow(data3)
data3$Y=rep(0,ncounty);data3$n=rep(0,ncounty);
for(i in 1:ncounty){
  tmp=data2%>%filter(County==as.character(data3$County[i]))
  data3$n[i]=nrow(tmp)
  if(nrow(tmp)!=0){
    data3$Y[i]=sum(tmp$Y)
  }
}
rownames(data3)=data3$County
data3$id=1:ncounty
data3$p=data3$Y/data3$n

ia$n=rep(0,nrow(ia));ia$p=rep(0,nrow(ia))
for(i in 1:nrow(ia)){
  ia$n[i]=data3[ia$subregion[i],]$n
  ia$p[i]=data3[ia$subregion[i],]$p
}

ggplot()+
  geom_polygon(aes(x=long,y=lat,group = group,fill=n2, colour=""),data=ia)+
  scale_fill_viridis_c(option=  "plasma" ,end=0.85,values=c(1,0.4,0.2,0),na.value='grey')+
  scale_color_manual(values=NA) +
  scale_x_longitude(-97,-89.5)+ scale_y_latitude(40, 44) +
  xlab('Logitide') + ylab('Latitude') +
  labs(fill='Number of Tests') +
  geom_text(aes(x=long,y=lat,label = subregion), data = seats, size = 3.5, angle = 45,colour='white')+
  guides(colour=guide_legend("Zero Test", override.aes=list(colour="grey"))) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),
        aspect.ratio=0.4695096,
        legend.text=element_text(size=14),legend.key.height=unit(1,"cm"),
        legend.title=element_text(size=16),
        panel.background=element_rect(fill="white", colour = "black"))



ggplot() +
  geom_polygon(data=ia, aes(long, lat, group=group, fill=p, colour="")) +
  scale_fill_viridis_c(option=  "viridis" ,end=0.9,values=c(0,0.15,0.3,1),na.value='grey')+
  scale_colour_manual(values=NA) +   
  scale_x_longitude(-97,-89.5)+ scale_y_latitude(40, 44) +
  xlab('Logitide') + ylab('Latitude') +
  labs(fill=expression(bold(hat(p)))) +
  geom_text(aes(x=long,y=lat,label = subregion), data = seats, size = 3.5, angle = 45,colour='white')+
  guides(colour=guide_legend("No data", override.aes=list(colour="grey"))) +
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),
        aspect.ratio=0.4695096,
        legend.text=element_text(size=14),
        legend.title=element_text(size=16),legend.key.height=unit(1,"cm"),
        panel.background=element_rect(fill="white", colour = "black"))

ggplot(aes(x=n,y=..count..),data=data3)+
  geom_histogram(stat='bin',fill="#69b3a2", breaks=seq(0,650,by=50),color="#e9ecef", alpha=0.9)+
  theme_bw()+
  geom_text(aes(label=as.character(..count..)),stat="bin",breaks=seq(0,650,by=50),boundary = 0,vjust=-0.5) +
  labs(x="the number of tests in each county",y="count")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),
        panel.background=element_rect(fill="white", colour = "black"))

n=data3$n
h=data.frame(n=c('[0,5)','[5,10)','[10,20)','[20,40)','[40,70)','[70,100)','[100,200)','[200,300)','[300,500)','[500,700)'),
             y=c(length(which(n<5 & n>=0)),length(which(n<10 & n>=5)),length(which(n<20 & n>=10)),length(which(n<40 & n>=20)),
               length(which(n<70 & n>=40)),length(which(n<100 & n>=70)),length(which(n<100 & n>=200)),length(which(n<300 & n>=200)),
               length(which(n<500 & n>=300)),length(which(n<700 & n>=500))))

ggplot(h, aes(x=factor(n,level=n), y=y)) + 
  geom_bar(stat="identity")+geom_text(aes(label=y), size=5,vjust=-0.2) +
  labs(x="the number of existing tests in an individual county",y="count")+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16),
        #axis.text.x=element_text(angle=15),
        panel.background=element_rect(fill="white", colour = "black"))


ggmap(iowamap,maprange=TRUE)+
  geom_point(data=data2,aes(x=lon,y=lat,col=as.factor(Y),shape=as.factor(Y)),size=1.25)+
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44) +
  labs(title="Water Quality (Arsenic)", x ="Longitude", y = "Latitude")+
  scale_colour_discrete(name  ="Arsenic",breaks=c(0,1),labels=c("< MCL", "Exceed MCL")) +
  scale_shape_discrete(name  ="Arsenic", breaks=c(0,1), labels=c("< MCL", "Exceed MCL"))+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia ,fill = NA, colour = "grey60",size=0.25,alpha=0.5) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
          aspect.ratio=0.6153846,
          legend.text=element_text(size=12),
          legend.title=element_text(size=14))




H=matrix(0,0,ncounty)
from=c();to=c();lon=c();lat=c()
for(i in 1:ncounty){
  v1=data3$County[i]
  for(v2 in NB[[as.character(v1)]]){
    j=data1[v2,]$id
    if(j>i){
      tmp=rep(0,ncounty)
      tmp[i]=1;tmp[j]=-1
      H=rbind(H,tmp)
      from=c(from,i)
      to=c(to,j)
    }
  }
  coord=attributes(IA_As$geometry[i])$bbox
  lon=c(lon,mean(coord[c(1,3)]))
  lat=c(lat,mean(coord[c(2,4)]))
}
data3$lon=lon;data3$lat=lat

es=data3$Y/data3$n; 
es[which(is.na(es))]=0
es=log((es+0.01)/(1.01-es))
w=exp(-abs(as.vector(H%*%es)))

nlambda=50;#lambda.list=10^seq(-5,5,length=nlambda)
lambda.list=10^seq(-4,4,length=nlambda)
fit=admm.MN(N=data3$n,y=data3$Y,H=H,lambda.list,maxiter=1000)


weight=1/(as.vector(H%*%fit$x[,16]))^2
fit_ad=admm.MN(N=data1$n,y=data1$Y,H=H,lambda.list,w=weight,maxiter=500)






IA_As %>%
  ggplot(aes(fill = 1/(1+exp(-fit$x[,20])))) + 
  geom_sf() + 
  coord_sf(crs = 26915) +  
  scale_fill_gradientn(colours = rainbow(5),limits=c(0,0.5),name=NULL)





## Candidate Well
rawwell=read.csv('/Users/yinlihao/Downloads/Iowa Private Well Water Quality/Wellinfo.csv'
                 ,na.strings = '')

welldata=rawwell%>%filter(!is.na(fLatitude) & !is.na(fLongitude))
testPoints <- welldata%>% dplyr::select(fLatitude,fLongitude)
testPoints <- st_as_sf(testPoints, coords = c("fLongitude", "fLatitude"), crs = st_crs(US_state))
welldata$State=st_join(testPoints, US_state)$ID
US_county <- st_as_sf(maps::map("county", plot = FALSE, fill = TRUE)) 
welldata$County=st_join(testPoints, US_county)$ID

welldata=welldata%>%filter(State=='iowa')%>%
  separate(County,c('State Name','County'),sep=',')%>%
  dplyr::select(fLatitude,fLongitude,State,County)


for(i in 1:nrow(welldata)){
  if(welldata$County[i]=="obrien"){
    welldata$County[i]="o'brien"
  }
}


library(spatstat)
library(maps)
library(maptools)

iamap <- maps::map('state',region='iowa',fill=TRUE, col="transparent", plot=FALSE)
iapoly <- map2SpatialPolygons(iamap, IDs=iamap$names)
spatstat.options(checkpolygons=FALSE)
iaowin <- as.owin.SpatialPolygons(iapoly)
pts <- as.ppp(cbind(welldata$fLongitude,welldata$fLatitude), W=iaowin)
marks(pts)<-NULL
Window(pts) <- iaowin
dsty <- density.ppp(pts,sigma=0.075,kernel='quartic',diggle=TRUE)

pts2 <- as.ppp(cbind(data2$lon,data2$lat), W=iaowin)
marks(pts2)<-NULL
Window(pts2) <- iaowin
dsty2 <- density.ppp(pts2,sigma=0.075,kernel='quartic',diggle=TRUE)

idy=cut(welldata$fLatitude,breaks=seq(dsty$yrange[1],dsty$yrange[2],length=129),labels=FALSE)
idx=cut(welldata$fLongitude,breaks=seq(dsty$xrange[1],dsty$xrange[2],length=129),labels=FALSE)
welldata$intensity=dsty$v[cbind(idy,idx)]
welldata$intensity2=dsty2$v[cbind(idy,idx)]
welldata=welldata%>%filter((!is.na(intensity))& (!is.na(intensity2)))

dxy1 = deldir(data2$lon,data2$lat)
welldata$clust=NA;welldata$id=NA
l=c()
for(i in 1:nrow(data2)){
  print(i)
  plyg=dxy1$dirsgs%>%filter(ind1==i | ind2==i)
  plyg$A=plyg$y2-plyg$y1;plyg$B=plyg$x1-plyg$x2
  plyg$C=(plyg$x2-plyg$x1)*plyg$y1-(plyg$y2-plyg$y1)*plyg$x1
  a=plyg$A*data2$lon[i]+plyg$B*data2$lat[i]+plyg$C
  test=0
  for(j in 1:length(a)){
    test=(((plyg$A[j]*welldata$fLongitude+plyg$B[j]*welldata$fLatitude+plyg$C[j])*a[j]) <= 0)+test
    
  }
  c=which(test==0)
  l=c(l,c)
  welldata$id[c]=i
}

welldata=welldata%>%filter(!is.na(id))
welldata$clust=data2$label[welldata$id]
welldata$IDs=1:nrow(welldata)
#welldata$isSelect=0




### Cluster boundary
lon=c();lat=c();group=c();count=0
for(i in 1:nrow(dxy1$dirsgs)){
  i1=dxy1$dirsgs$ind1[i];i2=dxy1$dirsgs$ind2[i]
  if(data2$label[i1]!=data2$label[i2]){
    count=count+1
    lon=c(lon,dxy1$dirsgs$x1[i],dxy1$dirsgs$x2[i])
    lat=c(lat,dxy1$dirsgs$y1[i],dxy1$dirsgs$y2[i])
    group=c(group,count,count)
  }
}

ia_state <- map_data("state", "iowa")
p_ia <- st_polygon(list(cbind(ia_state$long,ia_state$lat)))

long=lon[which(lon< -95.75)];lati=lat[which(lon< -95.75)]
it=23;lon3=c(long[it],long[it+1]);lat3=c(lati[it],lati[it+1]);
v1=long[it+1];v2=lati[it+1];
long=long[-c(it,it+1)];lati=lati[-c(it,it+1)];
while(length(long)>0){
  it=which(long==v1&lati==v2)
  if(it%%2==1){
    lon3=c(lon3,long[it+1]);lat3=c(lat3,lati[it+1]);
    v1=long[it+1];v2=lati[it+1];
    long=long[-c(it,it+1)];lati=lati[-c(it,it+1)];
  }else{
    lon3=c(lon3,long[it-1]);lat3=c(lat3,lati[it-1]);
    v1=long[it-1];v2=lati[it-1];
    long=long[-c(it-1,it)];lati=lati[-c(it-1,it)];
  }
}
p3 = st_polygon(list(cbind(c(lon3,lon3[1]),c(lat3,lat3[1]))))
p3 = st_intersection(p3,p_ia)

long=lon[which(lon> -95.75)];lati=lat[which(lon> -95.75)]
it=293;lon2=c(long[it],long[it+1]);lat2=c(lati[it],lati[it+1]);
v1=long[it+1];v2=lati[it+1];
long=long[-c(it,it+1)];lati=lati[-c(it,it+1)];
while(length(long)>0){
  it=which(long==v1&lati==v2)
  if(it%%2==1){
    lon2=c(lon2,long[it+1]);lat2=c(lat2,lati[it+1]);
    v1=long[it+1];v2=lati[it+1];
    long=long[-c(it,it+1)];lati=lati[-c(it,it+1)];
  }else{
    lon2=c(lon2,long[it-1]);lat2=c(lat2,lati[it-1]);
    v1=long[it-1];v2=lati[it-1];
    long=long[-c(it-1,it)];lati=lati[-c(it-1,it)];
  }
}
p2 = st_polygon(list(cbind(c(lon2,lon2[1]),c(lat2,lat2[1]))))
p2 = st_intersection(p2,p_ia)

p1 = st_difference(st_difference(p_ia,p2),p3)

a1=st_area(p1);a2=st_area(p2);a3=st_area(p3)
ClustMap <- data.frame(lon=c(p1[[1]][,1],p2[[1]][,1],p3[[1]][,1]),
                       lat=c(p1[[1]][,2],p2[[1]][,2],p3[[1]][,2]),
                       group=c(rep('Cluster 1',nrow(p1[[1]])),rep('Cluster 2',nrow(p2[[1]])),rep('Cluster 3',nrow(p3[[1]]))))
ClustMap1 <- ClustMap[which(ClustMap$group=='Cluster 1'),]
ClustMap2 <- ClustMap[which(ClustMap$group=='Cluster 2'),]
ClustMap3 <- ClustMap[which(ClustMap$group=='Cluster 3'),]


### Sampling
n1=12446;n2=1444;n3=743;
i1=1190;i2=305;i3=6145
welldata$intensity3=NA
welldata$intensity3[which(welldata$clust=='Cluster 1')]=i1
welldata$intensity3[which(welldata$clust=='Cluster 2')]=i2
welldata$intensity3[which(welldata$clust=='Cluster 3')]=i3

welldata$isSelect=0
samplewell=welldata%>%filter(clust=='Cluster 1')
isSelect=runif(nrow(samplewell)) <(i1/samplewell$intensity)
welldata$isSelect[samplewell$IDs[which(isSelect==1)]]=1

samplewell=welldata%>%filter(clust=='Cluster 2')
isSelect=runif(nrow(samplewell)) <(i2/samplewell$intensity)
welldata$isSelect[samplewell$IDs[which(isSelect==1)]]=1

samplewell=welldata%>%filter(clust=='Cluster 3')
isSelect=runif(nrow(samplewell)) <(i3/samplewell$intensity)
welldata$isSelect[samplewell$IDs[which(isSelect==1)]]=1

samplewell=welldata%>%filter(isSelect==1)
ggmap(iowamap,maprange=TRUE)+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia_state ,fill = NA, colour = "grey60",size=0.25,alpha=0.5) +
  geom_polygon(aes(x=lon,y=lat,group=1,size='Cluster 1'),fill=NA,data=ClustMap1,colour='green')+
  geom_polygon(aes(x=lon,y=lat,group=2,size='Cluster 2'),fill=NA,data=ClustMap2,colour='blue')+
  geom_polygon(aes(x=lon,y=lat,group=3,size='Cluster 3'),fill=NA,data=ClustMap3,colour='red')+
  scale_size_manual(name='Clusters',values=c(1,1.01,1.02),guide=guide_legend(override.aes = list(colour=c("green","blue","red")))) +
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44) +
  labs(title="Water Quality (Arsenic)", x ="Longitude", y = "Latitude")+
  geom_point(data=samplewell,aes(x=fLongitude,y=fLatitude,colour=isSelect),show.legend=FALSE,size=0.35,alpha=0.5)+
  #scale_colour_discrete(name  ="Wells",breaks=c(0,1),labels=c("Available Candidates", "New Samples")) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        aspect.ratio=0.6153846,
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))


welldata$intensity4=welldata$intensity3-welldata$intensity2
welldata$intensity4[which(welldata$intensity4<0)]=0

welldata$isSelect2=0
samplewell=welldata%>%filter(clust=='Cluster 1')
isSelect=runif(nrow(samplewell)) <(samplewell$intensity4/samplewell$intensity)
welldata$isSelect2[samplewell$IDs[which(isSelect==1)]]=1

samplewell=welldata%>%filter(clust=='Cluster 2')
isSelect=runif(nrow(samplewell)) <(samplewell$intensity4/samplewell$intensity)
welldata$isSelect2[samplewell$IDs[which(isSelect==1)]]=1

samplewell=welldata%>%filter(clust=='Cluster 3')
isSelect=runif(nrow(samplewell)) <(samplewell$intensity4/samplewell$intensity)
welldata$isSelect2[samplewell$IDs[which(isSelect==1)]]=1


samplewell=data.frame(lon=c(welldata$fLongitude[which(welldata$isSelect2==1)],data2$lon),
                      lat=c(welldata$fLatitude[which(welldata$isSelect2==1)],data2$lat),
                      label=c(rep(1,length(which(welldata$isSelect2==1))),rep(0,nrow(data2))))
  
ggmap(iowamap,maprange=TRUE)+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia ,fill = NA, colour = "grey60",size=0.25,alpha=0.5) +
  geom_polygon(aes(x=lon,y=lat,group=1,size='Cluster 1'),fill=NA,data=ClustMap1,colour='green')+
  geom_polygon(aes(x=lon,y=lat,group=2,size='Cluster 2'),fill=NA,data=ClustMap2,colour='blue')+
  geom_polygon(aes(x=lon,y=lat,group=3,size='Cluster 3'),fill=NA,data=ClustMap3,colour='red')+
  scale_size_manual(name='Clusters',values=c(1,1.01,1.02),guide=guide_legend(override.aes = list(colour=c("green","blue","red")))) +
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44) +
  labs(title="Well Sampling Example", x ="Longitude", y = "Latitude")+
  geom_point(data=samplewell,aes(x=lon,y=lat,colour=factor(label)),size=.35,alpha=0.5)+
  scale_colour_discrete(name  ="Wells",breaks=c(0,1),labels=c("Existing Samples", "New Samples")) +
  geom_text(aes(x=long,y=lat,label = subregion), data = seats, size = 3.5, angle = 45,colour='black')+
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        aspect.ratio=0.6153846,
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))  


ggmap(iowamap,maprange=TRUE)+
  geom_point(data=welldata,aes(x=fLongitude,y=fLatitude),size=0.15,alpha=0.25,colour='#56B4E9')+
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44) +
  labs(title="Candidate wells", x ="Longitude", y = "Latitude")+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia ,fill = NA, colour = "grey60",size=0.5,alpha=0.5) +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=14),
        aspect.ratio=0.6153846,
        legend.text=element_text(size=12),
        legend.title=element_text(size=14))


### intensity map
xrange=c(-97.5,-89.5);yrange=c(40, 44);
np=200;lx=(xrange[2]-xrange[1])/np;ly=(yrange[2]-yrange[1])/np
coords.grid=expand.grid(seq(1/(2*np),1-1/(2*np),1/np),seq(1/(2*np),1-1/(2*np),1/np))
coords.grid[,1]=coords.grid[,1]*(xrange[2]-xrange[1])+xrange[1]
coords.grid[,2]=coords.grid[,2]*(yrange[2]-yrange[1])+yrange[1]

idy=cut(coords.grid[,2],breaks=c(-1000,seq(dsty$yrange[1],dsty$yrange[2],length=129)[2:128],1000),labels=FALSE)
idx=cut(coords.grid[,1],breaks=c(-1000,seq(dsty$xrange[1],dsty$xrange[2],length=129)[2:128],1000),labels=FALSE)
intensity=dsty$v[cbind(idy,idx)]

lon=c();lat=c();group=c();v=c()

for(i in 1:nrow(coords.grid)){
  v1=coords.grid[i,1];v2=coords.grid[i,2]
  p4 = st_polygon(list(cbind(c(v1-lx/2,v1-lx/2,v1+lx/2,v1+lx/2,v1-lx/2),
                             c(v2-ly/2,v2+ly/2,v2+ly/2,v2-ly/2,v2-ly/2))))
  p4 = st_intersection(p4,p_ia)
  if(!is.empty(p4)){
    if(class(p4[[1]])=='list'){
      p4=p4[[1]]
    }
    lon=c(lon,p4[[1]][,1]);lat=c(lat,p4[[1]][,2])
    print(nrow(p4[[1]]))
    group=c(group,rep(i,nrow(p4[[1]])))
    if(is.na(intensity[i])){
      d=(v1-coords.grid[,1])^2+(v2-coords.grid[,2])^2
      insty=intensity[order(d)]
      v=c(v,rep(insty[min(which(!is.na(insty)))],nrow(p4[[1]])))
    }else{
      v=c(v,rep(intensity[i],nrow(p4[[1]])))
    }
  }
}

IntstyMap<-data.frame(lon=lon,lat=lat,group=group,v=v)

A=56272;B=a1+a2+a3

IntstyMap$v=IntstyMap$v/A*B

ggmap(iowamap,maprange=TRUE)+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia_state ,fill = NA, colour = "grey60",size=0.25,alpha=0.5) +
  geom_polygon(aes(x=lon,y=lat,group = group,fill=v),data=IntstyMap) +
  scale_fill_viridis_c(option= "plasma",name='Intensity',trans='log',breaks=c(1,4.5,20)) +
  geom_polygon(aes(x=lon,y=lat,group = group,col=v),data=IntstyMap,fill=NA) +
  scale_color_viridis_c(option= "plasma",name='Intensity',trans='log',breaks=c(1,4.5,20)) +
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44, 1) +
  labs(title='Intensity of candidate wells', x ="Longitude", y = "Latitude")+
  theme(title=element_text(size=16),
        axis.text=element_text(size=14),axis.title=element_text(size=16),
        aspect.ratio=0.6153846,
        legend.text=element_text(size=14),
        legend.key.height=unit(2,"cm"),
        legend.title=element_text(size=16))


xrange=c(-97.5,-89.5);yrange=c(40, 44);
np=200;lx=(xrange[2]-xrange[1])/np;ly=(yrange[2]-yrange[1])/np
coords.grid=expand.grid(seq(1/(2*np),1-1/(2*np),1/np),seq(1/(2*np),1-1/(2*np),1/np))
coords.grid[,1]=coords.grid[,1]*(xrange[2]-xrange[1])+xrange[1]
coords.grid[,2]=coords.grid[,2]*(yrange[2]-yrange[1])+yrange[1]

idy=cut(coords.grid[,2],breaks=c(-1000,seq(dsty$yrange[1],dsty$yrange[2],length=129)[2:128],1000),labels=FALSE)
idx=cut(coords.grid[,1],breaks=c(-1000,seq(dsty$xrange[1],dsty$xrange[2],length=129)[2:128],1000),labels=FALSE)
intensity=dsty2$v[cbind(idy,idx)]

lon=c();lat=c();group=c();v=c();vv1=c();vv2=c()
ii=c(i1,i2,i3)

for(i in 1:nrow(coords.grid)){
  v1=coords.grid[i,1];v2=coords.grid[i,2]
  p4 = st_polygon(list(cbind(c(v1-lx/2,v1-lx/2,v1+lx/2,v1+lx/2,v1-lx/2),
                             c(v2-ly/2,v2+ly/2,v2+ly/2,v2-ly/2,v2-ly/2))))
  if(!is.empty(st_intersection(p4,p1))){
    s2=ii[1]
  }else if(!is.empty(st_intersection(p4,p2))){
    s2=ii[2]
  }else{
    s2=ii[3]
  }
  p4 = st_intersection(p4,p_ia)
  if(!is.empty(p4)){
    if(class(p4[[1]])=='list'){
      p4=p4[[1]]
    }
    lon=c(lon,p4[[1]][,1]);lat=c(lat,p4[[1]][,2])
    print(s2)
    group=c(group,rep(i,nrow(p4[[1]])))
    if(is.na(intensity[i])){
      d=(v1-coords.grid[,1])^2+(v2-coords.grid[,2])^2
      insty=intensity[order(d)]
      s1=insty[min(which(!is.na(insty)))]
    }else{
      s1=intensity[i]
    }
    v=c(v,rep(max(s2-s1,0),nrow(p4[[1]])))
    vv1=c(vv1,rep(s1,nrow(p4[[1]])))
    vv2=c(vv2,rep(s2,nrow(p4[[1]])))
  }
}

IntstyMap2<-data.frame(lon=lon,lat=lat,group=group)
IntstyMap2$v1<-vv1
IntstyMap2$v2<-vv2
IntstyMap2$v3<-v

IntstyMap2$v1=IntstyMap2$v1/A*B
IntstyMap2$v2=IntstyMap2$v2/A*B
IntstyMap2$v3=IntstyMap2$v3/A*B

IntstyMap2$v1[which(IntstyMap2$v1<0)]=0
IntstyMap2$v2[which(IntstyMap2$v2<0)]=0
IntstyMap2$v3[which(IntstyMap2$v3<0)]=0

ggmap(iowamap,maprange=TRUE)+
  geom_polygon(aes(x=long,y=lat,group = group),data=ia_state ,fill = NA, colour = "grey60",size=0.25,alpha=0.5) +
  geom_polygon(aes(x=lon,y=lat,group = group,fill=v1),data=IntstyMap2) +
  scale_fill_viridis_c(option= "plasma",name='Intensity',trans='sqrt',breaks=c(0.03,0.35,1,2)) +
  geom_polygon(aes(x=lon,y=lat,group = group,col=v1),data=IntstyMap2,fill=NA) +
  scale_color_viridis_c(option= "plasma",name='Intensity',trans='sqrt',breaks=c(0.03,0.35,1,2)) +
  xlim(-97.5,-89.5) + ylim(40, 44) +
  scale_x_longitude(-97.5,-89.5)+ scale_y_latitude(40, 44, 1) +
  labs(title='Intensity of existing wells', x ="Longitude", y = "Latitude")+
  theme(title=element_text(size=16),
        axis.text=element_text(size=14),axis.title=element_text(size=16),
        aspect.ratio=0.6153846,
        legend.text=element_text(size=14),
        legend.key.height=unit(2,"cm"),
        legend.title=element_text(size=16))


## Delete duplicate
output=Delete_Duplicate(coords,Y=y)

g=nng(output$coords,k=5)
H=Get_H(g)

nlambda=25;#lambda.list=10^seq(-5,5,length=nlambda)
lambda.list=10^seq(-2,2,length=nlambda)
# fit=admm.Bern(y, H, lambda.list=lambda.list, w=weight, maxiter=500)
fit=admm.Bern(output$Y, H, lambda.list=lambda.list, maxiter=1000)

round.digits=2 #or 3
BIC=Get_BIC(fit,output$Y,round.digits=round.digits) 
plot(BIC$BIC)

it=which(BIC$BIC==min(BIC$BIC))
beta.hat=round(fit$x[,it],digits=round.digits)
beta.hat=1/(1+exp(-beta.hat))
data_real <- data.frame("beta" = beta.hat, "lon" = output$coords[,1], "lat" = output$coords[,2])
ggplot(data_real, aes(lat, lon)) + geom_point(aes(colour = beta)) +
  scale_color_gradientn(colours = rainbow(4),limits=c(0,1))




## test
# n=1000
# lat=runif(n);lon=runif(n);coords=cbind(lat,lon)
# ClustType=function(x,y,type=1){ 
#   if(type==1){
#     mclust=cut(x+y,breaks=c(-10000,1,10000),labels=FALSE);
#   }
#   return(mclust)
# }
# a=c(0.2,0.8)
# mclust=ClustType(lon,lat,type=1)
# beta=a[mclust]
# beta.true=beta
# y=rep(0,n)
# for(x in unique(beta)){
#   y[which(beta==x)]=rbinom(sum(beta==x),1,x)
# }
#  
# g <- nng(coords,para=4)
# H=Get_H(g)
# 
# nlambda=20;lambda.list=10^seq(-1,5,length=nlambda)
# fit1=admm.Bern(y, H, lambda.list=lambda.list, maxiter=500, beta0=beta.true)
# 
# mse=c()
# for(i in 1:length(fit1$lambda)){
#   beta.hat=fit1$x[,i]
#   beta.hat=1/(1+exp(-beta.hat))
#   mse=c(mse,mean((beta.hat-beta.true)^2))
# }
# 
# it=min(which(mse==min(mse)))
# beta.hat=fit1$x[,it]
# beta.hat=1/(1+exp(-beta.hat))
# data_hat <- data.frame("beta" = beta.hat, "lon" = lat, "lat" = lon)
# ggplot(data_hat, aes(lon, lat)) + 
#   geom_point(aes(colour = beta)) +
#   scale_fill_gradient(limits=c(0,1))
# 
# 
# data_true <- data.frame("beta" = beta.true, "lon" = lat, "lat" = lon)
# ggplot(data_true, aes(lon, lat)) + 
#   geom_point(aes(colour = beta)) +
#   scale_fill_gradient(limits=c(0,1))

IA_As <- get_acs(state = "IA",  geography = "county", 
                     variables = "B19013_001", geometry = TRUE)

Convert_Cty<-function(x){
  tolower(substr(x,start=1,stop=nchar(x)-13))
}
IA_As$County <- apply(matrix(IA_As$NAME,,1), MARGIN=1, FUN=Convert_Cty)
IA_As$count <- IA_As$estimate <-rep(0,nrow(IA_As))

for(i in 1:nrow(IA_As)){
  id=which(data$County==IA_As$County[i])
  IA_As$count[i]=length(id)
  if(length(id)!=0){
    IA_As$estimate[i]=mean(beta.hat[id])
  }
  print(c(i,length(id),IA_As$estimate[i]))
}

IA_As %>%
  ggplot(aes(fill = 1/(1+exp(-fit$x[,25])))) + 
  geom_sf() + 
  coord_sf(crs = 26915) +  
  scale_fill_gradientn(colours = rainbow(5),limits=c(0,0.5),name=NULL)


data_real <- data.frame("beta" = p.hat[id2], "lon" = data2$lon[id2], "lat" = data2$lat[id2])
ggplot(data=IA_As) +
  geom_sf() +
  #coord_sf(crs = 26915) +
  geom_point(data=data_real,aes(lon, lat,colour = beta),alpha=1/7)  +
  #coord_sf(xlim=c(-97,-90), ylim = c(40, 44)) +
  scale_color_gradientn(colours = rainbow(3),limits=c(0,0.3))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=18),
        legend.text=element_text(size=14),legend.key.height=unit(2.4,"cm"),
        legend.title=element_text(size=0))









