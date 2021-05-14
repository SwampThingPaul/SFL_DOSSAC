## 
## Proposed DO SSAC - Verification
##
## Code was compiled by Paul Julian
## contact info: pjulian@sccf.org

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
# Data Wrangling
library(AnalystHelper);
library(plyr)
library(reshape)
library(openxlsx)
library(zoo)

# GIS libraries 
library(rgdal)
library(rgeos)
library(raster)
library(tmap)
library(ceramic)

# GAM 
library(mgcv)
library(gratia)
library(EnvStats)

# Paths
wd="C:/Julian_LaCie/_GitHub/SFL_DOSSAC"
paths=c(paste0(wd,c("/Exports/","/Plots/")),"C:/Julian_LaCie/_Data/FDEP/iwr2020_run60_2020-10-05/")
# Folder.Maker(paths);#One and done. Creates folders in working directory.
export.path=paths[1]
plot.path=paths[2]
data.path=paths[3]

GIS.path="C:/Julian_LaCie/_GISData"

## add mapbox token here
## See Ceramic package for info
# public.token=" "
# Sys.setenv(MAPBOX_API_KEY=public.token)

# Helper variables
nad83.pro=CRS("+init=epsg:4269")
utm17=CRS("+init=epsg:26917")
wgs84=CRS(SRS_string = "EPSG:4326")

tmap_mode("view")

### Functions
# from descdist() function in fitdistrplus
moment <- function(data, k) {
        m1 <- mean(data)
        return(sum((data - m1)^k)/length(data))
}
skewness <- function(data) {
        sd <- sqrt(moment(data, 2))
        n <- length(data)
        gamma1 <- moment(data, 3)/sd^3
        unbiased.skewness <- sqrt(n * (n - 1)) * gamma1/(n - 
                                                                 2)
        return(unbiased.skewness)
}
kurtosis <- function(data) {
        n <- length(data)
        var <- moment(data, 2)
        gamma2 <- moment(data, 4)/var^2
        unbiased.kurtosis <- (n - 1)/((n - 2) * (n - 3)) * 
                ((n + 1) * gamma2 - 3 * (n - 1)) + 3
        return(unbiased.kurtosis)
}
standdev <- function(data) {
        sd(data)
}
dist.beta=function(kurtosis.boot){
        kurtmax <- max(10, ceiling(max(kurtosis.boot)))
        
        p <- exp(-100)
        lq <- seq(-100, 100, 0.1)
        q <- exp(lq)
        s2a <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 * 
                                                      p * q)
        ya <- kurtmax - (3 * (p + q + 1) * (p * q * (p + 
                                                             q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) * 
                                                                                              (p + q + 3)))
        p <- exp(100)
        lq <- seq(-100, 100, 0.1)
        q <- exp(lq)
        s2b <- (4 * (q - p)^2 * (p + q + 1))/((p + q + 2)^2 * 
                                                      p * q)
        yb <- kurtmax - (3 * (p + q + 1) * (p * q * (p + 
                                                             q - 6) + 2 * (p + q)^2)/(p * q * (p + q + 2) * 
                                                                                              (p + q + 3)))

        return(data.frame(s2a=s2a,s2b=s2b,ya=ya,yb=yb))
}
dist.lognorm=function(kurtmax){
lshape <- seq(-100, 100, 0.1)
shape <- exp(lshape)
es2 <- exp(shape^2)
s2 <- (es2 + 2)^2 * (es2 - 1)
y <- kurtmax - (es2^4 + 2 * es2^3 + 3 * es2^2 - 
                        3)
return(data.frame(s2=s2,y=y))
}

# GIS ---------------------------------------------------------------------
public.token="pk.eyJ1IjoicGp1bGlhbiIsImEiOiJjanllbmJ0eXkxMzV0M2dzNXh5NGRlYXdqIn0.g4weKGOt1WdNZLg2hxBz1w"
Sys.setenv(MAPBOX_API_KEY=public.token)


roads.all=spTransform(readOGR(paste0(GIS.path,"/FDOT"),"FDOT_Roads"),utm17)

wbids=spTransform(readOGR(paste0(GIS.path,"/FDEP"),"WBIDs"),wkt(utm17))

# -------------------------------------------------------------------------
dat.qual=data.frame(QUALIFIER=c(NA,"!","A","D","E","F","I","R","T","U","*","?","B","H","J","K","L","M","N","O","Q","V","Y","Z"),
                    FATALYN=c("N",rep("N",9),rep("Y",14)))


# WBID3240F - Daughtrey ---------------------------------------------------

dat3240F=read.xlsx(paste0(data.path,"periodOfRecordData3240F.xlsx"))
dat3240F$date=date.fun(convertToDate(dat3240F$date))

dat3240F=subset(dat3240F,year%in%seq(2006,2020,1))
unique(dat3240F$station.id)
dat3240F.sites=ddply(dat3240F,c("station.id","lat","long"),summarise,N.val=N.obs(year))
dat3240F.shp=SpatialPointsDataFrame(dat3240F.sites[,c("long","lat")],data=dat3240F.sites[,1:3],proj4string=nad83.pro)
dat3240F.shp=spTransform(dat3240F.shp,utm17)
tm_shape(dat3240F.shp)+tm_dots()

## Map
wbid.3240F=subset(wbids,WBID=="3240F")

###
roi=extent(spTransform(gBuffer(wbid.3240F,width=5000),wgs84))
im <- cc_location(roi,zoom=13)
im2=projectRaster(im,crs=wkt(utm17))
im2=setValues(im2,scales::rescale(values(im2), c(0,255)))

lower.sites=c("21FLEECOAB96009","112WRD  264140081494600")
roi=extent(spTransform(gBuffer(subset(dat3240F.shp,station.id%in%lower.sites),width=1000),wgs84))
im <- cc_location(roi)
im2.site=projectRaster(im,crs=wkt(utm17))
im2.site=setValues(im2.site,scales::rescale(values(im2.site), c(0,255)))

roi=extent(spTransform(gBuffer(dat3240F.shp,width=3000),wgs84))
im <- cc_location(roi)
im3.site=projectRaster(im,crs=wkt(utm17))
im3.site=setValues(im3.site,scales::rescale(values(im3.site), c(0,255)))

bbox.poly1=as(raster::extent(gBuffer(subset(dat3240F.shp,station.id%in%lower.sites),width=1000)),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly1)=utm17
# png(filename=paste0(plot.path,"WBID3240F_map.png"),width=4.25,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1,1,2,3),2,2,byrow=F),widths=c(1))

bbox.lims=bbox(gBuffer(wbid.3240F,width=2000))
bbox.poly=as(raster::extent(gBuffer(wbid.3240F,width=2000)),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

plotRGB(im2,ext=extent(gBuffer(wbid.3240F,width=2000)))
plot(wbid.3240F,border="grey",col=adjustcolor("white",0.1),lwd=1,add=T)
plot(crop(roads.all,gBuffer(bbox.poly,width=100)),col="grey",lwd=0.8,add=T)
plot(dat3240F.shp,add=T,pch=21,bg="indianred1",cex=1,lwd=0.1)
plot(bbox.poly1,lty=2,lwd=1.5,border="red",add=T)
mapmisc::scaleBar(utm17,"topleft",bty="n",cex=0.8,col="white")

plotRGB(im2.site,ext=extent(gBuffer(subset(dat3240F.shp,station.id%in%lower.sites),width=1000)))
plot(subset(dat3240F.shp,station.id%in%lower.sites),add=T,pch=21,bg="indianred1",cex=1.25)
text(subset(dat3240F.shp,station.id%in%lower.sites),"station.id",halo=T,font=2,pos=3,cex=0.5)

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend=c("WBID 3240F", "Major Roads","Monitoring Locations" ),
       pch=c(22,NA,21),
       lty=c(NA,1,NA),
       lwd=c(1,1,0.5),
       pt.bg=c(adjustcolor("white",0.1),NA,"indianred1"),
       col=c("grey","grey","black"),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()
# ##
dat3240F=subset(dat3240F,station.id!="112WRD 264140081494600")

unique(dat3240F$xcode)
# QA/QC qualifiers 
unique(dat3240F$rcode)
unique(dat3240F$xcode)
quals=as.character(unique(dat3240F$xcode))
spl=strsplit(quals,split="")
quals=data.frame(xcode=quals,q1=sapply(spl,"[",1),q2=sapply(spl,"[",2),q3=sapply(spl,"[",3))
quals$Fatal=with(quals,ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                                q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,"Y","N"))
dat3240F=merge(dat3240F,quals[,c("xcode","Fatal")],"xcode",all.x=T)
dat3240F=subset(dat3240F,Fatal=="N")
dat3240F$HalfMDL=with(dat3240F,ifelse(is.na(mdl)==T,result,ifelse(result<=mdl,result/2,result)))

head(dat3240F)
# head(subset(dat3240F,xcode=="U"))
ddply(dat3240F,c("master.code","parameter"),summarise,N.val=N.obs(HalfMDL))

params.keep=c("CHLAC","COLOR","COND","DO","DOSAT","TN","TP","TOC","TEMP")
dat3240F.xtab=reshape2::dcast(subset(dat3240F,master.code%in%params.keep),wbid+station.id+lat+long+date~master.code,value.var="HalfMDL",mean)
dat3240F.xtab$Sal=with(dat3240F.xtab,SalinityCalc(COND,TEMP))
dat3240F.xtab$DOsat.calc=with(dat3240F.xtab,DO_PerSat(TEMP,DO,Sal))
dat3240F.xtab$DoY=as.numeric(format(dat3240F.xtab$date,"%j"))
dat3240F.xtab$CY=as.numeric(format(dat3240F.xtab$date,"%Y"))
N.obs(dat3240F.xtab$DOsat.calc)
mean(dat3240F.xtab$DOsat.calc,na.rm=T)
sd(dat3240F.xtab$DOsat.calc,na.rm=T)
quantile(dat3240F.xtab$DOsat.calc,na.rm=T,probs=c(0.1,0.25,0.5,0.75))
qnorm(0.1,mean(dat3240F.xtab$DOsat.calc,na.rm=T),sd(dat3240F.xtab$DOsat.calc,na.rm=T))

##
dat3240F$time=with(dat3240F,ifelse(nchar(time)==3,paste0(0,time),time))
dat3240F$datetime=with(dat3240F,date.fun(paste(date,time),form="%F %H%M"))
dat3240F.xtab2=reshape2::dcast(subset(dat3240F,master.code%in%params.keep),wbid+station.id+lat+long+datetime~master.code,value.var="HalfMDL",mean)
dat3240F.xtab2$CY=as.numeric(format(dat3240F.xtab2$date,"%Y"))
dat3240F.xtab2$Sal=with(dat3240F.xtab2,SalinityCalc(COND,TEMP))
dat3240F.xtab2$DOsat.calc=with(dat3240F.xtab2,DO_PerSat(TEMP,DO,Sal))
dat3240F.xtab2$DO.TOD.WQS=with(dat3240F.xtab2,DO.TOD.WQS.stream(datetime))
dat3240F.xtab2$exceed=with(dat3240F.xtab2,ifelse(DOsat.calc<DO.TOD.WQS,1,0))

rslt.3240F=ddply(dat3240F.xtab2,c("station.id","CY"),summarise,N.exceed=sum(exceed,na.rm=T),N.val=N.obs(DOsat.calc))
rslt.3240F$PerExceed=with(rslt.3240F,N.exceed/N.val)*100
rslt.3240F$status=with(rslt.3240F,ifelse(PerExceed>10,1,0))
rslt.3240F
rslt.3240F=ddply(rslt.3240F,"station.id",summarise,sum.status=sum(status),n.val=N.obs(status))

plot(DOsat.calc~datetime,subset(dat3240F.xtab2,station.id=="21FLEECO20-29GR"))
with(subset(dat3240F.xtab2,station.id=="21FLEECO20-29GR"),points(datetime,DO.TOD.WQS,pch=2))
##

library(fitdistrplus)
descdist(subset(dat3240F.xtab,is.na(DOsat.calc)==F)$DOsat.calc,boot=1000)
normal_dist <- fitdist(subset(dat3240F.xtab,is.na(DOsat.calc)==F)$DOsat.calc, "norm")
plot(normal_dist)
gofstat(normal_dist)

tmp=subset(dat3240F.xtab,is.na(DOsat.calc)==F)$DOsat.calc
plotdist(tmp,"norm",para=list(mean=mean(tmp),sd=sd(tmp)))
swtest.rslt=shapiro.test(tmp)
kstest.rslt=ks.test(tmp,"pnorm",mean(tmp),sd(tmp))
adtest.rslt=nortest::ad.test(tmp)
lill.rslt=nortest::lillie.test(tmp)

datdist.3240F=data.frame(Test=c("Shapiro-Wilks","Kolmogorov-Smirnov","Anderson-Darling","Lilliefors"),
           Statistic=as.numeric(c(swtest.rslt$statistic,kstest.rslt$statistic,adtest.rslt$statistic,lill.rslt$statistic)),
           pval=c(swtest.rslt$p.value,kstest.rslt$p.value,adtest.rslt$p.value,lill.rslt$p.value))
datdist.3240F$pval.c=with(datdist.3240F,ifelse(pval<0.01,"$<$0.01",ifelse(pval<0.05,"$<$0.05",format(round(pval,2)))))
datdist.3240F$Statistic=round(datdist.3240F$Statistic,2)

# png(filename=paste0(plot.path,"WBID3240F_normDist.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:4,2,2,byrow=T))
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.5));

#hist
h=hist(tmp, plot=F)
xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
yhist <- dnorm(xhist,mean(tmp),sd(tmp))
xlim.val=range(h$breaks);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,max(yhist)+max(yhist)*0.1);by.y=0.005;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(tmp, freq = FALSE,ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,yaxs="i")
lines(xhist,yhist, lty = 1,col="red",lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"Denisty")
mtext(side=1,line=1.75,"DO (% Sat)")

#QQplot
x.val=qnorm(ppoints(sort(tmp)),mean(tmp),sd(tmp))
y.val=sort(tmp)
xlim.val=c(-20,110);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Quantile")
mtext(side=1,line=1.75,"Theoretical Quantile")

#CDFplot
x.val=y.val
y.val=ppoints(x.val)
xlim.val=c(0,120);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
sfin=seq(min(h$breaks),max(h$breaks),length.out=100)
lines(sfin,pnorm(sfin,mean(tmp),sd(tmp)),col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"CDF")
mtext(side=1,line=1.75,"DO (%Sat)")

#PPplot
x.val=pnorm(sort(tmp),mean(tmp),sd(tmp))
y.val=ppoints(sort(tmp))
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Probabilities")
mtext(side=1,line=1.75,"Theoretical Probabilities")
dev.off()

boot=1000
data=subset(dat3240F.xtab,is.na(DOsat.calc)==F)$DOsat.calc
res <- list(min = min(data), max = max(data), median = median(data), 
            mean = mean(data), sd = standdev(data), skewness = skewness(data), 
            kurtosis = kurtosis(data), method = "unbiased")
n <- length(data)
databoot <- matrix(sample(data, size = n * boot,replace = TRUE), nrow = n, ncol = boot)
s2boot <- sapply(1:boot, function(iter) skewness(databoot[,iter])^2)
kurtboot <- sapply(1:boot, function(iter) kurtosis(databoot[,iter]))
kurtmax <- max(10, ceiling(max(kurtboot)))

# png(filename=paste0(plot.path,"WBID3240F_desdist.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
layout(matrix(1:2,1,2),widths=c(1,0.5))
par(family="serif",mar=c(1,2,1,0.5),oma=c(2,1,0.25,0.5));

xmax <- max(4, ceiling(max(s2boot)))
xlim.val=c(0, xmax);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ymax <- kurtmax - 1
ylim.val=c(0, ymax)
plot(res$skewness^2, kurtmax - res$kurtosis, xlim=xlim.val, ylim=ylim.val,type="n",axes=F,ann=F)
abline(h=0:ymax,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
with(dist.beta(kurtboot),polygon(c(s2a,s2b), c(ya,yb), col = "lightgrey", border = "lightgrey"))# beta
with(dist.beta(kurtboot),lines(4/exp(seq(-100,100,0.1)),yb,lty=2))# gamma
with(dist.lognorm(kurtmax),lines(s2,y,lty=3))# lognormal
points(s2boot, kurtmax - kurtboot, pch = 21, col = "orange",cex = 0.5);#boot data
points(skewness(data)^2, kurtmax - kurtosis(data), pch = 21,cex = 1.5, bg = "darkblue")
points(0, kurtmax - 3, pch = 8, cex = 1.5, lwd = 2); #normal
points(0, kurtmax - 9/5, pch = 2, cex = 1.5, lwd = 2); #uniform
points(2^2, kurtmax - 9, pch = 7, cex = 1.5, lwd = 2); #exponential
points(0, kurtmax - 4.2, pch = 3, cex = 1.5, lwd = 2);# logistic
axis_fun(2,0:ymax,0:ymax,kurtmax-0:ymax)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5);box(lwd=1)
mtext(side=2,line=2,"Kurtosis")
mtext(side=1,line=1.5,"Square of Skewness")
mtext(side=3,"WBID 3240F (Dissolved Oxygen % Sat)",adj=0)

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(-0.25,0.8,legend=c("Observation",paste0("Bootstrapped values\n(N=",boot,")")),
       pch=21,
       lty=0,lwd=c(1,0.01),
       col=c("black","orange"),
       pt.bg=c("darkblue",NA),
       pt.cex=1.5,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0,yjust=0.5)

leg.labs=c("Normal","Uniform","Exponential","Logisitc","Beta","Lognormal","Gamma")#,"(Weibull is close to\nGamma and Lognormal)")
legend(-0.25,0.35,legend=leg.labs,
       pch=c(8,2,7,3,22,NA,NA,NA),
       lty=c(NA,NA,NA,NA,NA,3,2,NA),
       col="black",
       pt.bg=c(rep("black",4),"grey",rep("black",2)),
       pt.cex=1,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0,yjust=0.5,
       title.adj=0,title="Theoretical distributions")
dev.off()

##
dat3240F.xtab$month=as.numeric(format(dat3240F.xtab$date,'%m'))

dat3240F.xtab.mean=ddply(dat3240F.xtab,"date",summarise,mean.DOSat=mean(DOsat.calc,na.rm=T))
dat3240F.xtab.mean$DoY=as.numeric(format(dat3240F.xtab.mean$date,"%j"))
dat3240F.xtab.mean$CY=as.numeric(format(dat3240F.xtab.mean$date,'%Y'))
dat3240F.xtab.mean=subset(dat3240F.xtab.mean,is.na(mean.DOSat)==F)
dat3240F.xtab.mean$dum.val=1:nrow(dat3240F.xtab.mean)

plot(mean.DOSat~date,dat3240F.xtab.mean)
with(dat3240F.xtab.mean,cor.test(mean.DOSat,as.numeric(date),method="kendall"))


dat3240F.xtab.mean.sea=ddply(dat3240F.xtab,c("month","CY"),summarise,mean.DOSat=mean(DOsat.calc,na.rm=T))
tmp=kendallSeasonalTrendTest(mean.DOSat~month+CY,data=dat3240F.xtab.mean.sea)
print(tmp)
##
dep.trend=lm(mean.DOSat~date,dat3240F.xtab.mean)
summary(dep.trend)
gvlma::gvlma(dep.trend)
layout(matrix(1:4,2,2));plot(dep.trend)

lmtest::bgtest(dep.trend,order=3)
# lmtest::dwtest(dep.trend)
shapiro.test(residuals(dep.trend))

dep.trend.mblm=mblm::mblm(mean.DOSat~dum.val,dat3240F.xtab.mean)

ylim.val=c(0,139);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("2006-01-01","2021-01-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# png(filename=paste0(plot.path,"WBID3240F_trend.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));

plot(mean.DOSat~date,dat3240F.xtab.mean,type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(dat3240F.xtab.mean,pt_line(date,mean.DOSat,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5)))
mod.pred=cbind(dat3240F.xtab.mean,predict(dep.trend,dat3240F.xtab.mean,interval="confidence"))
with(mod.pred,shaded.range(date,lwr,upr,"indianred1",lty=0))
with(mod.pred,lines(date,fit,lwd=2,col="indianred1"))
mod.pred2=cbind(dat3240F.xtab.mean,predict(dep.trend.mblm,dat3240F.xtab.mean,interval="confidence"))
with(mod.pred2,shaded.range(date,lwr,upr,"forestgreen",lty=0))
with(mod.pred2,lines(date,fit,lwd=2,col="forestgreen"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Dissolved Oxygen (% Sat)")
mtext(side=1,line=1.75,"Date (Month-Year)")
mtext(side=3,adj=0,"WBID: 3240F")
legend("topleft",legend=c("Daily Data","FDEP Trend (Least Squared) \u00B1 95% CI","Thiel-Sen Trend \u00B1 95% CI"),
       pch=c(21,NA,NA),pt.bg=adjustcolor("dodgerblue1",0.5),col=c("black","indianred1","forestgreen"),
       lty=c(NA,1,1),lwd=c(0.1,1,1),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

# png(filename=paste0(plot.path,"WBID3240F_trendmodel_diag.png"),width=6,height=5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.75));
layout(matrix(1:4,2,2))

cols="grey"#rainbow(length(mod.TP.all$fitted.values))
ylim.val=c(-40,80);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(42,49);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(dep.trend$fitted.values,dep.trend$residuals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
with(dep.trend,points(fitted.values,residuals,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25))
fit.vals=dep.trend$fitted.values
res.vals=dep.trend$residuals
with(lowess(fit.vals,res.vals),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Residuals")
mtext(side=1,line=1.75,"Fitted Values")

ylim.val=c(-3,5);by.y=2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(-3,3);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
rstd=rstandard(dep.trend)
qq.x=qq.function(dep.trend$residuals)
plot(rstd~qq.x,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(qq.x,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
abline(0,1,lty=3)
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2.5,"Standardized Residuals")
mtext(side=1,line=1.75,"Theoretical Quantiles")

ylim.val=c(0,2.5);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(42,49);by.x=1;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
plot(sqrt(abs(rstd))~fit.vals,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
points(fit.vals,sqrt(abs(rstd)),pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(fit.vals,sqrt(abs(rstd))),lines(x,y,col="red",lwd=2))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj);axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,expression(sqrt("Standardized Residuals")))
mtext(side=1,line=1.75,"Fitted Values")

ylim.val=c(-3,3);by.y=1;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,0.01);by.x=0.01;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
lev=hatvalues(dep.trend)
plot(rstd~lev,ylim=ylim.val,xlim=xlim.val,type="n",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey80")
abline(h=0)
points(lev,rstd,pch=21,bg=adjustcolor(cols,0.5),col="grey",lwd=0.1,cex=1.25)
with(lowess(lev,rstd),lines(x,y,col="red",lwd=2))
inf=lm.influence(dep.trend)
hh=seq(min(range(inf)[1],range(inf)[2]/100),xlim.val[2]+(xlim.val[2]*0.5),length.out=101)
hh=hh[hh>0]
crit.cook=0.5
cl.h=sqrt(crit.cook*length(coef(mod.TP.all))*(1-hh)/hh)
lines(hh,cl.h,lty=2,col=2)
lines(hh,-cl.h,lty=2,col=2)
crit.cook=1
cl.h=sqrt(crit.cook*length(coef(mod.TP.all))*(1-hh)/hh)
lines(hh,cl.h,lty=2,col=2)
lines(hh,-cl.h,lty=2,col=2)
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj));axis_fun(2,ymaj,ymin,ymaj);box(lwd=1)
mtext(side=2,line=2,"Standardized Residuals")
mtext(side=1,line=1.75,"Leverage")
dev.off()

acf(residuals(dep.trend))
resid.val=dep.trend$residuals
acf(resid.val)

acf.dep.trend.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.dep.trend.rslt=rbind(acf.dep.trend.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
acf.dep.trend.rslt
points(acf.dep.trend.rslt$lag,acf.dep.trend.rslt$estimate)

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"WBID3240F_trendmodel_ACF.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.dep.trend.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(dep.trend$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.dep.trend.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.dep.trend.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()




samp.screen=ddply(dat3240F.xtab,"station.id",summarise,N.val=N.obs(DOsat.calc))
samp.screen=subset(samp.screen,N.val>20)
dat3240F.trend=subset(dat3240F.xtab,station.id%in%samp.screen$station.id)

site.trend=ddply(dat3240F.trend,"station.id",summarise,N.val=N.obs(DOsat.calc,"NaN"),est=cor.test(DOsat.calc,as.numeric(date),method="kendall")$estimate,
      pval=cor.test(DOsat.calc,as.numeric(date),method="kendall")$p.value)
site.trend$rslt=with(site.trend,ifelse(pval<0.05&est<0,"sig.decline",ifelse(pval<0.05&est>0,"sig.incline","no trend")))
tmp=merge(dat3240F.shp,site.trend,"station.id",all.x=T)

# png(filename=paste0(plot.path,"WBID3240F_trendmap.png"),width=4.25,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1,1,2,3),2,2,byrow=F),widths=c(1))

bbox.lims=bbox(gBuffer(wbid.3240F,width=2000))
bbox.poly=as(raster::extent(gBuffer(wbid.3240F,width=2000)),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

bbox.poly=as(raster::extent(im3.site),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

plotRGB(im2,ext=extent(gBuffer(wbid.3240F,width=2000)))
plot(wbid.3240F,border="grey",col=adjustcolor("white",0.1),lwd=1,add=T)
plot(crop(roads.all,gBuffer(bbox.poly,width=100)),col="grey",lwd=0.8,add=T)
plot(dat3240F.shp,add=T,pch=21,bg="indianred1",cex=1,lwd=0.1)
# plot(bbox.poly,lty=1,lwd=2,border="red",add=T)
plot(gBuffer(as(raster::extent(im2),"SpatialPolygons"),width=-150),add=T,border="black",lwd=1.5)
mapmisc::scaleBar(utm17,"topleft",bty="n",cex=0.8,col="white")

plotRGB(im3.site)
plot(crop(wbid.3240F,gBuffer(bbox.poly,width=-150)),border="grey",col=adjustcolor("white",0.1),lwd=1,add=T,xpd=F)
plot(subset(tmp,is.na(rslt)==T),add=T,pch=21,col="grey",bg="grey",cex=1,lwd=0.1)
plot(subset(tmp,rslt=="sig.decline"),add=T,pch=25,bg="red",cex=1,lwd=0.1)
plot(subset(tmp,rslt=="no trend"),add=T,pch=22,bg="yellow",cex=1,lwd=0.1)
plot(gBuffer(as(raster::extent(im3.site),"SpatialPolygons"),width=-150),add=T,border="black",lwd=1.5)

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend=c("Significant Increasing\n(\u03C4 >0 & \u03C1<0.05)","No Trend(\u03C1 >0.05)","Significant Decreasing\n(\u03C4 <0 & \u03C1<0.05)","Insufficient Data"),
       pch=c(24,22,25,21),
       lty=NA,
       lwd=0.5,
       pt.bg=c("green","yellow","red","grey"),
       col="black",
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj = 0,title="Kendall Trend Result")
dev.off()

### NNC 
## TP=0.12, TN=1.54; not to be exceeded more than once in any three calendar year period
nut.N=ddply(dat3240F.xtab,c("station.id","CY"),summarise,N.TP=N.obs(TP,"NaN"),N.TN=N.obs(TN,"NaN"))
nut.N$TP.scn=with(nut.N,ifelse(N.TP<4,0,1))
nut.N$TN.scn=with(nut.N,ifelse(N.TN<4,0,1))

vars=c("station.id","CY","TP","TN")
vars2=c("station.id","CY","TP.scn","TN.scn")
dat3240F.xtab2=merge(dat3240F.xtab[,vars],nut.N[,vars2],c("station.id","CY"))
dat3240F.AGM=ddply(dat3240F.xtab2,c("station.id","CY"),summarise,
                   TP.AGM=exp(mean(log(ifelse(TP.scn==1,TP,NA)),na.rm=T)),
                   TN.AGM=exp(mean(log(ifelse(TN.scn==1,TN,NA)),na.rm=T)))
dat3240F.AGM$TP.ann.ex=with(dat3240F.AGM,ifelse(TP.AGM>0.12,1,0))
dat3240F.AGM$TP_1_3=with(dat3240F.AGM,ave(TP.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))
dat3240F.AGM$TN.ann.ex=with(dat3240F.AGM,ifelse(TN.AGM>1.54,1,0))
dat3240F.AGM$TN_1_3=with(dat3240F.AGM,ave(TN.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))

subset(dat3240F.AGM,TN_1_3==1|TP_1_3==1)

# WBID3240Q - Popash Creek --------------------------------------------
dat3240Q=read.xlsx(paste0(data.path,"periodOfRecordData3240Q.xlsx"))
dat3240Q$date=date.fun(convertToDate(dat3240Q$date))
dat3240Q=subset(dat3240Q,year%in%seq(2006,2020,1))

dat3240Q.sites=ddply(dat3240Q,c("station.id","lat","long"),summarise,N.val=N.obs(year))
dat3240Q.shp=SpatialPointsDataFrame(dat3240Q.sites[,c("long","lat")],data=dat3240Q.sites[,1:3],proj4string=nad83.pro)
dat3240Q.shp=spTransform(dat3240Q.shp,utm17)
tm_shape(dat3240Q.shp)+tm_dots()

# QA/QC qualifiers 
unique(dat3240Q$rcode)
unique(dat3240Q$xcode)
quals=as.character(unique(dat3240Q$xcode))
spl=strsplit(quals,split="")
quals=data.frame(xcode=quals,
                 q1=sapply(spl,"[",1),
                 q2=sapply(spl,"[",2),
                 q3=sapply(spl,"[",3))
quals$Fatal=with(quals,
                 ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                          q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                          q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,
                        "Y","N"))
dat3240Q=merge(dat3240Q,quals[,c("xcode","Fatal")],"xcode",all.x=T)
dat3240Q=subset(dat3240Q,Fatal=="N")
dat3240Q$HalfMDL=with(dat3240Q,
                      ifelse(is.na(mdl)==T,result,
                             ifelse(result<=mdl,result/2,result)))

params.keep=c("CHLAC","COLOR","COND",
              "DO","DOSAT","TN","TP","TOC","TEMP")
dat3240Q.xtab=reshape2::dcast(subset(dat3240Q,master.code%in%params.keep),wbid+station.id+lat+long+date~master.code,value.var="HalfMDL",mean)
dat3240Q.xtab$Sal=with(dat3240Q.xtab,SalinityCalc(COND,TEMP))
dat3240Q.xtab$DOsat.calc=with(dat3240Q.xtab,DO_PerSat(TEMP,DO,Sal))
dat3240Q.xtab$DoY=as.numeric(format(dat3240Q.xtab$date,"%j"))
dat3240Q.xtab$CY=as.numeric(format(dat3240Q.xtab$date,"%Y"))
dat3240Q.xtab$month=as.numeric(format(dat3240Q.xtab$date,'%m'))

## WQS eval
dat3240Q$time=with(dat3240Q,ifelse(nchar(time)==3,paste0(0,time),time))
dat3240Q$datetime=with(dat3240Q,date.fun(paste(date,time),form="%F %H%M"))
dat3240Q.xtab2=reshape2::dcast(subset(dat3240Q,master.code%in%params.keep),wbid+station.id+lat+long+datetime~master.code,value.var="HalfMDL",mean)
dat3240Q.xtab2$CY=as.numeric(format(dat3240Q.xtab2$date,"%Y"))
dat3240Q.xtab2$Sal=with(dat3240Q.xtab2,SalinityCalc(COND,TEMP))
dat3240Q.xtab2$DOsat.calc=with(dat3240Q.xtab2,DO_PerSat(TEMP,DO,Sal))
dat3240Q.xtab2$DO.TOD.WQS=with(dat3240Q.xtab2,DO.TOD.WQS.stream(datetime))
dat3240Q.xtab2$exceed=with(dat3240Q.xtab2,ifelse(DOsat.calc<DO.TOD.WQS,1,0))

rslt.3240Q=ddply(dat3240Q.xtab2,c("station.id","CY"),summarise,N.exceed=sum(exceed,na.rm=T),N.val=N.obs(DOsat.calc,"NaN"))
rslt.3240Q$PerExceed=with(rslt.3240Q,N.exceed/N.val)*100
rslt.3240Q$status=with(rslt.3240Q,ifelse(PerExceed>38,1,0))
rslt.3240Q
ddply(rslt.3240Q,"station.id",summarise,sum.status=sum(status),n.val=N.obs(status))


## Spatially Average Data
dat3240Q.xtab.mean=ddply(dat3240Q.xtab,"date",summarise,mean.DOSat=mean(DOsat.calc,na.rm=T))
dat3240Q.xtab.mean$DoY=as.numeric(format(dat3240Q.xtab.mean$date,"%j"))
dat3240Q.xtab.mean$CY=as.numeric(format(dat3240Q.xtab.mean$date,'%Y'))
dat3240Q.xtab.mean=subset(dat3240Q.xtab.mean,is.na(mean.DOSat)==F)
dat3240Q.xtab.mean$dum.val=1:nrow(dat3240Q.xtab.mean) #time index

# Kendall Trend
ken.rslt=with(dat3240Q.xtab.mean,cor.test(mean.DOSat,as.numeric(date),method="kendall"))
ken.rslt
# Seasonal Kendall
dat3240Q.xtab.mean.sea=ddply(dat3240Q.xtab,c("month","CY"),summarise, mean.DOSat=mean(DOsat.calc,na.rm=T))
sea.rslt=kendallSeasonalTrendTest(mean.DOSat~month+CY,data=dat3240Q.xtab.mean.sea)
print(sea.rslt)

# Individual station Kendall Trend
samp.screen=ddply(dat3240Q.xtab,"station.id",summarise,N.val=N.obs(DOsat.calc))
samp.screen=subset(samp.screen,N.val>20)
dat3240Q.trend=subset(dat3240Q.xtab,station.id%in%samp.screen$station.id)
dat3240Q.trend$date.num=as.numeric(dat3240Q.trend$date)
site.trend=ddply(dat3240Q.trend,"station.id",summarise,
                 N.val=N.obs(DOsat.calc,"NaN"),
                 est=cor.test(DOsat.calc,date.num,
                              method="kendall")$estimate,
                 pval=cor.test(DOsat.calc,date.num,
                               method="kendall")$p.value)
site.trend$rslt=with(site.trend,ifelse(pval<0.05&est<0,"sig.decline",ifelse(pval<0.05&est>0,"sig.incline","no trend")))
dat3240Q.shp.trend=merge(dat3240Q.shp,site.trend,"station.id",all.x=T)

## Map
wbid.3240Q=subset(wbids,WBID=="3240Q")
###
roi=extent(spTransform(gBuffer(wbid.3240Q,width=5000),wgs84))
im <- cc_location(roi,zoom=13)
im2=projectRaster(im,crs=wkt(utm17))
im2=setValues(im2,scales::rescale(values(im2), c(0,255)))

roi=extent(spTransform(gBuffer(dat3240Q.shp,width=1000),wgs84))
im <- cc_location(roi)
im2.site=projectRaster(im,crs=wkt(utm17))
im2.site=setValues(im2.site,scales::rescale(values(im2.site), c(0,255)))

# png(filename=paste0(plot.path,"WBID3240Q_trendmap.png"),width=3.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1,1,2,3),2,2,byrow=F),widths=c(1))

bbox.lims=bbox(gBuffer(wbid.3240Q,width=2000))
bbox.poly=as(raster::extent(gBuffer(wbid.3240Q,width=2000)),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

bbox.poly=as(raster::extent(im2.site),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

plotRGB(im2,ext=extent(gBuffer(wbid.3240F,width=2000)))
plot(wbid.3240Q,border="grey",col=adjustcolor("white",0.1),lwd=1,add=T)
plot(crop(roads.all,gBuffer(bbox.poly,width=100)),col="grey",lwd=0.8,add=T)
plot(dat3240Q.shp,add=T,pch=21,bg="indianred1",cex=1,lwd=0.1)
# plot(gBuffer(as(raster::extent(im2),"SpatialPolygons"),width=-150),add=T,border="black",lwd=1.5)
mapmisc::scaleBar(utm17,"topleft",bty="n",cex=0.8,col="white")

plotRGB(im2,ext=extent(gBuffer(dat3240Q.shp.trend,width=2000)))
plot(wbid.3240Q,border="grey",lwd=1,add=T,xpd=F)
plot(subset(dat3240Q.shp.trend,is.na(rslt)==T),add=T,pch=21,col="grey",bg="grey",cex=1,lwd=0.1)
plot(subset(dat3240Q.shp.trend,rslt=="sig.decline"),add=T,pch=25,bg="red",cex=1,lwd=0.1)
plot(subset(dat3240Q.shp.trend,rslt=="no trend"),add=T,pch=22,bg="yellow",cex=1,lwd=0.1)
# plot(gBuffer(as(raster::extent(im3.site),"SpatialPolygons"),width=-150),add=T,border="black",lwd=1.5)

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend=c("Significant Increasing\n(\u03C4 >0 & \u03C1<0.05)","No Trend(\u03C1 >0.05)","Significant Decreasing\n(\u03C4 <0 & \u03C1<0.05)","Insufficient Data"),
       pch=c(24,22,25,21),
       lty=NA,
       lwd=0.5,
       pt.bg=c("green","yellow","red","grey"),
       col="black",
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj = 0,title="Kendall Trend Result")
dev.off()


# lm trend
dep.trend=lm(mean.DOSat~date,dat3240Q.xtab.mean)
summary(dep.trend)
gvlma::gvlma(dep.trend)
layout(matrix(1:4,2,2));plot(dep.trend)

lmtest::bgtest(dep.trend,order=3)
shapiro.test(residuals(dep.trend))

dep.trend.mblm=mblm::mblm(mean.DOSat~dum.val,dat3240Q.xtab.mean)

ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("2006-01-01","2021-01-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# png(filename=paste0(plot.path,"WBID3240Q_trend.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));

plot(mean.DOSat~date,dat3240Q.xtab.mean,type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(dat3240Q.xtab.mean,pt_line(date,mean.DOSat,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5)))
mod.pred=cbind(dat3240Q.xtab.mean,predict(dep.trend,dat3240Q.xtab.mean,interval="confidence"))
with(mod.pred,shaded.range(date,lwr,upr,"indianred1",lty=0))
with(mod.pred,lines(date,fit,lwd=2,col="indianred1"))
mod.pred2=cbind(dat3240Q.xtab.mean,predict(dep.trend.mblm,dat3240Q.xtab.mean,interval="confidence"))
with(mod.pred2,shaded.range(date,lwr,upr,"forestgreen",lty=0))
with(mod.pred2,lines(date,fit,lwd=2,col="forestgreen"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Dissolved Oxygen (% Sat)")
mtext(side=1,line=1.75,"Date (Month-Year)")
mtext(side=3,adj=0,"WBID: 3240Q")
legend("topleft",legend=c("Daily Data","FDEP Trend (Least Squared) \u00B1 95% CI","Thiel-Sen Trend \u00B1 95% CI"),
       pch=c(21,NA,NA),pt.bg=adjustcolor("dodgerblue1",0.5),col=c("black","indianred1","forestgreen"),
       lty=c(NA,1,1),lwd=c(0.1,1,1),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

acf(residuals(dep.trend))
resid.val=dep.trend$residuals
acf(resid.val)

acf.dep.trend.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.dep.trend.rslt=rbind(acf.dep.trend.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
acf.dep.trend.rslt
points(acf.dep.trend.rslt$lag,acf.dep.trend.rslt$estimate)

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"WBID3240Q_trendmodel_ACF.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.dep.trend.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(dep.trend$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.dep.trend.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.dep.trend.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

## Data distribution
tmp=subset(dat3240Q.xtab,is.na(DOsat.calc)==F)$DOsat.calc
shapiro.test(tmp)
nortest::ad.test(tmp)

plotdist(tmp,"norm",para=list(mean=mean(tmp),sd=sd(tmp)))

# png(filename=paste0(plot.path,"WBID3240Q_normDist.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:4,2,2,byrow=T))
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.5));

#hist
h=hist(tmp, plot=F)
xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
yhist <- dnorm(xhist,mean(tmp),sd(tmp))
xlim.val=range(h$breaks);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,max(yhist)+max(yhist)*0.1);by.y=0.005;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(tmp, freq = FALSE,ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,yaxs="i")
lines(xhist,yhist, lty = 1,col="red",lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"Denisty")
mtext(side=1,line=1.75,"DO (% Sat)")

#QQplot
x.val=qnorm(ppoints(sort(tmp)),mean(tmp),sd(tmp))
y.val=sort(tmp)
xlim.val=c(-10,100);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Quantile")
mtext(side=1,line=1.75,"Theoretical Quantile")

#CDFplot
x.val=y.val
y.val=ppoints(x.val)
xlim.val=c(0,120);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
sfin=seq(min(h$breaks),max(h$breaks),length.out=100)
lines(sfin,pnorm(sfin,mean(tmp),sd(tmp)),col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"CDF")
mtext(side=1,line=1.75,"DO (%Sat)")

#PPplot
x.val=pnorm(sort(tmp),mean(tmp),sd(tmp))
y.val=ppoints(sort(tmp))
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Probabilities")
mtext(side=1,line=1.75,"Theoretical Probabilities")
dev.off()


### NNC 
## TP=0.12, TN=1.54; not to be exceeded more than once in any three calendar year period
nut.N=ddply(dat3240Q.xtab,c("station.id","CY"),summarise,N.TP=N.obs(TP,"NaN"),N.TN=N.obs(TN,"NaN"))
nut.N$TP.scn=with(nut.N,ifelse(N.TP<4,0,1))
nut.N$TN.scn=with(nut.N,ifelse(N.TN<4,0,1))

vars=c("station.id","CY","TP","TN")
vars2=c("station.id","CY","TP.scn","TN.scn")
dat3240Q.xtab2=merge(dat3240Q.xtab[,vars],nut.N[,vars2],c("station.id","CY"))
dat3240Q.AGM=ddply(dat3240Q.xtab2,c("station.id","CY"),summarise,
                   TP.AGM=exp(mean(log(ifelse(TP.scn==1,TP,NA)),na.rm=T)),
                   TN.AGM=exp(mean(log(ifelse(TN.scn==1,TN,NA)),na.rm=T)))
dat3240Q.AGM$TP.ann.ex=with(dat3240Q.AGM,ifelse(TP.AGM>0.12,1,0))
dat3240Q.AGM$TP_1_3=with(dat3240Q.AGM,ave(TP.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))
dat3240Q.AGM$TN.ann.ex=with(dat3240Q.AGM,ifelse(TN.AGM>1.54,1,0))
dat3240Q.AGM$TN_1_3=with(dat3240Q.AGM,ave(TN.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))

subset(dat3240Q.AGM,TN_1_3==1|TP_1_3==1)

# WBID3235C - Cypress Creek --------------------------------------------
dat3235C=read.xlsx(paste0(data.path,"periodOfRecordData3235C.xlsx"))
dat3235C$date=date.fun(convertToDate(dat3235C$date))
dat3235C=subset(dat3235C,year%in%seq(2006,2020,1))

cypresscreeksites=c("21FLBABRCYPRESS_HEAD","21FLBABRCYPRESS_OUTFLOW","21FLFTM 28020237","21FLGW  56335","21FLFTM G3SD0084",'21FLEECOCYPRESSGR',"21FLFTM CYPRESSGR")

dat3235C.sites=ddply(dat3235C,c("station.id","lat","long"),summarise,N.val=N.obs(year))
dat3235C.shp=SpatialPointsDataFrame(dat3235C.sites[,c("long","lat")],data=dat3235C.sites[,1:3],proj4string=nad83.pro)
dat3235C.shp=spTransform(dat3235C.shp,utm17)
tm_shape(dat3235C.shp)+tm_dots()+
  tm_shape(subset(dat3235C.shp,!(station.id%in%cypresscreeksites)))+tm_dots("red")
tm_shape(subset(dat3235C.shp,station.id%in%c("21FLEECOCYPRESSGR","21FLWWCO23-5GR")))+tm_dots()

tm_shape(dat3235C.shp)+tm_dots()+
  tm_shape(subset(dat3235C.shp,station.id=="21FLFTM 28020237"))+tm_dots("red")

# dat3235C=subset(dat3235C,station.id%in%cypresscreeksites)
# QA/QC qualifiers 
unique(dat3235C$rcode)
unique(dat3235C$xcode)
quals=as.character(unique(dat3235C$xcode))
spl=strsplit(quals,split="")
quals=data.frame(xcode=quals,
                 q1=sapply(spl,"[",1),
                 q2=sapply(spl,"[",2),
                 q3=sapply(spl,"[",3))
quals$Fatal=with(quals,
                 ifelse(q1%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                          q2%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER|
                          q3%in%subset(dat.qual,FATALYN=="Y")$QUALIFIER,
                        "Y","N"))
dat3235C=merge(dat3235C,quals[,c("xcode","Fatal")],"xcode",all.x=T)
dat3235C=subset(dat3235C,Fatal=="N")
dat3235C$HalfMDL=with(dat3235C,
                      ifelse(is.na(mdl)==T,result,
                             ifelse(result<=mdl,result/2,result)))

params.keep=c("CHLAC","COLOR","COND",
              "DO","DOSAT","TN","TP","TOC","TEMP")
dat3235C.xtab=reshape2::dcast(subset(dat3235C,master.code%in%params.keep),wbid+station.id+lat+long+date~master.code,value.var="HalfMDL",mean)
dat3235C.xtab$Sal=with(dat3235C.xtab,SalinityCalc(COND,TEMP))
dat3235C.xtab$DOsat.calc=with(dat3235C.xtab,DO_PerSat(TEMP,DO,Sal))
dat3235C.xtab$DoY=as.numeric(format(dat3235C.xtab$date,"%j"))
dat3235C.xtab$CY=as.numeric(format(dat3235C.xtab$date,"%Y"))
dat3235C.xtab$month=as.numeric(format(dat3235C.xtab$date,'%m'))

## WQS eval
dat3235C$time=with(dat3235C,ifelse(nchar(time)==3,paste0(0,time),time))
dat3235C$datetime=with(dat3235C,date.fun(paste(date,time),form="%F %H%M"))

dat3235C.xtab2=reshape2::dcast(subset(dat3235C,master.code%in%params.keep),wbid+station.id+lat+long+datetime~master.code,value.var="HalfMDL",mean)
dat3235C.xtab2$CY=as.numeric(format(dat3235C.xtab2$date,"%Y"))
dat3235C.xtab2$month=as.numeric(format(dat3235C.xtab2$date,"%m"))
dat3235C.xtab2$Sal=with(dat3235C.xtab2,SalinityCalc(COND,TEMP))
dat3235C.xtab2$DOsat.calc=with(dat3235C.xtab2,DO_PerSat(TEMP,DO,Sal))
dat3235C.xtab2$DO.TOD.WQS=with(dat3235C.xtab2,DO.TOD.WQS.stream(datetime))
dat3235C.xtab2$exceed=with(dat3235C.xtab2,ifelse(DOsat.calc<DO.TOD.WQS,1,0))

rslt.3235C=ddply(dat3235C.xtab2,c("station.id","CY"),summarise,N.exceed=sum(exceed,na.rm=T),N.val=N.obs(DOsat.calc,"NaN"))
rslt.3235C$PerExceed=with(rslt.3235C,N.exceed/N.val)*100
rslt.3235C$status=with(rslt.3235C,ifelse(PerExceed>38,1,0))
rslt.3235C
ddply(rslt.3235C,"station.id",summarise,sum.status=sum(status),n.val=N.obs(status))

kendallSeasonalTrendTest(DOsat.calc~month+CY,data=subset(dat3235C.xtab2,station.id=="21FLEECOCYPRESSGR"))

## Spatially Average Data
dat3235C.xtab.mean=ddply(dat3235C.xtab,"date",summarise,mean.DOSat=mean(DOsat.calc,na.rm=T))
dat3235C.xtab.mean$DoY=as.numeric(format(dat3235C.xtab.mean$date,"%j"))
dat3235C.xtab.mean$CY=as.numeric(format(dat3235C.xtab.mean$date,'%Y'))
dat3235C.xtab.mean=subset(dat3235C.xtab.mean,is.na(mean.DOSat)==F)
dat3235C.xtab.mean$dum.val=1:nrow(dat3235C.xtab.mean) #time index

# Kendall Trend
ken.rslt=with(dat3235C.xtab.mean,cor.test(mean.DOSat,as.numeric(date),method="kendall"))
ken.rslt
# Seasonal Kendall
dat3235C.xtab.mean.sea=ddply(dat3235C.xtab,c("month","CY"),summarise, mean.DOSat=mean(DOsat.calc,na.rm=T))
sea.rslt=kendallSeasonalTrendTest(mean.DOSat~month+CY,data=dat3235C.xtab.mean.sea)
print(sea.rslt)

# Individual station Kendall Trend
samp.screen=ddply(dat3235C.xtab,"station.id",summarise,N.val=N.obs(DOsat.calc))
# samp.screen=subset(samp.screen,N.val>20)
dat3235C.trend=subset(dat3235C.xtab,station.id%in%subset(samp.screen,N.val>20)$station.id)
dat3235C.trend$date.num=as.numeric(dat3235C.trend$date)
site.trend=ddply(dat3235C.trend,"station.id",summarise,
                 N.val=N.obs(DOsat.calc,"NaN"),
                 est=cor.test(DOsat.calc,date.num,
                              method="kendall")$estimate,
                 pval=cor.test(DOsat.calc,date.num,
                               method="kendall")$p.value)
site.trend$rslt=with(site.trend,ifelse(pval<0.05&est<0,"sig.decline",ifelse(pval<0.05&est>0,"sig.incline","no trend")))
dat3235C.shp.trend=merge(dat3235C.shp,site.trend,"station.id",all.x=T)

## Map
wbid.3235C=subset(wbids,WBID=="3235C")
###
roi=extent(spTransform(gBuffer(wbid.3235C,width=5000),wgs84))
im <- cc_location(roi,zoom=13)
im2=projectRaster(im,crs=wkt(utm17))
im2=setValues(im2,scales::rescale(values(im2), c(0,255)))

roi=extent(spTransform(gBuffer(dat3235C.shp,width=1000),wgs84))
im <- cc_location(roi)
im2.site=projectRaster(im,crs=wkt(utm17))
im2.site=setValues(im2.site,scales::rescale(values(im2.site), c(0,255)))

# png(filename=paste0(plot.path,"WBID3235C_trendmap.png"),width=5,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.1,0.1,0.1,0.1),oma=c(0.5,0.5,0.5,0.5));
layout(matrix(c(1:2),1,2,byrow=F),widths=c(1))

bbox.lims=bbox(gBuffer(wbid.3235C,width=2000))
bbox.poly=as(raster::extent(gBuffer(wbid.3235C,width=2000)),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

bbox.poly=as(raster::extent(im2.site),"SpatialPolygons")#makes the polygon
proj4string(bbox.poly)=utm17

plotRGB(im2,ext=extent(gBuffer(wbid.3235C,width=2000)))
plot(wbid.3235C,border="grey",col=adjustcolor("white",0.1),lwd=1,add=T)
plot(crop(roads.all,gBuffer(bbox.poly,width=100)),col="grey",lwd=0.8,add=T)
#plot(dat3235C.shp,add=T,pch=21,bg="indianred1",cex=1,lwd=0.1)
# plot(gBuffer(as(raster::extent(im2),"SpatialPolygons"),width=-150),add=T,border="black",lwd=1.5)
plot(subset(dat3235C.shp.trend,is.na(rslt)==T),add=T,pch=21,col="grey",bg="grey",cex=1,lwd=0.1)
plot(subset(dat3235C.shp.trend,rslt=="sig.decline"),add=T,pch=25,bg="red",cex=1,lwd=0.1)
plot(subset(dat3235C.shp.trend,rslt=="no trend"),add=T,pch=22,bg="yellow",cex=1,lwd=0.1)
mapmisc::scaleBar(utm17,"topleft",bty="n",cex=0.8,col="white")

plot(0:1,0:1,type="n",axes=F,ylab=NA,xlab=NA)
legend(0.5,0.5,legend=c("Significant Increasing\n(\u03C4 >0 & \u03C1<0.05)","No Trend(\u03C1 >0.05)","Significant Decreasing\n(\u03C4 <0 & \u03C1<0.05)","Insufficient Data"),
       pch=c(24,22,25,21),
       lty=NA,
       lwd=0.5,
       pt.bg=c("green","yellow","red","grey"),
       col="black",
       pt.cex=1.5,ncol=1,cex=0.9,bty="n",y.intersp=1.5,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5,
       title.adj = 0,title="Kendall Trend Result")
dev.off()


# lm trend
dep.trend=lm(mean.DOSat~date,dat3235C.xtab.mean)
summary(dep.trend)
gvlma::gvlma(dep.trend)
layout(matrix(1:4,2,2));plot(dep.trend)

lmtest::bgtest(dep.trend,order=3)
shapiro.test(residuals(dep.trend))

dep.trend.mblm=mblm::mblm(mean.DOSat~dum.val,dat3235C.xtab.mean)

ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=date.fun(c("2006-01-01","2021-01-01"));xmaj=seq(xlim.val[1],xlim.val[2],"3 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")
# png(filename=paste0(plot.path,"WBID3235C_trend.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,2,0.75,1),oma=c(2,2,0.25,0.5));

plot(mean.DOSat~date,dat3235C.xtab.mean,type="n",ylim=ylim.val,xlim=xlim.val,axes=F,ann=F)
abline(h=ymaj,v=xmaj,lty=1,col=adjustcolor("grey",0.5))
with(dat3235C.xtab.mean,pt_line(date,mean.DOSat,2,adjustcolor("dodgerblue1",0.5),1,21,adjustcolor("dodgerblue1",0.5)))
mod.pred=cbind(dat3235C.xtab.mean,predict(dep.trend,dat3235C.xtab.mean,interval="confidence"))
with(mod.pred,shaded.range(date,lwr,upr,"indianred1",lty=0))
with(mod.pred,lines(date,fit,lwd=2,col="indianred1"))
mod.pred2=cbind(dat3235C.xtab.mean,predict(dep.trend.mblm,dat3235C.xtab.mean,interval="confidence"))
with(mod.pred2,shaded.range(date,lwr,upr,"forestgreen",lty=0))
with(mod.pred2,lines(date,fit,lwd=2,col="forestgreen"))
axis_fun(1,line=-0.5,xmaj,xmin,format(xmaj,"%m-%Y"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Dissolved Oxygen (% Sat)")
mtext(side=1,line=1.75,"Date (Month-Year)")
mtext(side=3,adj=0,"WBID: 3235C")
legend("topleft",legend=c("Daily Data","FDEP Trend (Least Squared) \u00B1 95% CI","Thiel-Sen Trend \u00B1 95% CI"),
       pch=c(21,NA,NA),pt.bg=adjustcolor("dodgerblue1",0.5),col=c("black","indianred1","forestgreen"),
       lty=c(NA,1,1),lwd=c(0.1,1,1),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

acf(residuals(dep.trend))
resid.val=dep.trend$residuals
acf(resid.val)

acf.dep.trend.rslt=data.frame()
for(h in 0:24){
  #demean
  #x=sweep(as.matrix(wq.dat.xtab.mon$mean.TP),2,colMeans(as.matrix(wq.dat.xtab.mon$mean.TP),na.rm=T))
  lagged=lag(as.zoo(resid.val),-h,na.pad=T)
  tmp.dat=as.zoo(resid.val)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  acf.dep.trend.rslt=rbind(acf.dep.trend.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
acf.dep.trend.rslt
points(acf.dep.trend.rslt$lag,acf.dep.trend.rslt$estimate)

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,24);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"WBID3235C_trendmodel_ACF.png"),width=4,height=3,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,acf.dep.trend.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(length(dep.trend$residuals))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-1,25,25,-1),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(acf.dep.trend.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(acf.dep.trend.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("ACF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()

## Data distribution
tmp=subset(dat3235C.xtab,is.na(DOsat.calc)==F)$DOsat.calc
shapiro.test(tmp)
nortest::ad.test(tmp)

plotdist(tmp,"norm",para=list(mean=mean(tmp),sd=sd(tmp)))

# png(filename=paste0(plot.path,"WBID3235C_normDist.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:4,2,2,byrow=T))
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.5));

#hist
h=hist(tmp, plot=F)
xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
yhist <- dnorm(xhist,mean(tmp),sd(tmp))
xlim.val=range(h$breaks);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,max(yhist)+max(yhist)*0.1);by.y=0.005;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(tmp, freq = FALSE,ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,yaxs="i")
lines(xhist,yhist, lty = 1,col="red",lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"Denisty")
mtext(side=1,line=1.75,"DO (% Sat)")

#QQplot
x.val=qnorm(ppoints(sort(tmp)),mean(tmp),sd(tmp))
y.val=sort(tmp)
xlim.val=c(-10,100);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Quantile")
mtext(side=1,line=1.75,"Theoretical Quantile")

#CDFplot
x.val=y.val
y.val=ppoints(x.val)
xlim.val=c(0,120);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
sfin=seq(min(h$breaks),max(h$breaks),length.out=100)
lines(sfin,pnorm(sfin,mean(tmp),sd(tmp)),col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"CDF")
mtext(side=1,line=1.75,"DO (%Sat)")

#PPplot
x.val=pnorm(sort(tmp),mean(tmp),sd(tmp))
y.val=ppoints(sort(tmp))
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Probabilities")
mtext(side=1,line=1.75,"Theoretical Probabilities")
dev.off()

## Cypress Creek only
## Data distribution
tmp=subset(dat3235C.xtab,station.id%in%cypresscreeksites&is.na(DOsat.calc)==F)$DOsat.calc
shapiro.test(tmp)
nortest::ad.test(tmp)

# png(filename=paste0(plot.path,"WBID3235C_normDist_CCsitesOnly.png"),width=6.5,height=5,units="in",res=200,type="windows",bg="white")
layout(matrix(1:4,2,2,byrow=T))
par(family="serif",mar=c(2,3,1,0.5),oma=c(2,1,0.25,0.5));

#hist
h=hist(tmp, plot=F)
xhist <- seq(min(h$breaks), max(h$breaks), length = 1000)
yhist <- dnorm(xhist,mean(tmp),sd(tmp))
xlim.val=range(h$breaks);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,max(yhist)+max(yhist)*0.15);by.y=0.005;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
hist(tmp, freq = FALSE,ylim=ylim.val,xlim=xlim.val,axes=F,ann=F,yaxs="i")
lines(xhist,yhist, lty = 1,col="red",lwd=2)
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"Denisty")
mtext(side=1,line=1.75,"DO (% Sat)")

#QQplot
x.val=qnorm(ppoints(sort(tmp)),mean(tmp),sd(tmp))
y.val=sort(tmp)
xlim.val=c(-10,100);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,120);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Quantile")
mtext(side=1,line=1.75,"Theoretical Quantile")

#CDFplot
x.val=y.val
y.val=ppoints(x.val)
xlim.val=c(0,120);by.x=20;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
sfin=seq(min(h$breaks),max(h$breaks),length.out=100)
lines(sfin,pnorm(sfin,mean(tmp),sd(tmp)),col="red")
axis_fun(1,xmaj,xmin,xmaj,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=3,"CDF")
mtext(side=1,line=1.75,"DO (%Sat)")

#PPplot
x.val=pnorm(sort(tmp),mean(tmp),sd(tmp))
y.val=ppoints(sort(tmp))
xlim.val=c(0,1);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
plot(x.val,y.val,xlim=xlim.val,ylim=ylim.val,axes=F,ann=F,type="n")
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
points(x.val,y.val,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=1.25)
abline(0,1,col="red")
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Empirical Probabilities")
mtext(side=1,line=1.75,"Theoretical Probabilities")
dev.off()

### NNC 
## TP=0.12, TN=1.54; not to be exceeded more than once in any three calendar year period
nut.N=ddply(dat3235C.xtab,c("station.id","CY"),summarise,N.TP=N.obs(TP,"NaN"),N.TN=N.obs(TN,"NaN"))
nut.N$TP.scn=with(nut.N,ifelse(N.TP<4,0,1))
nut.N$TN.scn=with(nut.N,ifelse(N.TN<4,0,1))

vars=c("station.id","CY","TP","TN")
vars2=c("station.id","CY","TP.scn","TN.scn")
dat3235C.xtab2=merge(dat3235C.xtab[,vars],nut.N[,vars2],c("station.id","CY"))
dat3235C.AGM=ddply(dat3235C.xtab2,c("station.id","CY"),summarise,
                   TP.AGM=exp(mean(log(ifelse(TP.scn==1,TP,NA)),na.rm=T)),
                   TN.AGM=exp(mean(log(ifelse(TN.scn==1,TN,NA)),na.rm=T)))
dat3235C.AGM$TP.ann.ex=with(dat3235C.AGM,ifelse(TP.AGM>0.12,1,0))
dat3235C.AGM$TP_1_3=with(dat3235C.AGM,ave(TP.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))
dat3235C.AGM$TN.ann.ex=with(dat3235C.AGM,ifelse(TN.AGM>1.54,1,0))
dat3235C.AGM$TN_1_3=with(dat3235C.AGM,ave(TN.ann.ex,station.id,FUN=function(x) ifelse(c(rep(NA,2),rollsum(x,k=3))>=3,1,0)))

subset(dat3235C.AGM,TN_1_3==1|TP_1_3==1)
