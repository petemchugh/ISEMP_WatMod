#### Code for summarizing LCM outputs
#library...
rm(list = ls(all=TRUE)) #Clean up shop before beginning
library(fanplot)
library(reshape)
library(tidyr)
library(colorRamps)
library(RColorBrewer)

this_is_the_place<-"C:\\Users\\Peter\\Desktop" #file.path(choose.dir())
sim.dat.base<-read.table(file.choose(),header=TRUE,sep=",") #opt

#Name for scenarios
base<-"Testing"

#Subset to the time horizon and chuck the burn-in
burnin<-10
years<-40
QET<-10

simsub.base<-subset(sim.dat.base,yr>burnin)
simsub.base<-subset(simsub.base,yr<(burnin+years+1))

simsub.3<-subset(simsub.3,yr<(burnin+years+1))

simsub.base$RBTfry.pct<-simsub.base$RBTfry/(simsub.base$RBTfry+simsub.base$Nfry) #convert to 1:50


#attach(simsub1)
simsub.base$yr<-simsub.base$yr-burnin #convert to 1:50


finalN.base<-simsub.base$Nsp[simsub.base$yr==years]
#hist(finalN.base,br=seq(0,2000,50))
QEprob.base<-length(finalN.base[finalN.base<QET])/length(finalN.base)
mdNbase<-median(finalN.base)
plot(density(finalN.base),xlim=c(-100,800))


help(hist)
base50<-median(finalN.base)
baseMN<-mean(finalN.base)
base25<-quantile(finalN.base,0.25)
base75<-quantile(finalN.base,0.75)




##########################################################################################################################
##########################################################################################################################
##########################################################################################################################
# EVALUATION OF QUASI-EXTINCTION THRESHOLDS AND A FUNCTION FOR EXTRACTING THEM FROM RUNS...
##########################################################################################################################
QET<-10




##################################################
# For tabulating spawner abundance through time
##################################################
yr1<-1:years

meanN.base<-rep(NA,years)
medianN.base<-yr1
tenth.base<-rep(NA,years)
ninetieth.base<-rep(NA,years)
twentyfifth.base<-rep(NA,years)
seventyfifth.base<-rep(NA,years)



y<-1
for(y in 1:years)
{
  NNN.base<-simsub.base$Nsp[simsub.base$yr==y]
  meanN.base[y]<-mean(NNN.base)
  medianN.base[y]<-median(NNN.base)
  tenth.base[y]<-quantile(NNN.base,0.10)
  ninetieth.base[y]<-quantile(NNN.base,0.90)
  twentyfifth.base[y]<-quantile(NNN.base,0.25)
  seventyfifth.base[y]<-quantile(NNN.base,0.75)
  

}
#######################################################





###############################################################################
###############################################################################
### Figures of time series with polygons for middle 50% of distn

# Base
dfFan<-spread(simsub.base[c(1,2,23)],yr,Nsp)
plot(simsub.base$yr,simsub.base$Nsp,xaxt="n",pch="",cex=1,cex.lab=1.5,
     xlab = "Simulation Year",ylab = "Spawners",bg="dodgerblue",
     ylim=c(0,3000),xlim=c(2,39))
pal<-colorRampPalette(rev(brewer.pal(5,"Blues")))
fan(dfFan[2:41],fan.col=pal,ln=c(0.1,0.9),
    llab=FALSE,ln.col="black",rlab=FALSE)
lines(medianN.base,col="black", lwd=6)
#lines(meanN.base,col="black", lwd=6)
#abline(h=QET,lwd=1.5,lty=2)
abline(h=1774,lty=2,lwd=2)
#text(41,950,paste("Base: ",base,sep=""),pos=2,cex=1.1)
axis(1,at=seq(1,length(rep(NA,years)),5), labels = seq(1,length(rep(NA,years)),5)-1)
box(lwd=2)  
