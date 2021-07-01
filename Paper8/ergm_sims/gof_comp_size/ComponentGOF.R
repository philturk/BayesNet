require(tidyr)
require(ggplot2)
require(ergm)
require(sna)
source("Network.R")

#ERGM for testing
#ergm<-ergm(GeneticNetwork~edges)

#Function for fit diagnostics for component distribution
ComponentGOF<-function(ergm,nsim=100){
  #ergm: output from ergm()
  #nsim: number of networks to simulate (default of 100)
  
  ncomp<-9
  #Store original network
  network<-ergm$network
  #Store component distribution of original network
  true<-t(component.dist(network)$cdist)
  #Count number of nodes
  nodes<-length(true)
  #Store component distribution of original network as wide data
  wide<-as.data.frame(true)
  #Simulate and store specified number of networks
  sims<-simulate(ergm,nsim)
  #Loop through all network simulations
  for(i in 1:nsim){
    #Store component distribution for simulation i
    wide[i+1,]<-component.dist(sims[[i]])$cdist
    #Combine largest component size considered with all others
    wide[i,ncomp]<-sum(wide[i,ncomp:nodes])
  }
  #Exclude all component sizes largere than ncomp
  wide<-wide[,1:ncomp]
  #Store network number
  wide$network<-c(0:nsim)
  #Convert data to long format
  long<-gather(wide,size,number,V1:paste("V",ncomp,sep=""))
  #Convert component size to factor
  long$size<-factor(long$size)
  #Relable component sizes
  levels(long$size)<-c("1","2","3","4","5","6","7","8","9+")
  ggplot(long,aes(x=size,y=number))+
    geom_boxplot(data=long[long$network!=0,])+
    geom_point(data=long[long$network==0,])+
    geom_line(data=long[long$network==0,],aes(group=network))+
    labs(x="Component Size",y="Number of Components")
}

#ComponentGOF(ergm)
