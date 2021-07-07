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
  
  ncomp<-25
  #Store original network
  network<-ergm$network
  #Store component distribution of original network
  obs<-t(component.dist(network)$cdist)
  #Count number of nodes
  nodes<-length(obs)
  #Combine component sizes past ncomp
  obs[ncomp]<-sum(obs[ncomp:nodes])
  #Store component distribution of original network as wide data
  wide<-as.data.frame(obs)
  #Simulate and store specified number of networks
  sims<-simulate(ergm,nsim)
  #Loop through all network simulations
  for(i in 1:nsim){
    #Store component distribution for simulation i
    wide[i+1,]<-component.dist(sims[[i]])$cdist
    #Combine largest component size considered with all others
    wide[i+1,ncomp]<-sum(wide[i,ncomp:nodes])
  }
  #Exclude all component sizes larger than ncomp
  wide<-wide[,1:ncomp]
  #Store network number
  wide$network<-c(0:nsim)
  #Store component sizes
  sizes<-c(as.character(1:(ncomp-1)),paste(ncomp,"+",sep=""))
  #Prepopulate vector of mean number of components by size
  means<-matrix(0,ncomp,1)
  #Loop through all component sizes considered
  for(i in 1:ncomp){
    means[i]<-mean(wide[2:nsim,i])
  }
  #Store component sizes
  summary<-as.data.frame(sizes)
  #Store true component distribution
  summary$obs<-obs[1:ncomp]
  #Store means for component distribution
  summary$mean<-round(means,2)
  #Print summary
  print.data.frame(summary,row.names=FALSE)
  #Convert data to long format
  long<-gather(wide,size,number,V1:paste("V",ncomp,sep=""))
  #Convert component size from character
  long$size<-rep(1:ncomp,each=nsim+1)
  #Convert component size to factor
  long$size<-factor(long$size)
  #Relabel component sizes
  levels(long$size)<-sizes
  ggplot(long,aes(x=size,y=number))+
    geom_boxplot(data=long[long$network!=0,])+
    geom_point(data=long[long$network==0,])+
    geom_line(data=long[long$network==0,],aes(group=network))+
    labs(x="Component Size",y="Number of Components")
}

#Test
#ComponentGOF(ergm,10)
