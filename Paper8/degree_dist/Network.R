#IMPORT & FORMAT DATA
file_location = dirname(rstudioapi::getActiveDocumentContext()$path)
# print("Running Network.R in file location ", as.character(file_location))

#import genetic distances CSV
#ID1, ID2, Distance
#ID1 is first node (j), ID2 is second node (k)
# Distances <- read.csv("distances_0.02.csv",
                      # stringsAsFactors=FALSE)

Distances <- read.csv(paste0(file_location, "\\distances_0.02.csv"),
                      stringsAsFactors=FALSE)

#import ID key CSV
#ID, Sequences
#ID is either ###v######## or ###-###
IDs <- read.csv(paste0(file_location, "\\ids.csv"))
#1259 unique IDs
#length(unique(IDs$ID))

#format ID1 to only have ID number, without extra text information
#9 sections total for each observation, first section is ID
Distances$ID1 <- matrix(unlist(strsplit(Distances$ID1,"|",fixed=TRUE)),
                        nrow=5558,ncol=9,byrow=TRUE)[,1]
#1258 unique ID1s
#length(unique(Distances$ID1))

#format ID2,  to only have ID number, without extra text information
#5558 rows, 9 sections total for each row, 1st section is ID
#last 3 rows are HXB2_prrt|01011983|lab-strain, only 3 sections
#add on filler NAs for missing sections
#replace last 2 rows of NA with "HXB2_prrt", correctly leaves last 3 ID2s as "HXB2_prrt"
Distances$ID2 <- matrix(c(unlist(strsplit(Distances$ID2,"|",fixed=TRUE)),rep(NA,18)),
                        nrow=5558,ncol=9,byrow=TRUE)[,1]
Distances$ID2[c(5557,5558)] <- c("HXB2_prrt","HXB2_prrt")
#1259 unique ID2s
#length(unique(Distances$ID2))


#remove any pairs that link to themselves
#replace distance of pair with NA if pair is to itsel
for (i in 1:length(Distances$ID1)) {
  Distances[i,3] <- ifelse(Distances$ID1[i]==Distances$ID2[i],NA,Distances$Distance[i])
}
#remove rows of data that have NAs (= links to self)
Distances <- Distances[complete.cases(Distances),]

#remove duplicate pairs
library(dplyr)
Distances <- Distances %>% distinct(ID1,ID2,.keep_all=TRUE)

#create EdgeList = links
#filter Distances to create 2-column matrix, where ID1 & ID2 appear on list if Distance <0.15
#remove Distance variable
EdgeList <- subset(Distances,Distances$Distance<0.015)[,-3]

#create VertexNames = nodes
#remove Sequences variable
#vector of unique ID values from IDs, character strings, not factored
VertexNames <- as.character(subset(IDs,IDs$ID %in% unique(IDs$ID))[,-2])


#INITIALIZE NETWORK

#intitialize network object with number of vertices, from VertexNames
#n=the number of vertices to initialize
#directed=F, should edges be interpreted as directed (no, edge indicated "close" genetic distance)
#loops=F, should loops be allowed (no, don't want connections to self)
library(network)
GeneticNetwork <- network.initialize(n=length(VertexNames),
                                     directed=FALSE,
                                     loops=FALSE)
GeneticNetwork


#VERTEX NAMES

#set vertex names with VertexNames
network.vertex.names(GeneticNetwork) <- VertexNames
network.vertex.names(GeneticNetwork)


#ADD EDGES

#check that both IDs in edge list pairs are in the ID list (vertex names)
setdiff(EdgeList[,1],VertexNames)
setdiff(EdgeList[,2],VertexNames)

#remove pair that contains lab strain (HXB2_prrt), not in vertex names
EdgeList <- EdgeList[-1022,]

#recheck IDs in edgelist
#check that both IDs in edge list pairs are in the ID list (vertex names)
setdiff(EdgeList[,1],VertexNames)
setdiff(EdgeList[,2],VertexNames)

#add edges to network
GeneticNetwork[as.matrix(EdgeList)] <- 1


#PLOT
plot(GeneticNetwork)
