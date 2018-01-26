rm(list=ls(all=TRUE))
getwd()

library(data.table)     # for fread function
library(RColorBrewer)   # for colour palette
library(scales)         # for points transparency on plots
#install.packages('igraph')
library('igraph')       # for bipartite plots and graph analysis
library(psych)          # for caclulation of trace for a matrix

setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
source("functions.R")    # some useful functions

### read data
links<-fread("data/Network_5_5_20_7.txt", header = F, stringsAsFactors = F)
nodes<-unique(as.vector(as.matrix(links)))

### create network
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = F) 

### set different colors for Carbon and Nitrogen
nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
V(net)$type <- nodes[,2] %in% "N" 
col <- c("red", "royal blue")

### plot network
plot(net, edge.arrow.size=1, vertex.label.cex=1, edge.curved=curve_multiple(net), 
     vertex.label.color="black", layout=layout.bipartite, 
     vertex.size=30, edge.color="black", edge.width=2, edge.arrow.size=1,
     vertex.color = col[as.numeric(V(net)$type) + 1])

#print(net)

### compare two alternate states:
links1<-fread("data/Network_5_5_3_2_alt1.txt", header = F, stringsAsFactors = F)
links2<-fread("data/Network_5_5_3_2_alt2.txt", header = F, stringsAsFactors = F)
merged<-rbind(links1,links2) # joint network 
nodes<-unique(as.vector(as.matrix(merged)))

### create network
net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 
net2 <- graph_from_data_frame(d=links2, vertices=nodes, directed=T) 
net <- graph_from_data_frame(d=merged, vertices=nodes, directed=T) 

### set different colors for Carbon and Nitrogen
nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
V(net)$type <- nodes[,2] %in% "N" 
col <- c("red", "royal blue")

### set defferent colors for bacteria in different states
ecol <- rep("black", ecount(net))
ecol[which(E(net)%in%E(net1))] <- "orange"

### plot network
plot(net, edge.arrow.size=1, vertex.label.cex=1, 
     #edge.curved=seq(-0.5, 0.5, length = ecount(net)), 
     edge.curved=autocurve.edges2(net, start=0.2),
     vertex.label.color="black", layout=layout.bipartite, 
     vertex.size=30, edge.color=ecol, 
     edge.width=2, edge.arrow.size=1,
     asp = 0.9,
     vertex.color = col[as.numeric(V(net)$type) + 1])

#### some exp with calculation of a loop length
links1<-fread("data/Network_5_5_3_2_alt1.txt", header = F, stringsAsFactors = F)
# add cycles
links1<-rbind(links1, t(c("C1", "N3"))) # length 2
#links1<-rbind(links1, t(c("C5", "N5"))) # length 4
links1<-rbind(links1, t(c("C2", "N4"))) # length 2
links1<-rbind(links1, t(c("C4", "N5"))) # length 2

links1<-rbind(links1, t(c("C5", "N4"))) # length 6 (1)
links1<-rbind(links1, t(c("C2", "N5"))) # length 6 (2)

net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 

V(net1)$type <- nodes[,2] %in% "N" 
plot(net1, edge.arrow.size=1, vertex.label.cex=1, #edge.curved=curve_multiple(net), 
     vertex.label.color="black", layout=layout.bipartite, 
     vertex.size=30, edge.color="black", edge.width=2, edge.arrow.size=1,
     vertex.color = col[as.numeric(V(net1)$type) + 1])


### make adjacency matrix
lm<-laplacian_matrix(net1) # matrix representation
nm<-as.matrix(lm)
nm[which(nm>=0)]<-0
nm[which(nm<0)]<-1

tr(nm%*%nm) #find trace of a matrix a^2
tr(nm%*%nm%*%nm%*%nm) #find trace of a matrix a^4
tr(nm%*%nm%*%nm%*%nm%*%nm%*%nm) #find trace of a matrix a^6
tr(nm%*%nm%*%nm%*%nm%*%nm%*%nm%*%nm%*%nm%) #find trace of a matrix a^8

