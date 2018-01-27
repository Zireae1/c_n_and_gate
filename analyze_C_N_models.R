rm(list=ls(all=TRUE))
getwd()

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
library(scales)         # for points transparency on plots
#install.packages('igraph')
library('igraph')       # for bipartite plots and graph analysis
library(psych)          # for caclulation of trace for a matrix

setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
source("functions.R")    # some custom functions

### read data
links<-fread("data/Network_5_5_34_3_alt2.txt", header = F, stringsAsFactors = F) #how to read just one
flist<-list.files("data",pattern="Network_")

### plot all single networks:
for (file in flist) {
  print(file)
  cairo_pdf(paste("graphs/", gsub(".txt", "", file), "_plot.pdf", sep=""), width = 5, height = 5) # save plot
  plot_network(paste("data/", file, sep=""))
  dev.off()
}

### compare two networks (alternate states)
links1<-fread("data/Network_5_5_34_3_alt1.txt", header = F, stringsAsFactors = F)
links2<-fread("data/Network_5_5_34_3_alt2.txt", header = F, stringsAsFactors = F)
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
ecol[which(E(net)%in%E(net1))] <- "firebrick"

### plot network
cairo_pdf(paste("graphs/", "Network_5_5_34_3_alt1_alt2.pdf", sep=""), width = 5, height = 5)

plot(net, edge.arrow.size=1, vertex.label.cex=1, 
     #edge.curved=seq(-0.5, 0.5, length = ecount(net)), 
     edge.curved=autocurve.edges2(net, start=0.2),
     vertex.label.color="black", #layout=layout.bipartite, 
     vertex.size=20, edge.color=ecol, 
     edge.width=2, edge.arrow.size=1,
     asp = 0.9,
     vertex.color = col[as.numeric(V(net)$type) + 1])
dev.off()

#### some exp with calculation of a loop length
links1<-fread("data/Network_5_5_34_3_alt1.txt", header = F, stringsAsFactors = F)
# add cycles
links1<-rbind(links1, t(c("C1", "N3"))) # length 2
links1<-rbind(links1, t(c("C5", "N5"))) # length 4
links1<-rbind(links1, t(c("C2", "N4"))) # length 2
links1<-rbind(links1, t(c("C4", "N5"))) # length 2

links1<-rbind(links1, t(c("C5", "N4"))) # length 6 (1)
links1<-rbind(links1, t(c("C2", "N5"))) # length 6 (2)
nodes<-unique(as.vector(as.matrix(links1)))
nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))

net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 

V(net1)$type <- nodes[,2] %in% "N" 

cairo_pdf(paste("graphs/", "Network_5_5_34_3_alt1_plot1.pdf", sep=""), width = 5, height = 5)
plot(net1, edge.arrow.size=1, vertex.label.cex=1, #edge.curved=curve_multiple(net), 
     vertex.label.color="black", #layout=layout.bipartite, 
     vertex.size=25, edge.color="black", edge.width=2, edge.arrow.size=1,
     vertex.color = col[as.numeric(V(net1)$type) + 1])
dev.off()

### make adjacency matrix
lm<-laplacian_matrix(net1) # matrix representation
A<-as.matrix(lm)
A[which(A>=0)]<-0
A[which(A<0)]<-1

tr(A%*%A) #find trace of a matrix a^2
tr(A%*%A%*%A%*%A) #find trace of a matrix a^4
tr(A%*%A%*%A%*%A%*%A%*%A) #find trace of a matrix a^6
tr(A%*%A%*%A%*%A%*%A%*%A%*%A%*%A) #find trace of a matrix a^8

tr(A%*%A%*%A%*%A-A%*%A)/4

count_cycles<-function(net){
  ### create adjacency matrix
  lm<-laplacian_matrix(net1) # matrix representation
  A<-as.matrix(lm)
  A[which(A>=0)]<-0
  A[which(A<0)]<-1
  
  # find max possible length of a cycle
  nmax<-length(V(net1))
  
  # create table for cycle counts 
  cycles<-matrix(0,nmax/2,1)
  
  while (i <= nmax){ # while we haven't exceeded number of nodes
    if (i%%2==0){ # for all even numbers
      A
    }
  }
  A2<-A%*%A
  A4<-A2%*%A2
  A6<-A4%*%A2
  n2<-A2/2
  n4<-(A4-A2)/4
  n6<-(A6-A2)/6
}

