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

### list files with network in data/x_x folder
dir<-c("7_7")
flist<-list.files(paste("data/", dir, "/", sep=""),pattern="Network_")

###################################
###   Single network analysis   ###
###################################

cnum<-as.numeric(strsplit(dir, "_")[[1]][1])
nnum<-as.numeric(strsplit(dir, "_")[[1]][2])
nodes<-vector()
for (i in 1:cnum){
  nodes<-rbind(nodes, paste("C",i, sep=""))
}
for (i in 1:nnum){
  nodes<-rbind(nodes, paste("N",i, sep=""))
}
### save diagrams for all networks:
for (file in flist) {
  print(file)
  links<-fread(paste("data/", dir,"/", file, sep=""), header = F, stringsAsFactors = F)
  #nodes<-unique(as.vector(as.matrix(links)))
  cairo_pdf(paste("graphs/", dir, "/", gsub(".txt", "", file), "_plot.pdf", sep=""), 
            width = 8, height = 8)                             # save plot
  plot_network(links, nodes, layout=layout_nicely, # default layout is bipartite
               idc=4, idn=5, vertex.size=20)                   # ids of C and N in filename
  dev.off()
}

### look at cycles statistics:
counts<-data.frame(matrix(NA, nrow = length(flist), ncol = 6))
i<-1
for (file in flist) {
  print(file)
  links<-fread(paste("data/", dir, "/", file, sep=""), header = F, stringsAsFactors = F)
  nodes<-unique(as.vector(as.matrix(links)))
  ### create net
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  net <- simplify(net, remove.multiple = F, remove.loops = F) 
  ncyc<-count_cycles(net, nmax=10)
  counts[i,1]<-file
  counts[i,2:6]<-ncyc[,2]
  i<-i+1
}
### check how many cycles of each type we have:
print(rbind(ncyc[,1],colSums(counts[,2:6])))
### save barplot for cycle counts
colnames(counts)<-c("netname", "2", "4", "6", "8", "10")
cairo_pdf(paste("graphs/", dir, "/","Cycle_counts_plot.pdf", sep=""), 
          width = 8, height = 8) # save plot
barplot_cyc<-barplot(as.matrix(counts[,3:6]),ylim=c(0,50), 
        xlab="Cycle length", ylab="Count", 
        main = "Number of cycles of length N (445 states)", 
        col="royal blue")
text(x=barplot_cyc,y = colSums(counts[,3:6]), label=colSums(counts[,3:6]), 
     pos = 3, cex = 1, col = "black")
barplot_cyc
dev.off()

############################################################## 
### comparison of alternative states for bistable systems  ###
##############################################################
tt<-flist[grep("Alt", flist)]
stats.flist<-list.files(dir,pattern="TotalNoOfSteadyStates_")
k<-1
for (i in seq(1,length(tt), 2)){
  st<-fread(paste("data/", dir, "/", stats.flist[k], sep=""), header = F)
  ### save diagram of joint network
  cairo_pdf(paste("graphs/", dir, "/", gsub(".txt", "", tt[i]), "_alt2_plot.pdf", sep=""), 
            width = 6, height = 6)
  plot_joint_network(paste("data/", dir, "/", tt[i], sep=""), 
                     paste("data/", dir, "/",tt[i+1], sep=""))
  legend("topleft", title = "Frequency of states",legend = st$V2, col = c("firebrick", "black"), 
         lty= 1, lwd = 2, bty="n")
  dev.off()
  k<-k+1
}

### find alternative states
alt<-gsub("(.*)_", "", flist)
alt<-as.numeric(gsub(".txt", "", alt))
alt<-data.frame(flist, alt)
str(alt)
k<-1
i<-1
while(i<=nrow(alt)){
  alt_files<-vector()
  if(alt[i+1,2]=alt[i,2]+1){
    alt_files[k]<-alt[i,1]
    k<-k+1
    i<-i+1
  }else{
    k<-1
    i<-i+1}
}

### compare two networks (alternate states) for 3-stable
links1<-fread("data/Network_16_1002_5_5_1_10_100_0.1_1.0_100_500_Alt_1.txt", header = F, stringsAsFactors = F)
links2<-fread("data/Network_16_1002_5_5_1_10_100_0.1_1.0_100_500_Alt_2.txt", header = F, stringsAsFactors = F)
links3<-fread("data/Network_16_1002_5_5_1_10_100_0.1_1.0_100_500_Alt_3.txt", header = F, stringsAsFactors = F)

merged<-rbind(cbind(links1, "1"),cbind(links2,"2"),cbind(links3,"3")) # joint network 
colnames(merged)<-c("from", "to", "state")
nodes<-unique(as.vector(as.matrix(merged[,1:2])))

str(nodes)
### create network
# net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 
# net2 <- graph_from_data_frame(d=links2, vertices=nodes, directed=T) 
net3 <- graph_from_data_frame(d=links3, vertices=nodes, directed=T) 

ncyc<-count_cycles(net3, nmax=10)
net <- graph_from_data_frame(d=merged, vertices=nodes, directed=T) 

### set different colors for Carbon and Nitrogen
nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
V(net)$type <- nodes[,2] %in% "N" 
col <- c("red", "royal blue") # gold, etc

### set defferent colors for bacteria in different states
ecol <- rep("black", ecount(net))
ecol[which(E(net)$state=="2")] <- "firebrick"
ecol[which(E(net)$state=="3")] <- "olivedrab3"

st<-c(123,118,259)
### plot network
cairo_pdf(paste("graphs/", "Network_16_1002_5_5_1_10_100_0.1_1.0_100_500_Alt1_alt2_alt3_plot.pdf", 
                sep=""), width = 6, height = 6)

plot(net, edge.arrow.size=1, vertex.label.cex=1, 
     #edge.curved=seq(-0.5, 0.5, length = ecount(net)), 
     edge.curved=autocurve.edges2(net, start=0.4),
     vertex.label.color="black", #layout=layout.bipartite, 
     vertex.size=20, edge.color=ecol, 
     edge.width=2, edge.arrow.size=1,
     asp = 0.9,
     vertex.color = col[as.numeric(V(net)$type) + 1])
legend("topleft", title = "Frequency of states",legend = st, col = c("black", "firebrick","olivedrab3"), 
       lty= 1, lwd = 2, bty="n")
dev.off()
### redo 5_5_6_1

# #### some exp with calculation of a loop length
# links1<-fread("data/Network_5_5_1_2_alt1.txt", header = F, stringsAsFactors = F)
# # add cycles
# links1<-rbind(links1, t(c("N3", "C1"))) # length 2
# 
# links1<-rbind(links1, t(c("N2", "C2"))) # length 4
# links1<-rbind(links1, t(c("N5", "C3"))) # length 2
# 
# links1<-rbind(links1, t(c("C4", "N5"))) # length 2
# 
# links1<-rbind(links1, t(c("C5", "N4"))) # length 6 (1)
# links1<-rbind(links1, t(c("C2", "N5"))) # length 6 (2)
# nodes<-unique(as.vector(as.matrix(links1)))
# nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
# 
# net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 
# 
# V(net1)$type <- nodes[,2] %in% "N" 
# 
# #cairo_pdf(paste("graphs/", "Network_5_5_34_3_alt1_plot1.pdf", sep=""), width = 5, height = 5)
# plot(net1, edge.arrow.size=1, vertex.label.cex=1, #edge.curved=curve_multiple(net), 
#      vertex.label.color="black", #layout=layout.bipartite, 
#      vertex.size=25, edge.color="black", edge.width=2, edge.arrow.size=1,
#      vertex.color = col[as.numeric(V(net1)$type) + 1])
# #dev.off()

# for (i in seq(1,length(tt), 2)){
#   links1<-fread(paste("data/", tt[i], sep=""), header = F, stringsAsFactors = F)
#   links2<-fread(paste("data/", tt[i+1], sep=""), header = F, stringsAsFactors = F)
#   merged<-rbind(links1,links2) # joint network 
#   nodes<-unique(as.vector(as.matrix(merged)))
#   
#   ### create network
#   net1 <- graph_from_data_frame(d=links1, vertices=nodes, directed=T) 
#   net2 <- graph_from_data_frame(d=links2, vertices=nodes, directed=T) 
#   net <- graph_from_data_frame(d=merged, vertices=nodes, directed=T) 
#   
#   ### set different colors for Carbon and Nitrogen
#   nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
#   V(net)$type <- nodes[,2] %in% "N" 
#   col <- c("red", "royal blue") # gold, etc
#   
#   ### set defferent colors for bacteria in different states
#   ecol <- rep("black", ecount(net))
#   ecol[which(E(net)%in%E(net1))] <- "firebrick"
#   
#   ### plot network
#   cairo_pdf(paste("graphs/", gsub(".txt", "", tt[i]), "_alt2_plot.pdf", sep=""), width = 6, height = 6)
#   
#   plot(net, edge.arrow.size=1, vertex.label.cex=1, 
#        #edge.curved=seq(-0.5, 0.5, length = ecount(net)), 
#        edge.curved=autocurve.edges2(net, start=0.2),
#        vertex.label.color="black", #layout=layout.bipartite, 
#        vertex.size=20, edge.color=ecol, 
#        edge.width=2, edge.arrow.size=0.4,
#        asp = 0.9,
#        vertex.color = col[as.numeric(V(net)$type) + 1])
#   legend("topleft", legend = c("state 1", "state 2"), col = c("firebrick", "black"), 
#          lty= 1, lwd = 2, bty="n")
#   dev.off()
# }


count_cycles(net1)

