rm(list=ls(all=TRUE))
getwd()
#update.packages(ask = FALSE, dependencies = c('Suggests')) #update all installed R packages

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
library(scales)         # for points transparency on plots
#install.packages('igraph')
library('igraph')       # for bipartite plots and graph analysis
library(psych)          # for caclulation of trace for a matrix

setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
source("functions.R")    # some custom functions

### list files with network in data/x_x folder
dir<-c("3_3")
flist<-list.files(paste("data/", dir, "/", sep=""),pattern="Network_")

###################################
###   Single network analysis   ###
###################################
### construct set of nodes (number of C and N) from name of the input directory
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
  ### try to use abundance and concentration data
  # name<-gsub("Network_", "", file)
  # resdata<-paste("data/", dir,"/ResourceData_", name, sep="")
  # if(file.exists(resdata)){
  #   resources<-fread(resdata, header = F, stringsAsFactors = F)
  #   nodes<-cbind(nodes, resources$V2)
  # }
  # bacdata<-paste("data/", dir,"/SpeciesData_", name, sep="")
  # if(file.exists(bacdata)){
  #   bacteria<-fread(bacdata, header = F, stringsAsFactors = F)
  #   bacteria<-bacteria[which(bacteria[,3]>0),]
  #   links<-cbind(links, round(log10(bacteria[,3])))
  # }
  #nodes<-unique(as.vector(as.matrix(links)))
  cairo_pdf(paste("graphs/", dir, "/", gsub(".txt", "", file), "_plot.pdf", sep=""), 
            width = 8, height = 8)                                 # save plot
  plot_network(links, nodes, layout=layout_nicely, vertex.size=20) # default layout is bipartite
               
  dev.off()
}

### look at cycles statistics:
counts<-data.frame(matrix(NA, nrow = length(flist), ncol = 2+length(V(net))/2+length(nodes)*2))
colnames(counts)<-c("netname", as.character(seq(2,length(V(net)), by=2)), 
                    sub("", "d_", nodes), sub("", "ind_", nodes), "mean_path")
i<-1
for (file in flist) {
  print(file)
  links<-fread(paste("data/", dir, "/", file, sep=""), header = F, stringsAsFactors = F)
  #nodes<-unique(as.vector(as.matrix(links)))
  ### create net
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  net <- simplify(net, remove.multiple = F, remove.loops = F) 
  ncyc<-count_cycles(net)
  k<-1
  counts[i, k]<-file
  k<-k+length(file)
  counts[i, k:(k+length(ncyc[,2])-1)]<-ncyc[,2]
  k<-k+length(ncyc[,2])
  counts[i, k:(k+length(degree(net))-1)]<-degree(net)
  k<-k+length(degree(net))
  counts[i, k:(k+length(degree(net))-1)]<-degree(net, mode="in")
  counts[i, ncol(counts)]<-average.path.length(net, directed=T, unconnected=T)
  #mean(degree(net, mode="out")) # Outdegree for net
  i<-i+1
}
str(counts)
### check how many cycles of each type we have:
print(rbind(ncyc[,1],colSums(counts[,2:4])))
### save barplot for cycle counts
cairo_pdf(paste("graphs/", dir, "/","Cycle_counts_plot.pdf", sep = ""), 
          width = 8, height = 8) # save plot
barplot_cyc <- barplot(as.matrix(counts[, 3:(3+length(V(net))/2-2)]), ylim = c(0,50), 
        xlab = "Cycle length", ylab="Count", 
        main = paste("Number of cycles of length N (", nrow(counts), " states)", sep=""), 
        col="royal blue")
text(x = barplot_cyc, y = colSums(counts[, 3:(3+length(V(net))/2-2)]), 
     label = paste(round(colSums(counts[, 3:(3+length(V(net))/2-2)])/nrow(counts), 
                         digits = 4)*100, "%", sep=""), 
     pos = 3, cex = 1, col = "black")
#barplot_cyc
dev.off()

### create plot for degree of nodes, path, max indegree (hubs), etc 
### TO DO: coord in row auto
cairo_pdf(paste("graphs/", dir, "/","General_stats_plot.pdf", sep = ""), 
          width = 8, height = 8) # save plot
par(mfrow = c(2,2)) # Set up a 2x2 display
hist(rowMeans(counts[7:16]), xlab = "Mean node degree", 
     main = "Mean Degree", #prob = TRUE, 
     col = "royal blue", breaks = 10)
# hist(rowMeans(counts[17:26]), xlab = "Mean indegree", 
#      main = "Mean Indegree Distribution", #prob = TRUE,  
#      col =" royal blue", breaks = 10)
hist(apply(counts[17:26], 1, max) , xlab = "Max node indegree", 
     main = "Max Indegree", #prob = TRUE, 
     col = "royal blue", breaks = 10)
hist(counts$mean_path, xlab = "Mean path length", 
     main = "Mean Path Length", #prob = TRUE, 
     col = "royal blue", breaks = 40)
hist(counts$mean_path/rowSums(counts[17:26])*10, 
     xlab = "Mean path length / Number of species", 
     main = "Normalized Mean Path Length", #prob = TRUE, 
     col = "royal blue", breaks = 20)
dev.off()
par(mfrow = c(1,1)) # Restore display



plot(density(counts$mean_path/rowSums(counts[17:26])*10))

### look how many bacterial species survived:
counts$mean_path/rowSums(counts[17:26])


plot(c(3,4,5,7),c(13/200, 8/100, 41/400, 17/100), 
     xlab="CN number", ylab="Fraction of bistable states")


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


