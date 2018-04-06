rm(list=ls(all=TRUE))
getwd()
#update.packages(ask = FALSE, dependencies = c('Suggests')) #update all installed R packages

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
library(scales)         # for points transparency on plots
setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
source("functions.R")    # some custom functions

# load data (use binary format)
states<-fread("data_states/FinalStats_1_3_3_full_new_01.dat", header = F, stringsAsFactors = F)
head(states)

# drop half of the flux points
# states<-states[seq(1,nrow(states), 2),] 

# drop flux values for now
fluxes<-states[,1:6]
states<-states[,7:ncol(states)]
states<-as.data.frame(states)
str(states[1:20,]) # check that all values are of type numeric

colnames(states)<-c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09",
                    "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18",
                    "S19", "S20", "S21", "S22")

# prepare table for heatmap - calc intersections for two states 
counts<-matrix(0, ncol(states), ncol(states))
colnames(counts)<-colnames(states)
rownames(counts)<-colnames(states)
for (i in 1:ncol(states)){
  for (j in i:ncol(states)){
    counts[i,j]<-length(which(states[,i]==1 & states[,j]==1))
    counts[j,i]<-counts[i,j]
    print(i)
  }
}

# plot heatmap of pairwise intersections
library(gplots)
my_palette <- colorRampPalette(rev(brewer.pal(11,"PRGn")))(100)

cairo_pdf(paste("graphs/","heatmap_states_all.pdf", sep = ""), 
          width = 10, height = 10) # save plot
heatmap.2(counts, trace="none", col=my_palette, 
          dendrogram="column", #margins =c(4,4), 
          symm=F, symkey=F, symbreaks=F, srtCol=360,
          colsep=c(0,ncol(counts)), rowsep=c(0,nrow(counts)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") #dendrogram="column",  key.par=list(mar=c(12,12,14,14)), lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 )
dev.off()

#save counts files
write.table(counts, "output/counts_full.txt", sep='\t', row.names = T)
write.csv2(counts,"output/counts_full.csv",row.names = F)
counts_norm<-counts
for (i in 1:ncol(counts)){
  print(i)
  for (j in i:ncol(counts)){
    counts_norm[i,j]<-counts[i,j]/min(counts[i,i], counts[j,j])
    counts_norm[j,i]<-counts_norm[i,j]
  }
}
#save normalized counts files
write.table(counts_norm, "output/counts_full_norm.txt", sep='\t', row.names = T)
write.csv2(counts_norm,"output/counts_full_norm.csv",row.names = F)

# plot heatmap of pairwise intersections for lognormal values
counts<-counts+1
cairo_pdf(paste("graphs/","heatmap_states_all_log.pdf", sep = ""), 
          width = 10, height = 10) # save plot
heatmap.2(log10(counts), trace="none", col=my_palette, 
          dendrogram="column", #margins =c(4,4), 
          symm=F, symkey=F, symbreaks=F, srtCol=360,
          colsep=c(0,ncol(counts)), rowsep=c(0,nrow(counts)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") #dendrogram="column",  key.par=list(mar=c(12,12,14,14)), lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 )
dev.off()

counts_norm<-counts
for (i in 1:ncol(counts)){
  print(i)
  for (j in i:ncol(counts)){
    counts_norm[i,j]<-counts[i,j]/min(counts[i,i], counts[j,j])
    counts_norm[j,i]<-counts_norm[i,j]
  }
}

cairo_pdf(paste("graphs/","heatmap_states_all_normalized.pdf", sep = ""), 
          width = 10, height = 10) # save plot
heatmap.2(log10(counts_norm), trace="none", col=my_palette, 
          dendrogram="column", #margins =c(4,4), 
          symm=F, symkey=F, symbreaks=F, srtCol=360,
          colsep=c(0,ncol(counts)), rowsep=c(0,nrow(counts)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") #dendrogram="column",  key.par=list(mar=c(12,12,14,14)), lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 )
dev.off()

# look at state abundance distribution
vv<-vector()
for(i in 1:ncol(counts)){
  vv[i]<-counts[i,i]
}
vv<-t(as.data.frame(vv))
colnames(vv)<-colnames(counts)
vv<-vv[,order(vv[1,], decreasing = T)]
labs<-gsub("S", "", names(vv))
cairo_pdf(paste("graphs/","State_abundance_distribultion.pdf", sep = ""), 
          width = 9, height = 5) # save plot
plot(log10(vv), xaxt = 'n', ylab="log10(state frequency)", 
     xlab = "States", col="royal blue", pch=16, 
     main="Log distribution of state abundances")
axis(1, at=1:22, labels=labs)
dev.off()

### check how many times we have 1 or 2 or 3 etc states simultaneously
sums<-rowSums(states)
hist(sums)
length(which(sums==1))
length(which(sums==2)) # why this one is 0????
length(which(sums==3))
length(which(sums%in%c(3:22)))
states[1:100,]

rowSums(t(states))

# scatter plot for corr between number of states and the phase volume
species<-fread("data_states/FinalStats_1_3_3_100000.txt", header = F, stringsAsFactors = F)
cairo_pdf(paste("graphs/","Cor_number_of_species_vs_state_abund.pdf", sep = ""), 
          width = 6, height = 6) # save plot
plot(species$V10, log10(vv[order(names(vv))]), xlab = "number of species", xlim=c(2,7), 
     ylab="Log10 state abundance", main = "State phase volume vs number of species")
legend("topright", pch=4, legend = 
         paste("cor = ",round(cor(species$V10, log10(vv[order(names(vv))])), digits = 3), sep=""),bty="n")
dev.off()

# construct the network
nodes<-colnames(counts)
links<-data.frame()
library(reshape2)

counts_x<-counts-1
for(i in 1:ncol(counts)){
  counts_x[i,i]<-0
}
xxx<-melt(counts_x)
#xxx$value<-log10(xxx$value)
links<-xxx[which(xxx$value>1),]
links$value<-log10(links$value)

# drop single nodes
nodes<-sort(intersect(unique(links$Var1)[1:11], colnames(counts)))

library(igraph)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T) 
E(net)$width<-(round(links$value)-3)
E(net)$color<-"black"
vv_1<-vv[which(names(vv)%in%nodes)]
V(net)$size<-(round(log10(vv_1[order(names(vv_1))])*10/3, digits = 0)-13)*3

# plot in a nice way without singletons
cairo_pdf(paste("graphs/","Net_of_states_nice.pdf", sep = ""),
          width = 7, height = 7) # save plot

plot(net, edge.arrow.size=.0, vertex.label.cex=0.8, veredge.curved=.8, edge.width=E(net)$width,
     vertex.label.color="black", vertex.label.cex=0.5, layout=layout_nicely, vertex.size=V(net)$size)

dev.off()

# plot circular diagram
cairo_pdf(paste("graphs/","Net_of_states_circle_v2.pdf", sep = ""), 
          width = 8, height = 8) # save plot
plot(net, edge.arrow.size=.0, vertex.label.cex=0.8, veredge.curved=.8, edge.width=E(net)$width,
     edge.color=E(net)$color,
       vertex.label.color="black", vertex.label.cex=0.5, layout=layout_in_circle, vertex.size=V(net)$size)

dev.off()

### save network table for Sergey:
vv_1<-vv[order(names(vv))]
colnames(links)<-c("state1", "state2", "V_intersection")
nodes<-cbind(nodes, vv_1)
write.table(nodes, "output/network_nodes.txt", sep='\t', row.names = F, quote = F)

### remove duplicates from links network
links_x<-links
links_x$state1<-as.numeric(gsub("S", "", links_x$state1))
links_x$state2<-as.numeric(gsub("S", "", links_x$state2))
sub.links<-data.frame()
for (i in 1:nrow(links_x)){
  #i=1
  if(links_x$state1[i]<links_x$state2[i]){
    print(i)
    sub.links<-rbind(sub.links, links[i,])
  }
}
write.table(sub.links, "output/network_links.txt", sep='\t', row.names = F, quote = F)

### prepare networks for individual states
nodes<-c("C1", "C2", "C3", "N1", "N2", "N3")
list.states<-list()
for (i in 1:nrow(species)){
  print(i)
  #i=1
  sp.mat<-t(matrix(unlist(species[i, 1:9]), nrow=3, ncol = 3))
  colnames(sp.mat)<-c("C1", "C2", "C3")
  rownames(sp.mat)<-c("N1", "N2", "N3")
  sp.mat<-melt(sp.mat)
  sp.mat$Var1<-as.character(sp.mat$Var1)
  sp.mat$Var2<-as.character(sp.mat$Var2)
  str(sp.mat)
  sp.mat<-sp.mat[which(sp.mat$value!=0),]
  for (j in which(sp.mat$value==1)){
    #j=1
    tmp<-sp.mat$Var1[j]
    sp.mat$Var1[j]<-sp.mat$Var2[j]
    sp.mat$Var2[j]<-tmp
  }
  
  list.states[[i]]<-sp.mat
}
#list.states[[22]]

for(i in 1:nrow(species)){
  #i=1
  cairo_pdf(paste("graphs/Netw_",i, "_plot.pdf", sep=""), 
            width = 8, height = 8)                                 # save plot
  plot_network(list.states[[i]], nodes, layout=layout_nicely, vertex.size=20) # default layout is bipartite
  
  dev.off()
}

### joined networks for multistable cases
multi.states<-states[which(sums>1),]
comb.states<-unique(multi.states)

for (i in 1:nrow(comb.states)){
  #i=1
  ids<-which(comb.states[i,]>0)
  merged<-data.frame()
  k<-1
  for (j in ids){
    merged<-rbind(merged, cbind(list.states[[j]][,1:2],k))
    k<-k+1
    print(j)
  }
  rowSums(comb.states)
  colnames(merged)<-c("from", "to", "state")
  vv<-vv[order(names(vv))]
  st<-as.data.frame(round(vv[ids]/85766121, digits = 2))
  cairo_pdf(paste("graphs/Netw_joint_", ids[1],"_", ids[2],"_", ids[3],"_", "joint_plot.pdf", 
                  sep=""), width = 9, height = 9)
  plot_joint_network(merged, nodes, st, layout=layout_nicely, vertex.size=15, start=0.4)
  dev.off()
}


#### TO DO: network of bac intersections between states
merged
d<-as.matrix(dist(species[,1:9], method = "manhattan" ))
colnames(d)<-names(vv)
rownames(d)<-names(vv)

# check whether two states have an intersection or not
intersect<-matrix(0, ncol(d), ncol(d))
for (i in 1:ncol(d)){
  for (j in j:ncol(d)){
    i=1
    j=1
    intersect[i,j]<- comb.states[,i] comb.states[,j]
  }
  
}
comb.states
group_col<-character(nrow(d))

########## test piece
which(comb.states)
group_col[which(rownames(ProjGenTop)%in%alco)]<-"orangered"
group_col[which(rownames(ProjGenTop)%in%cirr)]<-"olivedrab3"


cairo_pdf(paste("graphs/Heatmap_states_structure_manhattan.pdf", 
                sep=""), width = 9, height = 9)
heatmap.2(d, trace="none", col=my_palette, 
          dendrogram="column", #margins =c(4,4), 
          symm=F, symkey=F, symbreaks=F, srtCol=360, ColSideColors=group_col,
          colsep=c(0,ncol(d)), rowsep=c(0,nrow(d)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") #dendrogram="column",  key.par=list(mar=c(12,12,14,14)), lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 )
legend("bottomleft", legend = c("ADS", "ALC"), col = c("orangered", "olivedrab3"), 
       lty= 1, lwd = 10, bty="n")
dev.off()

############################### Leftovers
#states[1:20,]
#states[1,1]
#states[which(states==states[1,1])]<-NA
#unique(as.vector(as.matrix(states)))

# convert to 0/1
#states[!is.na(states)]<-1
#states[is.na(states)]<-0

# for (i in 1:ncol(states)){ ## make all columns of numeric type
#   #i=5
#   tt<-states[,i]
#   tt<-as.numeric(unlist(tt))
#   states[,i]<-tt
# }
# head(states)

#### interactive vizNetwork package (not good)
# colnames(nodes)<-c("id")
# 
# nodes<-as.data.frame(nodes)
# nodes$shape <- "dot"  
# nodes$shadow <- TRUE # Nodes will drop shadow
# nodes$title <- nodes$id # Text on click
# nodes$size <- 15
# nodes$font.size <- 25
# nodes$borderWidth <- 2 # Node border width
# 
# links<-as.data.frame(links)
# colnames(links)[1:2]<-c("from", "to")
# links$width <- round(links$value)-3
# links$color <- "gray"    # line color  
# #links$arrows <- "from" # arrows: 'from', 'to', or 'middle'
# links$smooth <- FALSE    # should the edges be curved?
# links$shadow <- FALSE    # edge shadow
# 
# library(visNetwork)
# network <- visNetwork(nodes, links, height = "800px", width = "100%") %>% 
#   #visPhysics(stabilization = list(iterations = 10)) %>%
#   visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE,
#              manipulation = TRUE) %>% visLegend()
# network
# 
# visSave(network, file = "interactive_network_of_states.html")


#sum(vv)

# mydata_hist <- hist(sums, breaks=10, plot=FALSE)
# tmp<-cbind(c(1:6), mydata_hist$count+1)
# bp <- barplot(mydata_hist$count+1, log="y", col="white")
# text(bp, mydata_hist$counts, labels=c(1:6), pos=1)


# frequency(sums)
# head(states[which(sums==6),])
# hist(log10(vv), breaks=20)

# plot cluster dendrogram which shows how similar are the states
# d <- dist(t(states), method = "euclidean")
# hc <- hclust(d, method = "ward.D2")
# cairo_pdf(paste("graphs/","clusters.pdf", sep = ""), 
#           width = 6, height = 6) # save plot
# plot(hc, cex = 0.6, hang = -1, 
#      xlab="States")
# dev.off()







#source("https://bioconductor.org/biocLite.R")
#biocLite("OmicCircos")

#library(OmicCircos)
# my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(200)
# library(gplots)
# heatmap.2(t(as.matrix(t(states))), col=my_palette, density.info="none", trace="none", 
#           Rowv=NA, dendrogram="row", margins=c(5,9),ColSideColors=group_col, #srtCol=60, 
#           symm=F,symkey=F,symbreaks=F, 
#           colsep=c(0,nrow(ProjGenTop)),rowsep=c(0,ncol(ProjGenTop)), 
#           sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0",
#           key.xlab = "", key.title = "Relative abundance, %")


