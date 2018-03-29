rm(list=ls(all=TRUE))
getwd()
#update.packages(ask = FALSE, dependencies = c('Suggests')) #update all installed R packages

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
library(scales)         # for points transparency on plots
setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
#source("functions.R")    # some custom functions

# load data
states<-fread("data_states/FinalStats_1_3_3_full_new_01.dat", header = F, stringsAsFactors = F)
head(states)

# drop half of the flux points
states<-states[seq(1,nrow(states), 2),] 

# drop flux values for now
states<-states[,7:ncol(states)]
states<-as.data.frame(states)
#states[1:20,]
#states[1,1]
#states[which(states==states[1,1])]<-NA
#unique(as.vector(as.matrix(states)))

# convert to 0/1
#states[!is.na(states)]<-1
#states[is.na(states)]<-0
str(states[1:20,])
head(states)

colnames(states)<-c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09",
                    "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18",
                    "S19", "S20", "S21", "S22")

# for (i in 1:ncol(states)){ ## make all columns of numeric type
#   #i=5
#   tt<-states[,i]
#   tt<-as.numeric(unlist(tt))
#   states[,i]<-tt
# }
# head(states)

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
plot(species$V3, log10(vv[order(names(vv))]), xlab = "number of species", xlim=c(2,7), 
     ylab="Log10 state abundance", main = "State phase volume vs number of species")
legend("topright", pch=4, legend = 
         paste("cor = ",round(cor(species$V3, log10(vv[order(names(vv))])), digits = 3), sep=""),bty="n")
dev.off()

# construct the network
nodes<-colnames(counts)
links<-data.frame()
library(reshape2)

counts_x<-counts
for(i in 1:ncol(counts)){
  counts_x[i,i]<-0
}
xxx<-melt(counts_x)
#xxx$value<-log10(xxx$value)
links<-xxx[which(xxx$value>1),]
links$value<-log10(links$value)
library(igraph)
net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net <- simplify(net, remove.multiple = F, remove.loops = T) 
E(net)$width<-round(links$value)-3
V(net)$size<-round(log10(vv[order(names(vv))]))*2.5
# cairo_pdf(paste("graphs/","Net_of_states_nice.pdf", sep = ""), 
#           width = 10, height = 10) # save plot
# 
# plot(net, edge.arrow.size=.0, vertex.label.cex=0.8, veredge.curved=.8, edge.width=E(net)$width,
#      vertex.label.color="black", vertex.label.cex=0.5, layout=layout_nicely, vertex.size=0.1)
# 
# dev.off()
cairo_pdf(paste("graphs/","Net_of_states_circle.pdf", sep = ""), 
          width = 7, height = 7) # save plot

plot(net, edge.arrow.size=.0, vertex.label.cex=0.8, veredge.curved=.8, edge.width=E(net)$width,
       vertex.label.color="black", vertex.label.cex=0.5, layout=layout_in_circle, vertex.size=V(net)$size)

dev.off()

#### TO DO: network of bac intersections between states





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


