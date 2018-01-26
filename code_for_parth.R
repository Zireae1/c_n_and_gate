rm(list=ls(all=TRUE))
getwd()

library(data.table) # for fread function
library(RColorBrewer) # for colour palette
library(scales) # for points transparency on plots
#library(plyr) 
#install.packages('igraph')
library('igraph')
#install.packages('visNetwork')
library('visNetwork') 

setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory

# links<-fread("ResContrToSpeciesBiomass_1_0_1_0.0_20_10000_10000.txt", header = F, stringsAsFactors = F)
# links<-fread("ResContrToSpeciesBiomass_1_0_1_0.001_20_10000_10000.txt", header = F, stringsAsFactors = F)
# #colnames(links)<-c("microbe", "metab", "contribution")
# links$microbe<-sub("", "b", links$microbe)
# links$metab<-sub("", "c", links$metab)


### new tables
links<-fread("ResContrToSpeciesBiomass_1_0.1_1_0.001_0_20000_20000.txt", header = T, stringsAsFactors = F)
links<-fread("ResContrToSpeciesBiomass_1_0.1_2_0.0_0_20000_20000.txt", header = T, stringsAsFactors = F)

links<-fread("ResContrToSpeciesBiomass_4_0.1_4_0.0_0_20000_20000.txt", header = T, stringsAsFactors = F)


### load chemistry structure
chem<-fread("chemistry_tree.txt", header=T, stringsAsFactors = F)
chem$level<-6-chem$level
chem$resource<-sub("", "r", chem$resource)
chem$`by-products`<-sub("", "r", chem$`by-products`)

#links<-fread("AbundAndGrowthRates_2_0_1_0.0_20_10000_10000.txt", header = F, stringsAsFactors = F)
links$SpeciesID<-sub("", "s", links$SpeciesID)
links$ResourceID<-sub("", "r", links$ResourceID)

bac<-sort(unique(links$SpeciesID))
comp<-sort(unique(links$ResourceID))
nodes<-rbind(cbind(bac, "bac"), cbind(comp, "comp"))

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
net

#plot(net)
net <- simplify(net, remove.multiple = F, remove.loops = T) 

colrs <- c("tomato", "gold")
V(net)$type<-as.numeric(as.factor(nodes[,2]))
V(net)$type<-nodes[,2]

V(net)$color <- colrs[as.numeric(V(net)$type)]
deg <- degree(net, mode="all")
V(net)$size <- deg*0.5

V(net)

#cairo_pdf('graphs/test_firmicutes.pdf', width = 11, height = 10) # save plot
#cairo_pdf('test_set.pdf', width = 13, height = 13) # save plot

plot(net, edge.arrow.size=.4, vertex.label.cex=0.8, veredge.curved=.5, 
     vertex.label.color="black", vertex.label.cex=0.5, layout=layout_nicely)

#dev.off()

### for bipartite visualization:
V(net)$type <- FALSE
V(net)$type[V(net)$name%in%nodes[which(nodes[,2]=="comp"),1]] <- TRUE

l<-layout_as_bipartite(net)
#l <- cbind(1:vcount(net), c(1, nodes[which(nodes[,2]=="bac"),1], nodes[which(nodes[,2]=="comp"),1]))
#plot(net, edge.arrow.size=.4, vertex.label.cex=0.8, veredge.curved=.5, 
#     vertex.label.color="black", vertex.label.cex=0.5, layout=l)
#is.bipartite(net)
#length(which(nodes[,2]=="comp"))
# V(net)
# plot(g,
#      layout = l,
#      vertex.color = col[V(g)$type],
#      vertex.shape = shape[V(g)$type],
#      edge.width = E(g)$weights * 5 # optional, plot edges width proportional to weights.
# )


#### vizNetwork package
colnames(nodes)<-c("id", "type")

nodes<-as.data.frame(nodes)
nodes$shape <- "dot"  
nodes$shadow <- TRUE # Nodes will drop shadow
nodes$title <- nodes$bac # Text on click
#nodes$label <- nodes$type.label # Node label
#nodes$size <- degree(net, mode="all")*5
nodes$size <- 1
for (i in nodes$id[nodes$type=="bac"]){
  nodes$size[nodes$id==i]<-as.numeric(unique(links$SpeciesBiomass[links$from==i]))
}
nodes$size<-((3+log10(nodes$size))*3)

for (i in chem$resource){
  nodes$size[nodes$id==i]<-as.numeric(unique(chem$level[chem$resource==i]))*4
}
#nodes$size<-((3+log10(nodes$size))*3)

nodes$font.size <- 25
nodes$borderWidth <- 2 # Node border width
nodes$color.background <- c("tomato", "gold")[as.numeric(as.factor(nodes[,2]))]
nodes$color.border <- "black"
nodes$color.highlight.background <- "orange"
nodes$color.highlight.border <- "darkred"

links<-as.data.frame(links)
colnames(links)[1:2]<-c("from", "to")

#links$width <- as.vector(round(bps_perc, digits = 1)*10) # line width
links$width <- links$Affinity/10
#links$width <-(8+log10(links$type))*0.5
links$color <- "gray"    # line color  
links$arrows <- "from" # arrows: 'from', 'to', or 'middle'
links$smooth <- FALSE    # should the edges be curved?
links$shadow <- FALSE    # edge shadow


#network<-visNetwork(nodes, links)
network <- visNetwork(nodes, links, height = "800px", width = "100%") %>% 
  #visPhysics(stabilization = list(iterations = 10)) %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE,
             manipulation = TRUE) %>% visLegend()
network

visSave(network, file = "network_without_cost.html")


### check how rank of resource is correlated with abundance
subs<-links[which(links$to!="r1"),]

subs<-links


chem$level<-6-chem$level
subs$torank<-5
for(i in intersect(subs$to,chem$resource)){
  #i="r10"
  subs$torank[subs$to==i]<-as.numeric(unique(chem$level[chem$resource==i]))
}

count<-vector()
for(i in 1:5){
  count[i]<-length(which(subs$torank==i))/(2^(i-1))
}
count

plot(log10(c(1:5)), log10(count), pch=16)
plot(log10(count), pch=16)
model <- lm(log10(count) ~ log10(c(1:5)))
coef(model)
paste('y =', coef(model)[[2]], '* x', '+', coef(model)[[1]])

plot(subs$torank, log10(subs$SpeciesBiomass), pch=16, xlim = c(1,5))

### contribution to biomass uncorrelated with rank of the resource
subs[which(subs$torank==5),]
plot(subs$torank, subs$ContributionToBiomass, pch=16, xlim = c(1,5))


#hist(8+log10(links$width))
#library(rTRM) # somehow not working =(
#source("https://bioconductor.org/biocLite.R")
#biocLite("rTRM")
#l <- layout.concentric(net, concentric = list(bac, comp))
