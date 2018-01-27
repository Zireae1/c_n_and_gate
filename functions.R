### fuction for plotting single network
plot_network<-function(filename){ # paste("data/",flist[1], sep="")
  links<-fread(filename, header = F, stringsAsFactors = F)
  nodes<-unique(as.vector(as.matrix(links)))
  
  ### create net
  net <- graph_from_data_frame(d=links, vertices=nodes, directed=T) 
  net <- simplify(net, remove.multiple = F, remove.loops = F) 
  
  ### set different colors for Carbon and Nitrogen
  nodes<-cbind(nodes,lapply(nodes, function(x) strsplit(x, split = "")[[1]][1]))
  V(net)$type <- nodes[,2] %in% "N" 
  col <- c("red", "royal blue")
  
  ### plot network
  plot_net<-plot(net, edge.arrow.size=1, vertex.label.cex=1, edge.curved=curve_multiple(net), 
                 vertex.label.color="black", layout=layout.bipartite, 
                 vertex.size=30, edge.color="black", edge.width=2, edge.arrow.size=1,
                 vertex.color = col[as.numeric(V(net)$type) + 1])
  return(plot_net)
}

### fixed autocurve function for curving multiple edges on graph
autocurve.edges2 <-function (graph, start = 0.5)
{
  cm <- count.multiple(graph)
  mut <-is.mutual(graph)  #are connections mutual?
  el <- apply(get.edgelist(graph, names = FALSE), 1, paste,
              collapse = ":")
  ord <- order(el)
  res <- numeric(length(ord))
  p <- 1
  while (p <= length(res)) {
    m <- cm[ord[p]]
    mut.obs <-mut[ord[p]] #are the connections mutual for this point?
    idx <- p:(p + m - 1)
    if (m == 1 & mut.obs==FALSE) { #no mutual conn = no curve
      r <- 0
    }
    else {
      r <- seq(-start, start, length = m)
    }
    res[ord[idx]] <- r
    p <- p + m
  }
  res
}
