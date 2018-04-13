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
states[1:20,]

# drop part of the flux points
# states<-states[seq(1,nrow(states), 2),] 
ids<-sample(85766121, size=800000, replace = F)
ids[1:10]
fluxes<-states[ids,1:6]
sub.states<-states[ids, 7:ncol(states)]
fluxes<-states[,1:6]
sub.states<-states[, 7:ncol(states)]


colnames(sub.states)<-c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09",
                    "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18",
                    "S19", "S20", "S21", "S22")
sub.states<-as.data.frame(sub.states)

fc.means<-vector()
fn.means<-vector()
for (i in 1:ncol(sub.states)){
  print(i)
  st.ids<-which(sub.states[,i]>0)
  sub.flux<-fluxes[st.ids,]
  fc<-rowMeans(sub.flux[,1:3])
  fn<-rowMeans(sub.flux[,4:6])
  fc.means[i]<-mean(fc)
  fn.means[i]<-mean(fn)
}

### save plot for Centroids of all states in <Fc> <Fn>
cairo_pdf('graphs/state_centroids_3x3.pdf', width = 7, height = 6) # save plot
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # change graph params for legend 
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(fc.means))
plot(fc.means, fn.means, col = my_palette, pch=16, xlim=c(0,491), ylim=c(0,491), 
     cex=1.3, xlab = "Mean(F_c)", ylab="Mean(F_n)")
legend("topright", inset=c(-0.2,-0.05), 
       bty="n", legend = colnames(sub.states), 
        col = my_palette, lty= 0, pch = 16)
dev.off()

# calc mean coord for each state
f.means<-data.frame()
for (i in 1:ncol(sub.states)){
  print(i)
  #i=1
  st.ids<-which(sub.states[,i]>0)
  sub.flux<-fluxes[st.ids,]
  f.means<-rbind(f.means, colMeans(sub.flux))
}

# compute PCA
pca<-prcomp(f.means)
#str(expr)
summary(pca) # check - 1,2 PC -15% variance 

# plot PCA
tmp<-as.data.frame(pca$x[,1:2])             # get first two comp
#tmp$group<-gsub( ".*-", "", subs$Sample)    # group by normal/not
tmp$names<-colnames(sub.states) # if we want to display sample ids..

# colours by tissue type
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(ncol(sub.states))
tmp$color<-my_palette

cairo_pdf('graphs/pca_means_3x3.pdf', width = 7, height = 7) # save plot
#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # change graph params for legend 
plot(x=tmp$PC1, y=tmp$PC2, col=tmp$color, pch=16, cex=1, 
     xlab="PC1", ylab="PC2") #, label=tmp$names,
text(x=tmp$PC1, y=tmp$PC2, labels = tmp$names, pos = 1)
# legend("topright", inset=c(-0.2,-0.05), 
#        bty="n", legend = colnames(sub.states), 
#        col = my_palette, lty= 0, pch = 16)
abline(-50, -1)
dev.off()

### create hetmap in <Fc> <Fn> for a given state
st.ids<-which(sub.states[,13]>0)
sub.flux<-fluxes[st.ids,]
fc<-rowMeans(sub.flux[,1:3])
fn<-rowMeans(sub.flux[,4:6])
f.means<-as.data.frame(cbind(fc,fn))
br <- seq(1,481,by=4.8)
freq<-matrix(0, length(br)-1, length(br)-1)
#f.means$fc
for (i in 2:length(br)){
  print(i)
  #freq[i,1]<-length(which(f.means$fc<br[i]&f.means$fc>=br[i-1]))
  for (j in 2:length(br)){
    freq[i-1,j-1]<-length(which((f.means$fc<br[i]&f.means$fc>=br[i-1])&(f.means$fn<br[j]&f.means$fn>=br[j-1])))#,]
  }
}
sum(rowSums(freq)) # check that all states are in some bins

# save freqs
write.table(freq, "output/freq_bins_state_13.txt", sep='\t', row.names = T)

library(gplots)
my_palette <- colorRampPalette(rev(brewer.pal(11,"PRGn")))(100)
cairo_pdf('graphs/Heatmap_state_13.pdf', width = 7, height = 7) # save plot

heatmap.2(freq, trace="none", col=my_palette, Rowv=FALSE, Colv=FALSE,
          dendrogram="none", #margins =c(4,4), 
          symm=T, symkey=F, symbreaks=F, srtCol=360,
          colsep=c(0,ncol(freq)), rowsep=c(0,nrow(freq)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") #dendrogram="column",  key.par=list(mar=c(12,12,14,14)), lmat=rbind( c(0, 3, 4), c(2,1,0 ) ), lwid=c(1.5, 4, 2 )
dev.off()

nrow(f.means)
# tmp$color<-""
# k<-0
# for(i in sort(unique(meta$V2))){
#   k<-k+1
#   tmp$color[which(tmp$type==i)]<-my_palette[k]
# }

