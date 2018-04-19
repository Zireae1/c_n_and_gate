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
st.ids<-which(sub.states[,19]>0)
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
write.table(freq, "output/freq_bins_state_19.txt", sep='\t', row.names = T)

library(gplots)
my_palette <- colorRampPalette(rev(brewer.pal(11,"PRGn")))(100)
cairo_pdf('graphs/Heatmap_state_19.pdf', width = 7, height = 7) # save plot

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

# check star diagram
st.ids<-which(sub.states[,19]>0)
sub.flux<-fluxes[st.ids,]

for (i in 1:ncol(sub.flux)){
  print(i)
  i=1
  #freq[i,1]<-length(which(f.means$fc<br[i]&f.means$fc>=br[i-1]))
  colMaxs
  apply(sub.flux, 2, max)
  apply(sub.flux, 2, min)
}

length(which(sub.flux[,1]==1))
length(which(sub.flux[,2]==1))
length(which(sub.flux[,6]==1))
sub.flux[which(sub.flux[,2]==1),]

colnames(sub.flux)<-c("C1", "C2", "C3", "N1", "N2", "N3")
col.lines<-alpha(rep(my_palette[19], nrow(sub.flux)), 0.05)
#test<-t(sub.flux)
# stars(sub.flux, locations = c(0, 0), radius  =  FALSE, 
#       key.loc = c(0, 0), main = "Test", #col.lines = my_palette[19], 
#       lwd=0.05)


### star plot for mean coords
colnames(f.means)<-c("C1", "C2", "C3", "N1", "N2", "N3")
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(ncol(sub.states))

palette(rainbow(12, s = 0.6, v = 0.75))
#max(f.means)
#min(f.means)

f.means<-f.means/max(f.means)
cairo_pdf('graphs/Stars_for_states_v1.pdf', width = 10, height = 10) # save plot

# plot in <Fc> <Fn> centroid space
stars(f.means, len = 6, 
      location = cbind(fc.means,fn.means), axes = TRUE, xlab = "Mean(F_c)", ylab="Mean(F_n)",
      #key.loc = c(13, 1),  
      key.loc = c(290, 330), labels = colnames(sub.states),
      main = "States mean fluxes", draw.segments = T, #F
      frame.plot = TRUE, nrow = 4, cex = 0.8, scale = T)
dev.off()

f.means[11,]

# exp with one state
# calc sd from mean coord for each state
f.sds<-data.frame()
for (i in 1:ncol(sub.states)){
  print(i)
  #i=1
  st.ids<-which(sub.states[,i]>0)
  sub.flux<-fluxes[st.ids,]
  f.sds<-rbind(f.sds, apply(sub.flux, 2, sd))
}

#f.means<-f.means*386.3987
id<-2

cairo_pdf('graphs/Stars_for_states_v2.pdf', width = 10, height = 10) # save plot

par(mfrow=c(4,6))
#par(mfrow=c(1,1))
for(id in 1:22){
  print(id)
  stars(rbind(f.means[id,], (f.means[id,]-f.sds[id,]), (f.means[id,]+f.sds[id,])), #len = 6, 
        #location = cbind(fc.means,fn.means), axes = TRUE, xlab = "Mean(F_c)", ylab="Mean(F_n)",
        #locations = c(100, 100),
        #key.loc = c(0, 0), #labels = colnames(sub.states),
        main = paste("S", id, sep=""),
        draw.segments = T, scale = F, 
        col.lines = c("red", "blue", "blue"), lwd=2)
  
}
dev.off()

#######    C1N1   C1N2  C1N3  C2N1  C2N2  C2N3  C3N1  C3N2  C3N3
lambda.c<-c(88.6,	12.0,	10.6,	30.5,	66.5,	85.5,	78.4,	46.4,	82.5)	
lambda.n<-c(33.6,	82.9,	26.8,	94.2,	44.2,	90.3,	24.4,	78.0,	18.3)

bac.names<-c("C1N1", "C1N2", "C1N3", "C2N1", "C2N2", "C2N3", "C3N1", "C3N2", "C3N3")
rank.c<-c(1, 2, 3, 3, 2, 1, 2, 3, 1)
rank.n<-c(2, 1, 2, 1, 3, 1, 3, 2, 3)

bac.ranks<-rbind(rank.c, rank.n)
colnames(bac_ranks)<-bac.names

species<-fread("data_states/FinalStats_1_3_3_100000.txt", header = F, stringsAsFactors = F)
species<-species[,1:9]
colnames(species)<-bac.names
### wrong way of calculation
# rc.mean<-vector()
# rn.mean<-vector()
# for (i in 1:nrow(species)){
#   #i=1
#   print(i)
#   rc.mean[i]<-mean(bac.ranks[1,which(species[i,]!=0)])
#   rn.mean[i]<-mean(bac.ranks[2,which(species[i,]!=0)])
# }

tmp.ranks.c<-t(matrix(rank.c, 3, 3))
#colnames(tmp.ranks.c)<-c("N1", "N2", "N3")
#rownames(tmp.ranks.c)<-c("C1", "C2", "C3")
tmp.ranks.n<-t(matrix(rank.n, 3, 3))
#colnames(tmp.ranks.n)<-c("N1", "N2", "N3")
#rownames(tmp.ranks.n)<-c("C1", "C2", "C3")
rc.mean<-vector()
rn.mean<-vector()
for (s in 1:nrow(species)){
  #i=1
  print(s)
  tmp<-t(matrix(species[s,], 3, 3))
  #colnames(tmp)<-c("N1", "N2", "N3")
  #rownames(tmp)<-c("C1", "C2", "C3")
  # calculate carbon ranks properly
  sum.rc<-0
  for (i in 1:nrow(tmp)){
    #i=2
    if (length(which(tmp[i,]==-1))>0){
      sum.rc<-sum.rc+tmp.ranks.c[i,which(tmp[i,]==-1)]
    }else{
      sum.rc<-sum.rc+4
    }
  }
  rc.mean[s]<-sum.rc/nrow(tmp)
  # calculate nitrogen ranks properly
  sum.rn<-0
  for (i in 1:ncol(tmp)){
    #i=1
    if (length(which(tmp[,i]==1))>0){
      sum.rn<-sum.rn+tmp.ranks.n[which(tmp[,i]==1),i]
    }else{
      sum.rn<-sum.rn+4
    }
  }
  rn.mean[s]<-sum.rn/ncol(tmp)
}


par(mfrow=c(1,1))
plot(rc.mean, rn.mean, xlim=c(1,3), ylim=c(1,3))

cairo_pdf('graphs/Stars_for_states_v3.pdf', width = 10, height = 10) # save plot

#f.means<-f.means/max(f.means) #need this for scale
# plot in <Fc> <Fn> centroid space
stars(f.means, len = 0.13, 
      location = cbind(rc.mean, rn.mean), axes = TRUE, xlim=c(1.2,4), ylim=c(0.9,4), 
      xlab = "Mean(r_c)", ylab="Mean(r_n)",
      key.loc = c(4, 3.8),  
      #key.loc = c(290, 330), 
      labels = colnames(sub.states),
      main = "States mean fluxes", draw.segments = T, #F
      frame.plot = TRUE, nrow = 4, cex = 0.8, scale = T)

dev.off()

cairo_pdf('graphs/Rank_vs_flux_plot.pdf', width = 6, height = 6) # save plot

plot(rc.mean, fc.means, col="royal blue", pch =16, 
     xlab = "Average rank of carbon utilization", 
     ylab = "Average carbon flux", 
     xlim=c(0.9,4))
text(rc.mean, fc.means,labels=colnames(sub.states), pos = 2, cex=0.8)
legend("topleft", paste("cor = ", round(cor(rc.mean, fc.means), digits = 2),sep=""),
       col = c("royal blue"), 
       pch=16, bty="n")

dev.off()

cor(rc.mean, rowMeans(f.means[,1:3]))
cor(rc.mean, fc.means)
plot(rn.mean, fn.means, col="royal blue", pch =16)
cor(rn.mean, fn.means)

plot(rc.mean, fc.means, col="royal blue", pch =16, 
     xlab = "Average rank of carbon utilization", 
     ylab = "Average carbon flux", 
     xlim=c(0.9,2.3))
text(rc.mean, fc.means,labels=colnames(sub.states), pos = 2, cex=0.8)
legend("topleft", paste("cor = ", round(cor(rc.mean, fc.means), digits = 2),sep=""),
       col = c("royal blue"), 
       pch=16, bty="n")




### nice pieces of code
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
x <- rnorm(n)
x4 <- x + .1*rnorm(n)
x5 <- x + .1*rnorm(n)
y <- 1 + x1 + x4 + rnorm(n)
d <- data.frame(y,x1,x2,x3,x4,x5)
library(ellipse)
plotcorr(cor(d), main = "plotcorr (in the \"ellipse\" package)")
