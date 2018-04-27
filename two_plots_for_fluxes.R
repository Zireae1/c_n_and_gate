rm(list=ls(all=TRUE))   # clean environment
getwd()
#update.packages(ask = FALSE, dependencies = c('Suggests')) #update all installed R packages

library(data.table)     # for fread function
library(RColorBrewer)   # for nice color palettes
#library(scales)         # for points transparency on plots
library(gplots)         # for heatmap.2 function

setwd("/Users/Zireael/Desktop/Maslov/CommunMetab") # replace with your working directory
#source("functions.R")    # some custom functions

### load data (use binary format for states)
# pick correct example:
# load from a single file:
states<-fread("data_states/FinalStats_1_3_3_full_new_01.dat", header = F, stringsAsFactors = F)

# load data from several files in the given directory and combine them:
dir<-"data_states/Data-MonteCarlo-PhaseVolume-April20"
files<-list.files(dir)
states<-data.frame()
for (f in files){
  states<-rbind(states, fread(paste(dir, "/", f, sep=""), header = F, stringsAsFactors = F))
}

# look at the loaded table:
states[1:10,]

# separate flux and states tables:
fluxes<-states[,1:6]
sub.states<-states[, 7:ncol(states)]

# column names for states
colnames(sub.states)<-c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09",
                        "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18",
                        "S19", "S20", "S21", "S22")
sub.states<-as.data.frame(sub.states) # this data type conversion is crucial

# calculate centroids for all states:
fc.means<-vector()
fn.means<-vector()
for (i in 1:ncol(sub.states)){
  print(i)
  #i=1
  st.ids<-which(sub.states[,i]==1)
  sub.flux<-fluxes[st.ids,]
  fc<-rowMeans(sub.flux[,1:3])
  fn<-rowMeans(sub.flux[,4:6])
  fc.means[i]<-mean(fc)
  fn.means[i]<-mean(fn)
}

### save plot for Centroids of all states in <Fc> <Fn>
cairo_pdf('State_centroids.pdf', width = 7, height = 6) # save plot
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE) # change graph params for legend 
my_palette <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(fc.means))
plot(fc.means, fn.means, col = my_palette, pch=16, xlim=c(0,1000), ylim=c(0,1000), 
     cex=1.3, xlab = "Mean(F_c)", ylab="Mean(F_n)")
legend("topright", inset=c(-0.2,-0.05), 
       bty="n", legend = colnames(sub.states), 
       col = my_palette, lty= 0, pch = 16)
dev.off()

### create heatmap in <Fc> <Fn> for a given state
id<-2 # pick a state
st.ids<-which(sub.states[,id]>0) # find all corresponding rows for the state 
sub.flux<-fluxes[st.ids,]        # create subset of fluxes for the state
fc<-rowMeans(sub.flux[,1:3])     # calc average fluxes
fn<-rowMeans(sub.flux[,4:6])

br <- seq(1,max(max(c(fc, fn))),by=max(max(c(fc, fn)))/100)  # create breaks for flux intervals on heatmap
freq<-matrix(0, length(br)-1, length(br)-1)    # initalize empty matrix for counts
for (i in 2:length(br)){                       # count how many times we have this state for a given range of fluxes
  print(i)
  for (j in 2:length(br)){
    freq[i-1,j-1]<-length(which((fc<br[i]&fc>=br[i-1])&(fn<br[j]&fn>=br[j-1])))
  }
}
sum(rowSums(freq)) # check that all instances are in some bins (there is a small difference due to the last interval)

# save plot
my_palette <- colorRampPalette(rev(brewer.pal(11,"PRGn")))(100)
cairo_pdf(paste("Heatmap_for_state_", id, ".pdf", sep=""), width = 7, height = 7) 
heatmap.2(freq, trace="none", col=my_palette, Rowv=FALSE, Colv=FALSE,
          dendrogram="none", #margins =c(4,4), 
          symm=T, symkey=F, symbreaks=F, srtCol=360,
          colsep=c(0,ncol(freq)), rowsep=c(0,nrow(freq)), 
          sepwidth=c(0.01, 0.01), sepcolor="#e0e0e0") 
dev.off()
