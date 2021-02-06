library(igraph)
library(Matrix)


# set working directory to current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../posterior_analysis.R",chdir=T)
source("../diagnostics.R",chdir=T)
source("../hsbm_infer_module.R",chdir=T)

##############
fnames = lapply(1:3, function(i) paste("data/A",i,".txt",sep=""))
A <- list()
for (j in 1:3) {
  A[[j]] <- Matrix( as.matrix(read.table(fnames[[j]])) )
}
gender <- read.csv('data/gender.txt',header=T,sep="")$x+1
ITRmax <- 10000
fname <- paste("GlasgowEstim_",ITRmax,'.RData',sep='')

RE_ESTIMATE <- FALSE
DIAGNOSTICS <- FALSE

#set.seed(1400) # this perhaps has no effect on Rcpp code

set.seed(2400) # this perhaps has no effect on Rcpp code

if (RE_ESTIMATE) {
  
  out <- hsbm_infer(A, beta0=5, gam0=5, ITRmax=ITRmax, Kcap=25, Gcap=25)
  zh <- out$zb
  diff(out$times)/ITRmax
  save(A,zh,ITRmax,file=paste("GlasgowEstim.RData",ITRmax,sep='_') )
} else {
  load(fname)  
}



out <- compMAPLabel(zh, burnin=5000, consecutive = T)
d <- out$MAPlabs
conf <- out$conf

numCom <- max(unlist(d))

numCom
d

png('conf_hist.png', width = 9, height=3, units = 'in', res = 150)
par(mfrow=c(1,3))
for (j in 1:3)   hist(conf[[j]],15, xlab="confindence",main = "") 
dev.off()


load('data/coords.RData')
par(mfrow=c(1,3), mar=rep(0, 4))
for (j in 1:3)   plotG(A[[j]], d[[j]], numCom, coord=coord[[j]], shape_group = gender)


png('conf_graph.png', width = 9, height=3, units = 'in', res = 150)
par(mfrow=c(1,3), mar=rep(0, 4))
for (j in 1:3) plotG(A[[j]], d[[j]], numCom, coord=coord[[j]], confidence = conf[[j]], shape_group = gender)
dev.off()

##### Diagnostics
if (DIAGNOSTICS) {
  library(plotly)
  seq_nmi <- compute_sequential_nmi(zh)
  plot_ly(x=2:ITRmax, y=seq_nmi,type="scatter",mode="markers",name="seq. NMI")
  # out <- compute_nmi_corr(zh)
  #  plot_ly(x=out$lag_vec, y=out$nmi_corr[,1],type="scatter",mode="markers",name="beginning")   %>% 
  #    add_trace(y=out$nmi_corr[,2],type="scatter", mode="markers", name="end")%>%  
  #    layout(yaxis = list(title = "NMI correlation"), 
  #           xaxis = list (title = "lag"))
}
##
# load("data/coords.RData")
layer_names <- c('1995','1996','1997')
coord <- list()
for (j in 1:3) {
  pdf(paste(layer_names[j],'_bygender.pdf',sep=''))
  out <- plotG(A[[j]], c(2,4)[gender], 5, shape_group = gender)
  #out <- plotG(A[[j]], gender, 2, shape_group = gender)
  coord[[j]] <- out$coord
  dev.off()
  
}

# pdf('1996_bygender.pdf')
# out <- plotG(A[[2]], gender, 2, shape_group = gender)
# coord[[2]] <- out$coord
# dev.off()
# 
# pdf('1997_bygender.pdf')
# out <- plotG(A[[3]], gender, 2, shape_group = gender)
# coord[[3]] <- out$coord
# dev.off()

##### 
#numCom <- max(unlist(d))
#png(filename = '1995.png')
for (j in 1:3) {
  pdf(paste(layer_names[j],'.pdf',sep=''))
  # fixed layout of layer 1
  out <- plotG(A[[j]], d[[j]], numCom, coord=coord[[1]], confidence = conf[[j]], shape_group = gender)
  dev.off()
  if (j > 1) {
    # free layout
    pdf(paste(layer_names[j],'_free.pdf',sep=''))
    out <- plotG(A[[j]], d[[j]], numCom,  coord=coord[[j]], confidence = conf[[j]], shape_group = gender)
    dev.off()
  }
}



