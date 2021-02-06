
# http://moreno.ss.uci.edu/data.html#krebs
library(igraph)
library(Matrix)
library(dplyr)
library(multiplex)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../posterior_analysis.R",chdir=T)
source("../diagnostics.R",chdir=T)
source("../hsbm_infer_module.R",chdir=T)


PLOT <- T
RE_ESTIMATE <- F
ITRmax <- 10000
netname <- 'kerbs'
result_filename <- paste(netname,'_',ITRmax,".RData",sep='')
fname <- 'krebs.dat'

Adj <- read.dl(file = fname)
Nlayers <- dim(Adj)[3]

groups <- rep(1:5, c(7,12,17,17,3))
# dim(na.omit(Adj[,,5]))

A <- list()
G <- list()
for (j in 1:Nlayers) {
  # E <- as.matrix( data %>% filter(layer == j) %>% select(e1,e2) )
  # E <- as.matrix( data %>% filter(layer == j) %>% filter(is.element(e1,idx) | is.element(e2,idx)) %>% select(e1,e2) )
  G[[j]] <- simplify( as.undirected(graph_from_adjacency_matrix(Adj[,,j])) )
  A[[j]] <- as_adj(G[[j]])
  # g <- as.undirected( graph_from_edgelist(E, directed = T) )
  # g <- delete_vertices(g, degree(g) == 0)
}

if (PLOT) {
  par(mfrow=c(1,Nlayers), mar=rep(0, 4))
  for (j in 1:Nlayers) {
    plotG(G[[j]])
    #if (j == 1)   out <- plotG(G[[j]]) else     plotG(G[[j]], coord = out$coord)
  }
}


# shapes() # all shapes
source("../additional_vertex_shapes.R", chdir = T)
shapes <- c("sphere","circle","square","triangle","star")

set.seed(2400) # this perhaps has no effect on Rcpp code
if (RE_ESTIMATE) {
  zh_raw <- hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb
  # out <- hsbm_infer(A, beta0=5, gam0=5, ITRmax=ITRmax, Kcap=25, Gcap=25)
  
  save(A,zh_raw,ITRmax,file=result_filename)
} else {
  load(result_filename)  
}

#zh <- compMAPLabel(zh_raw, burnin = burnin, consecutive = T)
out <- compMAPLabel(zh_raw, burnin=floor(ITRmax/2), consecutive = T)
zMAP <- out$MAPlabs
conf <- out$conf

numCom <- max(unlist(zMAP))

layer_names <- 1:Nlayers
for (j in 1:Nlayers) {
  pdf(paste(layer_names[j],'.pdf',sep=''))
  # fixed layout of layer 1
  if (j==1) {
    out <- plotG(A[[j]], zMAP[[j]], numCom, confidence = conf[[j]], shape_group = groups, shapes = shapes)
    coord_1 <- out$coord
  } else {
    out <- plotG(A[[j]], zMAP[[j]], numCom, coord=coord_1, confidence = conf[[j]], shape_group = groups, shapes = shapes)
  }
  
  dev.off()
  if (j > 1) {
    # free layout
    pdf(paste(layer_names[j],'_free.pdf',sep=''))
    out <- plotG(A[[j]], zMAP[[j]], numCom, confidence = conf[[j]], shape_group = groups, shapes = shapes)
    dev.off()
  }
}

par(mfrow=c(1,1))
degs <- sapply(1:Nlayers, function(i) degree(G[[i]]))
colnames(degs) <- layer_names
boxplot(degs)




# Estimate connectivity matrices
 B <- lapply(1:Nlayers, function(t) round(Matrix(estimConnMatrix(A[[t]], zMAP[[t]])),2))
# B
for (j in 1:Nlayers){
  pdf(paste('B',j,'.pdf',sep=''))
  print( image(B[[j]],xlab="",ylab="", sub="",axes=F, yaxt="n", xaxt="n") )
  dev.off()
}



