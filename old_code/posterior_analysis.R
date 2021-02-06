# here zh has the structure zh[[itr]][[j]]: the labels for layer "j" in iteration "itr"
source("network_commons.R")

Modes <- function(x) {
  # find the modes of a vector (the values with largest frequency)
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

#######  Compute NMI functions #######
compute_nmi <- function (y, z){
  # This computes aggregate NMI if y and z are lists, i.e., unlists them first
  # otherwise usual NMI between vector
  # require(clue)
  # y_labs <- as.cl_hard_partition( unlist(y) )
  # z_labs <- as.cl_hard_partition( unlist(z) )
  
  # cl_agreement(y_labs, z_labs, method = 'NMI') 
  compute_mutual_info(unlist(y), unlist(z))  # our own NMI since the one from "clue" is buggy
}

# compute_agg_nmi <- function (zlist, ylist){
#   # compute aggregate nmi
#   # zlist and ylist are two lists of label vector which will be unlisted before computing their NMI
#   compute_nmi(unlist(zlist), unlist(ylist))
# }

compute_slice_nmi <- function (zlist, ylist){
  # compute mean slicewise nmi
  Tn <- length(zlist)

  if (length(ylist) != Tn)  print("zlist and ylist should be lists of the same length.")
  
  
  #nmi <- numeric(Tn)
  #for (t in 1:Tn)   nmi[t] <- compute_nmi(zlist[[t]], ylist[[t]])
  #mean(nmi)
  sapply(1:Tn, function(t) compute_nmi(zlist[[t]], ylist[[t]]) )
}


############### Computes the MAP estimated of the labels ##################
compMAPLabel <- function(zh, burnin=NULL, consecutive=TRUE){
  ITRmax <- length(zh)
  nlayers <- length(zh[[1]])
  nn <- sapply(1:nlayers, function(j) length(zh[[1]][[j]]) )
  
  if ( is.null(burnin) ) {
    burnin <- round(ITRmax/2)
  }
  
  d <- list()
  zhConf <- list()
  for (j in 1:nlayers){
    nj <- nn[j]
    # Zh <- sapply((burnin+1):ITRmax, function(itr) zh[[itr]][[j]])
    Zh <-  do.call(cbind, lapply((burnin+1):ITRmax, function(itr) zh[[itr]][[j]])) 
    zhMAP <- sapply(1:nj, function(i)  Modes( Zh[i,] )[1] )
    
    d[[j]] <- zhMAP
    zhConf[[j]] <- sapply(1:nj, function(i) sum(Zh[i,]==zhMAP[i])/length(Zh[i,]) )
    
  }
  
  if (consecutive) {
    # make community labels consecutive from 1 to "number of communities"
    dflat <- unlist(d)
    dflat_new <- dflat
    comm_labs <- sort(unique(dflat))
    for (i in 1:length(comm_labs)){
      dflat_new[dflat==comm_labs[i]] <- i
    }
    d <- split(dflat_new,  unlist(lapply(1:nlayers, function(i) rep(i,nn[i]))))
  }
  
  list(MAPlabs=d, conf=zhConf, Zh=Zh)
}


########### Cumulative MAP
compCommMAP <- function(zh, burnin=NULL, burnout=NULL ) {
  # zh is a list containing estimated labels (themselves lists) at each iteration
  ITRmax <- length(zh)
  
  if (is.null(burnin)) burnin <- ceiling(ITRmax/5)
  if (is.null(burnout)) burnout <- floor(4*ITRmax/5)
  
 
  zhMAP_list <- lapply(burnin:burnout, function(bin) compMAPLabel(zh_raw, burnin = bin)$MAPlabs )
  consec_nmi <- sapply(2:length(zhMAP_list), function(t) compute_nmi(zhMAP_list[[t-1]], zhMAP_list[[t]]))
 
 return(list(zhMAPs = zhMAP_list, consec_nmi = consec_nmi)) 
}

################### Plot the graphs ########################
# assuming labels are 1-based and consecutive
plotG <- function(A, d=NULL, MAX_LABEL=NULL, confidence=NULL,
                  shape_group=NULL, shapes = c("square","circle"), 
                  coord=NULL, deg_func = function(deg) log(deg+3)*3, ...){
  # A : adj matrix
  # d : 0-based labels
  # MAX_LABEL: maximum label
  
  if (!is(A, 'igraph')) {
    if (is(A, 'sparseMatrix')) adj = A else adj = Matrix(A)  
    g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  } else{
      g <- A # A is in fact g
      adj <- get.adjacency(g)
  }
  
  
  if (is.null(d)) d <- rep(1,vcount(g))
  V(g)$community <- d
  
  if (is.null(MAX_LABEL)) MAX_LABEL <- max(d)
  
  if (is.null(coord)) {  
    #set.seed(1400)
    coord <- layout_with_fr(g, niter = 1000)  
  }
  
  if (!is.null(shape_group)) {  V(g)$shape <- shapes[shape_group]  }
  
  colfunc<-colorRampPalette(c("red","yellow","springgreen","royalblue"))
  colors <- colfunc(MAX_LABEL)
  
  vcolors <- colors[V(g)$community]  
  if ( !is.null(confidence) ) {
    vcolors <- sapply(1:length(vcolors), function(i) adjustcolor(vcolors[i], alpha.f = confidence[i]))
  }
  
  # plot(rep(1,5),col=colors, pch=19,cex=2)
  #vsize <- log(degree(g)+3)*3
  vsize <- deg_func(degree(g))
  plot(g, 
       layout=coord, 
       vertex.color = vcolors,
       vertex.label =NA, 
       vertex.size=vsize, ... )
  
  list(coord=coord,graph=g,adj=adj)
}

estimConnMatrix <- function(A,z,issym=T) {
  K <- max(z)
  nc <- tabulate(z)
  B <- matrix(0,nrow=K,ncol=K)
  for (k in 1:K) {
    #for (l in 1:K) {
    for (l in k:K) {
      if (k==l) nn <-nc[k]*(nc[k]-1) else nn <- nc[k]*nc[l]
      B[k,l] <- B[l,k] <- sum(A[z==k,z==l]) / nn
    }
  }
  B
}

