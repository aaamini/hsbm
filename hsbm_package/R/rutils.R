#' @export 
get_agg_nmi <- function (y, z){
  # This computes aggregate NMI if y and z are lists, i.e., unlists them first
  # otherwise usual NMI between vector
  # require(clue)
  # y_labs <- as.cl_hard_partition( unlist(y) )
  # z_labs <- as.cl_hard_partition( unlist(z) )
  
  # cl_agreement(y_labs, z_labs, method = 'NMI') 
  nett::compute_mutual_info(unlist(y), unlist(z))  # our own NMI since the one from "clue" is buggy
}


#' @export
get_slice_nmi <- function (zlist, ylist){
  # compute mean slicewise nmi
  Tn <- length(zlist)

  if (length(ylist) != Tn)  print("zlist and ylist should be lists of the same length.")
  
  # sapply(1:Tn, function(t) compute_nmi(zlist[[t]], ylist[[t]]) )
  temp = sapply(1:Tn, function(t) nett::compute_mutual_info(zlist[[t]], ylist[[t]]))
  mean(temp)
}


################### Plot the graphs ########################
# assuming labels are 1-based and consecutive
#' @export
plotG <- function(A, d=NULL, MAX_LABEL=NULL, confidence=NULL,
                  shape_group=NULL, shapes = c("square","circle"), 
                  coord=NULL, deg_func = function(deg) log(deg+3)*3, ...){
  # A : adj matrix
  # d : 0-based labels
  # MAX_LABEL: maximum label
  
  if (!is(A, 'igraph')) {
    if (is(A, 'sparseMatrix')) adj = A else adj = Matrix(A)  
    g <- graph_from_adjacency_matrix(adj, mode = "undirected")
    g = simplify(g)
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

