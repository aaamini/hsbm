# here zh has the structure zh[[itr]][[j]]: the labels for layer "j" in iteration "itr"
# source("network_commons.R")

Modes <- function(x) {
  # find the modes of a vector (the values with largest frequency)
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}



############### Computes the MAP estimated of the labels ##################
#' @export 
get_map_labels <- function(zh, burnin=NULL, consecutive=TRUE){
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
  
  list(labels=d, conf=zhConf, Zh=Zh)
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

