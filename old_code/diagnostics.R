# here zh has the structure zh[[itr]][[j]]: the labels for layer "j" in iteration "itr"
source("posterior_analysis.R")

compute_nmi_corr <- function(zh, lag_vec=1:50, blk_size=1000) {
  ITRmax <- length(zh)
  big_Zh <- sapply(1:ITRmax, function(itr) unlist(zh[[itr]]))

  begin_idx <- c(1,        ITRmax-blk_size)
  end_idx   <- c(blk_size, ITRmax)
  R <- length(begin_idx)
  nmi_corr <- matrix(0, nrow=length(lag_vec), ncol=R)
  for (l in 1:length(lag_vec)) {
    for (r in 1:R) {
      lag <- lag_vec[l]
      
      nmi_corr[l,r] <- compute_nmi(c(big_Zh[,begin_idx[r]:(end_idx[r]-lag)]), c(big_Zh[,(begin_idx[r]+lag):end_idx[r]]))
    }
  }
  
  # temp <- sapply(1:(ITRmax-lag), function(t) compute_nmi(big_Zh[,t],big_Zh[,t+lag])) 
  list(nmi_corr=nmi_corr, lag_vec=lag_vec)
}
  
compute_sequential_nmi <- function(zh, lag=1) {
  ITRmax <- length(zh)
  
  sapply( (1+lag):ITRmax, function(t) compute_nmi( unlist(zh[[t]]), unlist(zh[[t-lag]]) ) )
}
