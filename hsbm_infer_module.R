library(Rcpp)
Rcpp::sourceCpp('hsbm_inference.cpp')

hsbm_infer <- function(A, 
                       beta0=3,  gam0=1,  alpha_eta = 1, beta_eta = 1,
                       ITRmax=50, Kcap=20, Gcap=20) {
  
  res <- hsbm_inferC(A, 
                     beta0=beta0, gam0=gam0,
                     alpha_eta, beta_eta,
                     ITRmax=ITRmax, Kcap=Kcap, Gcap=Gcap)
  n <- sapply(A, function(At) dim(At)[1])
  
  Tn <- length(A)
  # only 1:n[t] entries in each zb[t,:] is valid
  zb <- lapply(1:ITRmax, function(itr) lapply(1:Tn, function(t) res$zb_list[[itr]][t,1:n[t]]))
  times <- res$times
  list(zb=zb, times=times)
}


dpsbm_infer <- function(A, 
                        gam0=1, alpha_eta = 1, beta_eta = 1, 
                        ITRmax=50, Zcap=20) {
  # A is a single adjacency matrix
  res <- dpsbm_inferC(A, 
              gam0=gam0, alpha_eta, beta_eta,
              ITRmax=ITRmax, Zcap=Zcap)
  res$zb_list
}

dpsbm_slice_infer <- function(A, 
                              gam0=1, alpha_eta = 1, beta_eta = 1, 
                              ITRmax=50, Zcap=20) {
  # A is a list of adjacency matrices
  Tn <- length(A)
  zh_raw_rev <- lapply(1:Tn, function(t) dpsbm_infer(A[[t]], 
                                                      gam0=gam0, alpha_eta, beta_eta, 
                                                      ITRmax=ITRmax, Zcap=Zcap) )
  
  # turn from indexing [[t]][[itr]] to [[itr]][[t]]
  zh_raw <- lapply(1:ITRmax, function(itr) lapply(1:Tn, function(t) zh_raw_rev[[t]][[itr]] ))
  
  zh_raw
}