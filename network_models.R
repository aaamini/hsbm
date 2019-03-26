################## Data-set example simulation###############
sample_sbm2 <- function(n, eta, z) {
  tab <- table(z)
  idx <- as.integer(names(tab))
  
  sample_sbm(n, pref.matrix = eta[idx,idx], block.sizes = tab)
}
#############################################################

markovian_label <- function(z,p,K) {
  all_labels <- 1:K
  sapply(1:length(z), function(i) ifelse( runif(1) < p, z[i], sample(setdiff(all_labels,z[i])) ))  
}

markov_label_process <- function(zinit, Tn, p, K) {
  z = list()
  z[[1]] <- zinit
  for (t in 2:Tn) {
    z[[t]] <- markovian_label(z[[t-1]],p,K)
  }
  z
}
