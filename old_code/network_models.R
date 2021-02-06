################## Data-set example simulation###############
sample_sbm2 <- function(n, eta, z) {
  tab <- table(z)
  idx <- as.integer(names(tab))
  
  sample_sbm(n, pref.matrix = eta[idx,idx], block.sizes = tab)
}
#############################################################

markovian_label <- function(z, trans_prob = 0.5, K) {
  # all_labels <- 1:K
  # sapply(1:length(z), function(i) ifelse( runif(1) < p, z[i], sample(setdiff(all_labels,z[i])) ))
  n = length(z)
  ifelse( runif(n) <= 1-K*trans_prob/(K-1), z, sample(1:K, n, T))  
}

markov_label_process <- function(zinit, Tn, p, K) {
  z = list()
  z[[1]] <- zinit
  for (t in 2:Tn) {
    z[[t]] <- markovian_label(z[[t-1]], p, K)
  }
  z
}


nlayers = 1000
z = sample(1:3, nlayers, T)
zz = markov_label_process(z, nlayers, .2, 3)
mean(sapply(1:(nlayers-1), function(t) mean(zz[[t]] == zz[[t+1]])))

microbenchmark::microbenchmark(zz = markov_label_process(z, 100, .2, 3))           
