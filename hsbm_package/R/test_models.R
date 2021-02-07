# `trans_prob` is the Markov transition probability:
# 1 corresponds to completely random labels for each layer,
# 0 corresponds to the same labels for all layers
#' @export
sample_personality_net = function(n, nlayers, trans_prob = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  Ktru = 3
  eta = matrix(c(0.9, 0.75, 0.5, 0.75, 0.60, 0.25, 0.5, 0.25, 0.10), Ktru, Ktru)
  pri = c(0.45,0.35,0.25)

  if (trans_prob == 1) {
    zb = lapply(1:nlayers, function(t) sample(Ktru, n, replace=T, prob=pri))
  } else if (trans_prob == 0) {
    zb = rep(list(sample(Ktru, n, replace=T, prob=pri)), nlayers)
  } else {
    zinit = sample(Ktru, n, T)
    zb =  sample_markov_labels(zinit, nlayers, trans_prob, Ktru)
  }

  A = lapply(1:nlayers, function(t) nett::fast_sbm(zb[[t]], eta))
  list(A = A, zb = zb, eta = eta, pri = pri)
}


markovian_label <- function(z, trans_prob, K) {
  if (trans_prob == 0) {
    return(z)
  }
  n = length(z)
  # ntrans = rbinom(1, n, K*trans_prob/(K-1))
  ntrans = rbinom(1, n, trans_prob)
  idx_trans = sample(n, ntrans, F)
  z[idx_trans] = sample(K, ntrans, T)
  z
}

#' @export
sample_markov_labels <- function(zinit, nlayers, trans_prob = 0.5, K) {
  z = list()
  z[[1]] <- zinit
  for (t in 2:nlayers) {
    z[[t]] <- markovian_label(z[[t-1]], trans_prob, K)
  }
  z
}

#' @export
test_markov_labels = function(z) {
  mean(sapply(seq_along(z)[-1], function(t) mean(z[[t-1]] == z[[t]])))
}

