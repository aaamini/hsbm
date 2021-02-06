
#' @export
spec_clust_sliced = function(A, K, tau = 0.25) {
  lapply(seq_along(A), function(t) nett::spec_clust(A[[t]], K, tau=tau))
}

#' @export
spec_clust_avg = function(A, K, tau = 0.25) {
  nlayers = length(A)
  Abar = Reduce(`+`,A)/nlayers
  zh = nett::spec_clust(Abar, K, tau=tau)
  rep(list(zh), nlayers) # repeat the same label vector for all layers
}

