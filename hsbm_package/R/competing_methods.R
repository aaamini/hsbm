
#' @export
spec_clust_sliced = function(A, K, tau = 0.25, ...) {
  lapply(seq_along(A), function(t) nett::spec_clust(A[[t]], K, tau=tau), ...)
}

#' @export
spec_clust_avg = function(A, K, tau = 0.25, ...) {
  nlayers = length(A)
  Abar = Reduce(`+`,A)/nlayers
  zh = nett::spec_clust(Abar, K, tau=tau, ...)
  rep(list(zh), nlayers) # repeat the same label vector for all layers
}

#' @export
spec_clust_bias_adj = function(A, K, nstart = 20, niter = 10) {
  nlayers = length(A)
  A_bias_adj = lapply(seq_along(A), function(t) {
    # A[[t]] %*% A[[t]] - diag(Matrix::rowSums(A[[t]]))
    A[[t]] %*% A[[t]] - Matrix::Diagonal(x=Matrix::rowSums(A[[t]]))
  })
  A2bar_bias_adj = Reduce(`+`, A_bias_adj)/nlayers
  eig_res <- RSpectra::eigs_sym(A2bar_bias_adj, K)
  U <- eig_res$vectors[ , 1:K]
  zh = kmeans(U, K, nstart = nstart, iter.max = niter)$cluster
  rep(list(zh), nlayers)
}

#' @export
spec_clust_omnibus = function(A, K, nstart = 20, niter = 10) {
  layer_seq = seq_along(A)
  Aomni = do.call(rbind, lapply(layer_seq, function(j)
    do.call(cbind, lapply(layer_seq, function(i) (A[[i]] + A[[j]])/2))
  ))
  eig_res <- RSpectra::eigs_sym(Aomni, K)
  U <- eig_res$vectors[ , 1:K]
  # S <- diag(eig_res$values)  # won't work in general due to negative eigenvalues
  zh = kmeans(U, K, nstart = nstart, iter.max = niter)$cluster

  net_sizes = sapply(A, nrow)
  cnet_sizes = c(0, cumsum(net_sizes))
  idx_list = lapply(layer_seq, function(t) cnet_sizes[t] + 1:net_sizes[t])
  lapply(idx_list, function(idx) zh[idx])
}

# zh = spec_clust_bias_adj(A, 3)
# get_agg_nmi(zh, zb)

