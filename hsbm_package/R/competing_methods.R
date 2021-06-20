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


# spec_proj = function(A, K, which = "LM", ...) {
#   # A should a sparse matrix
#   U = RSpectra::eigs_sym(A, K, which, ...)$vectors
#   U %*% t(U)
# }

#' @export
pisces = function(A, K, alpha = 0.1, niter = 50, tol = 1e-6,
                  verb = T, shared_kmeans_init = F) {
  nlayers = length(A)
  V0 = lapply(1:nlayers, function(r) nett::spec_repr(A[[r]], K))
  U = c(0, lapply(1:nlayers, function(r) V0[[r]] %*% t(V0[[r]])), 0)
  P = U
  Unew = vector("list", nlayers+2)
  V = vector("list", nlayers)
  for (itr in 1:niter) {
    for (t in 2:(nlayers+1)) {
      # Unew[[t]] = spec_proj(P[[t]] + alpha*(U[[t-1]] + U[[t+1]]), K)
      temp = P[[t]] + alpha*(U[[t-1]] + U[[t+1]])
      V[[t-1]] = RSpectra::eigs_sym(temp, K, which = "LM")$vectors
      Unew[[t]] = V[[t-1]] %*% t(V[[t-1]])
    }
    Unew[[1]] = Unew[[nlayers+2]] = 0
    err = max( sapply(2:(nlayers+1), function(t) norm(U[[t]] - Unew[[t]])) )
    U = Unew
    if (verb) nett::printf("%3.2e\n", err)
    if (err < tol) break
  }
  if (shared_kmeans_init) {
    zinit = sample(nrow(A[[1]]), K, F)
    zh = lapply(1:nlayers, function(t) kmeans(V[[t]], K, centers = V[[t]][zinit, ], iter.max = 30)$cluster)
  } else {
    zh = lapply(1:nlayers, function(t) kmeans(V[[t]], K, nstart = 20)$cluster)
  }
  zh
}

# get_agg_nmi(pisces(A, 3, shared_kmeans_init = F),zb)

