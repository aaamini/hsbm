#' @export
sample_personality_net = function(n, nlayers, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  Ktru = 3  
  eta = matrix(c(0.9, 0.75, 0.5, 0.75, 0.60, 0.25, 0.5, 0.25, 0.10), Ktru, Ktru)
  pri = c(0.45,0.35,0.25)
  zb = lapply(1:nlayers, function(t) sample(1:Ktru, n, replace=T, prob=pri))
  A = lapply(1:nlayers, function(t) fast_sbm(zb[[t]], eta))
  list(A = A, zb = zb, eta = eta, pri = pri)
}
