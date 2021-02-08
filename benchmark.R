library(nett)
library(hsbm)

n = 200 # number of nodes in each layer
nlayers = 5 # number of layers

# `trans_prob` is the Markov transition probability: 
# 1 corresponds to completely random labels for each layer, 
# 0 corresponds to the same labels for all layers
out = sample_personality_net(n, nlayers, trans_prob = 1, seed=1400) 
Ktru = nrow(out$eta)
A = out$A
zb = out$zb

niter <- 100 # number of Gibbs iterations 
burnin <- ceiling(niter/2)

tau = 0.0
methods = list()
methods[["HSBM"]] = function(A) {
  zh = fit_hsbm(A, beta0=0.1, gam0=.5, niter=niter, Kcap=10, Gcap=10, seq_g_update = F, verb = F)$zb
  get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["DP-SBM"]] =  function(A) {
   zh = fit_mult_dpsbm(A, gam0=.5, niter=niter, Zcap=10, verb = F)$zb
   get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["SC-sliced"]] = function(A) spec_clust_sliced(A, Ktru, tau = tau)
methods[["SC-avg"]] = function(A) spec_clust_avg(A, Ktru, tau = tau)
methods[["SC-ba"]] = function(A) spec_clust_bias_adj(A, Ktru)
methods[["SC-omni"]] = function(A) spec_clust_omnibus(A, Ktru)

mtd_names = names(methods)

res = do.call(rbind, lapply(seq_along(methods), function(j) {
  dt = as.numeric(system.time( zh <- methods[[j]](A) )["elapsed"])
  data.frame(method = mtd_names[j], 
             aggregate_nmi = get_agg_nmi(zb, zh), 
             slicewise_nmi = get_slice_nmi(zb, zh) , elapsed_time = dt)
}))
  
print( knitr::kable(res, digits = 4, format="pipe") )
