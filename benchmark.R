library(nett)
library(hsbm)

n = 200 # number of nodes in each layer
nlayers = 5 # number of layers

out = sample_personality_net(n, nlayers, seed=1400)
Ktru = nrow(out$eta)
A = out$A
zb = out$zb

niter <- 100
burnin <- ceiling(niter/2)

tau = 0.0
methods = list()
methods[["HSBM"]] = function(A) {
  # zh = hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb
  zh = fit_hsbm(A, beta0=0.1, gam0=.5, niter=niter, Kcap=10, Gcap=10, seq_g_update = F, verb = F)$zb
  get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["DP-SBM"]] =  function(A) {
   zh = fit_mult_dpsbm(A, gam0=.5, niter=niter, Zcap=10, verb = F)$zb
   get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["SC"]] = function(A) spec_clust_sliced(A, Ktru, tau = tau)
methods[["SC-Abar"]] = function(A) spec_clust_avg(A, Ktru, tau = tau)

mtd_names = names(methods)

res = do.call(rbind, lapply(1:length(methods), function(j) {
  dt = as.numeric(system.time( zh <- methods[[j]](A) )["elapsed"])
  agg_nmi = get_agg_nmi(zb, zh)
  slice_nmi = get_slice_nmi(zb, zh) # mean(compute_slice_nmi(zb, zh))
  data.frame(method_name = mtd_names[j], aggregate_nmi = agg_nmi, 
             slicewise_nmi = slice_nmi, elapsed_time = dt)
}))
  
print( knitr::kable(res, digits = 4, format="pipe") )
