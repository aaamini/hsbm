library(nett)
library(hsbm)

ncores = max(parallel::detectCores() - 1, 20) # number of cores to use to parallel
nreps = 25 # number of replications
ntrans = 8
seq_g_update = F
n = 200 # number of nodes in each layer
nlayers = 5 # number of layers

niter <- 100 # number of Gibbs iterations 
burnin <- ceiling(niter/2)

tau = 0.0
methods = list()
methods[["HSBM"]] = function(A, K) {
  # K is not used
  zh = fit_hsbm(A, beta0=0.1, gam0=.5, niter=niter, Kcap=10, Gcap=10, 
                seq_g_update = seq_g_update, verb = F)$zb
  get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["DP-SBM"]] =  function(A, K) {
   # K is not used
   zh = fit_mult_dpsbm(A, gam0=.5, niter=niter, Zcap=10, verb = F)$zb
   get_map_labels(zh, burnin = burnin, consecutive = T)$labels
}
methods[["SC-sliced"]] = function(A, K) spec_clust_sliced(A, K, tau = tau)
methods[["SC-avg"]] = function(A, K) spec_clust_avg(A, K, tau = tau)

mtd_names = names(methods)


runs = expand.grid(
  mtd_idx = seq_along(methods), 
  trans_prob = seq(0, 1, length.out = ntrans),
  rep = 1:nreps  
)

total_time = system.time(
  res <- do.call(rbind, parallel::mclapply(1:nrow(runs), function(j) {
  # res = do.call(rbind, lapply(1:nrow(runs), function(j) {
    mi = runs[j, "mtd_idx"]
    trans_prob = runs[j, "trans_prob"]
    out = sample_personality_net(n, nlayers, trans_prob = trans_prob) # , seed=1400) 
    Ktru = nrow(out$eta)
    A = out$A
    zb = out$zb
    
    dt = as.numeric(system.time( zh <- methods[[mi]](A, Ktru) )["elapsed"])
    data.frame(method = mtd_names[mi], 
               aggregate_nmi = get_agg_nmi(zb, zh), 
               slicewise_nmi = get_slice_nmi(zb, zh) , 
               elapsed_time = dt, trans_prob = trans_prob)
  # }))
  }, mc.cores = ncores))
)["elapsed"]
nett::printf("Total simulation time = %3.2f (s)\n" , total_time)

library(ggplot2)
library(dplyr)
res %>% 
  group_by(method, trans_prob) %>% 
  summarise_all(mean) %>% 
 # mutate(method = as.factor(method)) %>% 
  ggplot(aes(x = trans_prob, y = aggregate_nmi, color = method)) + 
  geom_line() +
  theme_minimal() +
  xlab("Transition Prob.") +
  ylab("Aggregate NMI")

ggsave(sprintf("markov_exp_n%d_nla%d_nre%d_ntr%d_seq%d.pdf", n, nlayers, nrep, ntrans, seq_g_update))
  
