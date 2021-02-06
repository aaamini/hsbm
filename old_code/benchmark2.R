set.seed(1400)
# set.seed(1234)  
library(igraph)
library(Matrix)

library(nett)

# set working directory to current folder
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("posterior_analysis.R")
source("hsbm_infer_module.R")
source("network_models.R")

n = 100 # number of nodes in each layer
Tn = 5 # number of layers
eta = matrix(c(0.9, 0.75, 0.5, 
               0.75, 0.60, 0.25, 
               0.5, 0.25, 0.10), 3, 3)

Ktru = 3
zb <- list()
for (i in 1:Tn) {
  zb[[i]] <- sample(c(1,2,3), n, replace=T, prob=c(0.45,0.35,0.25))
}

zb <- lapply(zb, sort) 
G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )

ITRmax <- 50
burnin <- ceiling(ITRmax/10)


tau = 0.0
methods = list()
methods[["HSBM"]] = function(A) {
  zh = hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb
  compMAPLabel(zh, burnin = burnin, consecutive = T)$MAPlabs
}
methods[["DP-SBM"]] =  function(A) { 
  zh = dpsbm_slice_infer(A, gam0=1, ITRmax=ITRmax, Zcap=20) 
  compMAPLabel(zh, burnin = burnin, consecutive = T)$MAPlabs
}
methods[["SC"]] = function(A) lapply(1:Tn, function(t) spec_clust(A[[t]], Ktru, tau=tau))
methods[["SC-Abar"]] = function(A) {
  Abar = Reduce(`+`,A)/Tn
  zh1 = spec_clust(Abar, Ktru, tau=tau)
  rep(list(zh1), Tn) # repeat the same label vector for all layers
}

mtd_names = names(methods)

res = do.call(rbind, lapply(1:length(methods), function(j) {
  dt = as.numeric(system.time( zh <- methods[[j]](A) )["elapsed"])
  nmi = compute_nmi(zb, zh)
  data.frame(method_name = mtd_names[j], aggregate_nmi = nmi, elapsed_time = dt)
}))
  
print( knitr::kable(res, digits = 4, format="pipe") )
