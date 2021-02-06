set.seed(1400)  
library(igraph)
library(Matrix)

printf <- function(...) invisible(cat(sprintf(...)))

# set working directory to current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("posterior_analysis.R")
source("hsbm_infer_module.R")
source("network_models.R")

n = 50 # number of nodes in each layer
Tn = 5 # number of layers
eta = matrix(c(0.9,0.75,0.5,0.75,0.60,0.25,0.5,0.25,0.10),3,3)

zb <- list()
for (i in 1:Tn) {
  zb[[i]] <- sample(c(1,2,3),n,replace=TRUE,prob=c(0.45,0.35,0.25))
}

zb <- lapply(zb, sort) 
G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )

#plotG( A[[2]], zb[[2]], deg_func = function(deg) log(deg+3)*2 );

ITRmax <- 100
fname <- paste("Bench_",ITRmax,'.RData',sep='')
RE_ESTIMATE <- T
if (RE_ESTIMATE) {
  dt_hsbm <- system.time( zh_raw <- hsbm_infer(A, beta0=1, gam0=1, ITRmax=ITRmax, Kcap=20, Gcap=20)$zb )
  
  dt_dpsbm <- system.time( zhi_raw <- dpsbm_slice_infer(A, gam0=1, ITRmax=ITRmax, Zcap=20) )
  
  save(A,zb,eta,zh_raw,ITRmax,file=fname)
} else {
  load(fname)  
}

burnin <- ceiling(ITRmax/10)
zh <- compMAPLabel(zh_raw, burnin = burnin, consecutive = T)$MAPlabs
hsbm_nmi <- compute_nmi(zb, zh)
#hsbm_slice_nmi <- compute_slice_nmi(zb,zh)

zhi <- compMAPLabel(zhi_raw, burnin = burnin, consecutive = T)$MAPlabs
dpsbm_nmi <- compute_nmi(zb, zhi)
#dpsbm_slice_nmi <- compute_slice_nmi(zb,zhi)

printf("\n--- Aggregate NMI, eplapsed time (s) ---\n%10s = %1.5f, %3.2f\n%10s = %1.5f, %3.2f\n",
       "HSBM",hsbm_nmi,dt_hsbm["elapsed"],
       "DP-SBM",dpsbm_nmi,dt_dpsbm["elapsed"])

