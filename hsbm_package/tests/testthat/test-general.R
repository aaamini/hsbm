set.seed(1400)
#set.seed(78465)
# set.seed(784)
# library(igraph)
library(Matrix)
library(nett)

# source("posterior_analysis.R")
# source("hsbm_infer_module.R")

n = 200 # number of nodes in each layer
Tn = 5 # number of layers
eta = matrix(c(0.9, 0.75, 0.5,
               0.75, 0.60, 0.25,
               0.5, 0.25, 0.10), 3, 3)

# source("competing_methods.R")

Ktru = 3
zb <- lapply(1:Tn, function(t) sample(c(1,2,3), n, replace=T, prob=c(0.45,0.35,0.25)))
# zb <- lapply(zb, sort)
# G <- lapply(1:Tn, function(t) sample_sbm2(n, eta, zb[[t]]))
# A <- lapply(1:Tn, function(t) as_adj(G[[t]]) )
A <- lapply(1:Tn, function(t) fast_sbm(zb[[t]], eta))

out = fit_mult_dpsbm(A)
out2 = fit_hsbm(A)
str(out$zb_list)
compute_mutual_info(out$zb_list[[50]][[5]], zb[[5]])
get_agg_nmi(out$zb_list[[50]], zb)
get_slice_nmi(out$zb_list[[50]], zb)
get_agg_nmi(out2$zb_list[[50]], zb)

cbind(out$zb_list[[50]][[5]]+1,zb[[5]])

zb_mat = hsbm::fit_dpsbm(A[[1]])
zb_mat2 = hsbm::fit_dpsbm(A[[2]])
lapply(1:niter, function(t) lapply(1:Tn, function(j) zb_mat[,t]))

zb = sample_markov_labels(sample(3, 200, T), 200, 0.2, 3)
out = sample_personality_net(20, 200, 0.75)
out = sample_personality_net(20, 200, 0)
test_markov_labels(out$zb)
microbenchmark::microbenchmark(zz = sample_markov_labels(sample(3, 200, T), 100, .2, 3))


test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})
