library(nett)
library(hsbm)
library(ggplot2)
library(dplyr)
library(patchwork)


n = 200 # number of nodes in each layer
Tn = 3 # number of layers
Ktru = 3
niter = 300

ma_flag = T
if (ma_flag) {
  pri = c(0.45,0.35,0.25)
  # pri = c(1,1,1)
  eta = matrix(c(0.9, 0.75, 0.5, 0.75, 0.60, 0.25, 0.5, 0.25, 0.10), 3, 3)
  # lambda = 60
  # scale = lambda/nett::get_dcsbm_exav_deg(n, pri, eta)
  # eta = eta * scale
  zb <- lapply(1:Tn, function(t) sample(c(1,2,3), n, replace=T, prob=pri))
  zb <- lapply(zb, sort)
} else {
  oir = 0.1
  lambda = 15
  pri = rep(1, Ktru)
  eta = pp_conn(n, oir, lambda, pri=pri)$B
  zb = lapply(1:Tn, function(t) sample(1:Ktru, n, replace=T, prob=pri))
  zb = lapply(zb, sort)
}

A <- lapply(1:Tn, function(t) fast_sbm(zb[[t]], eta))
zb_vec = do.call(rbind, zb)
zb_idx = as.vector(do.call(rbind, lapply(1:Tn, function(j) rep(j, n))))

row_to_list = function(x) lapply(seq_len(nrow(x)), function(i) x[i,])

methods = list()
methods[["hsbm-par"]] = function(A) {
  zall_h = fit_hsbm(A, 
                                   beta0 = .1, gam0 = .5,
                                   alpha_eta = 1, beta_eta = 5,
                                   niter=niter, Kcap=10, Gcap=10, verb = F, 
                                   rand_init = T, seq_g_update = F)$zb_list
  sapply(1:niter, function(t) as.vector(do.call(rbind, zall_h[[t]])))
}

methods[["hsbm-seq"]] = function(A) {
  zall_h = fit_hsbm(A, 
                             beta0 = .1, gam0 = .5,
                             alpha_eta = 1, beta_eta = 5,
                             niter=niter, Kcap=10, Gcap=10, verb = F, 
                             rand_init = T, seq_g_update = T)$zb_list
  sapply(1:niter, function(t) as.vector(do.call(rbind, zall_h[[t]])))
}


mtd_names = names(methods)
nreps = 15
runs = expand.grid(rep = 1:nreps, mi = 1:length(methods))

RNGkind("L'Ecuyer-CMRG")
set.seed(1400)
res = do.call(rbind, parallel::mclapply(1:nrow(runs), function(j) {
  mi = runs[j,"mi"] 
  rep = runs[j, "rep"]
  z_h = methods[[mi]](A)
  
  data.frame(iter = 1:niter, 
             nmi = sapply(1:niter, function(it) compute_mutual_info(z_h[, it], zb_vec)),
             method = mtd_names[mi],
             rep = rep, 
             Kest = sapply(1:niter, function(it) length(unique(z_h[, it])))
  )
}, mc.cores = 4))


res = res %>% mutate(rep = as.factor(rep)) 
p1 = res %>% 
  ggplot(aes(x = iter, y = nmi, color = method, alpha = rep)) + geom_line() +
  theme_minimal(base_size = 12) +
  ggplot2::scale_x_continuous(trans="log10")  +
  scale_alpha_discrete(range = c(0.1, .9)) +
  ylab("Aggregate NMI") + xlab("Iteration")  
  
  

p2 = res %>% mutate(sel_err = abs(Kest-Ktru)) %>% 
  group_by(method, iter) %>% 
  summarize(avg_sel_err = mean(sel_err)) %>%  
  ggplot(aes(x = iter, y = avg_sel_err, color = method)) + geom_line() +
  theme_minimal(base_size = 12)+
  ggplot2::scale_x_continuous(trans="log10") +
  ylab("Average selection error") + xlab("Iteration")


print( knitr::kable(res %>% 
                      filter(iter == niter | iter==2) %>% 
                      arrange(iter), digits = 4, format="pipe") )

print(p1+p2)
ggsave("hsbm_stable.png", width = 10, height=4)


